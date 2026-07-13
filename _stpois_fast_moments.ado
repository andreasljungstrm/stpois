*! _stpois_fast_moments  version 0.4.0  13jul2026
*! Exact Poisson EHA via iteratively reweighted cell moments
*!
*! Full Newton on (gamma, delta) where gamma are continuous-covariate
*! coefficients and delta are cell-level (categorical) coefficients.
*! Each iteration makes one pass over the microdata to accumulate
*! exponentially tilted within-cell moments (sum mu, sum mu*X, sum mu*XX');
*! all remaining algebra is at the cell level. Converges to the exact
*! individual-level MLE with exact OIM/robust/cluster standard errors.
*! Cost per iteration is O(n p^2) with p = number of continuous covariates,
*! versus O(n (p+k)^2) for poisson with k categorical parameters.
*!
*! Called by stpois when fast(moments) is specified. Not intended for direct use.
program _stpois_fast_moments, eclass
    version 14

    syntax varlist(fv ts) [aw iw fw pw], ///
        touse(varname)                    ///
        exposure(varname)                 ///
        [tol(real 1e-8)                   ///
        maxiter(integer 100)              ///
        nolog                             ///
        noconstant                        ///
        vce(string)                       ///
        robust                            ///
        cluster(varname)                  ///
        SKEWness                          ///
        *]

    if "`skewness'" != "" {
        di as txt "  (note: skewness is obsolete — fast(moments) is now exact; option ignored)"
    }

    // ── Parse varlist: categorical (i.) vs continuous (bare) ───────────────
    local cont_vars ""
    local cat_vars_raw ""
    local cat_terms ""

    foreach tok of local varlist {
        if strpos("`tok'", "#") > 0 {
            di as error "Interactions not supported with fast(moments)"
            exit 198
        }
        if regexm("`tok'", "^[0-9]*b?i\.") | substr("`tok'", 1, 2) == "i." {
            local rawname = regexr("`tok'", "^[0-9]*b?i\.", "")
            local cat_vars_raw `cat_vars_raw' `rawname'
            local cat_terms    `cat_terms' `tok'
        }
        else {
            if !regexm("`tok'", "\.") local rawname `tok'
            else local rawname = regexr("`tok'", "^[a-zA-Z]+\.", "")
            if "`rawname'" == "" local rawname `tok'
            local cont_vars `cont_vars' `rawname'
        }
    }

    if "`cont_vars'" == "" {
        di as error "fast(moments) requires at least one continuous covariate (bare variable name)"
        exit 198
    }
    if "`cat_vars_raw'" == "" {
        di as error "fast(moments) requires at least one categorical covariate (i.varname)"
        exit 198
    }

    local p : word count `cont_vars'
    local hascons = cond("`constant'" == "noconstant", 0, 1)

    // ── Parse VCE type ─────────────────────────────────────────────────────
    local ctype    "oim"
    local clusterv ""
    if "`robust'" != "" local ctype "robust"
    if `"`vce'"' != "" {
        local vcetype = lower(word("`vce'", 1))
        if inlist("`vcetype'", "r", "robust") {
            local ctype "robust"
        }
        else if inlist("`vcetype'", "cl", "clu", "clus", "clust", "cluste", "cluster") {
            local ctype    "cluster"
            local clusterv = word("`vce'", 2)
        }
        else if !inlist("`vcetype'", "oim", "") {
            di as error "fast(moments) supports vce(oim), vce(robust), vce(cluster varname)"
            exit 198
        }
    }
    if "`cluster'" != "" {
        local ctype    "cluster"
        local clusterv "`cluster'"
    }
    if "`ctype'" == "cluster" & "`clusterv'" == "" {
        di as error "cluster VCE requires a cluster variable"
        exit 198
    }
    if "`clusterv'" != "" {
        capture confirm numeric variable `clusterv'
        if _rc {
            di as error "cluster variable must be numeric with fast(moments)"
            exit 198
        }
    }

    // ── Count obs before preserve ──────────────────────────────────────────
    qui count if `touse'
    local N_orig = r(N)

    // ── Load microdata into Mata, then collapse for the cell design ────────
    preserve
    qui keep if `touse'

    tempvar cellid ptime
    qui egen long `cellid' = group(`cat_vars_raw')
    qui gen double `ptime' = `exposure'
    sort `cellid'

    if "`nolog'" == "" di as txt "  (fast moments) loading microdata..."
    mata: _stpm_load("_d", "`cont_vars'", "`ptime'", "`cellid'", "`clusterv'")

    // Collapse to cells (D_j and the categorical values); sorted by cellid
    collapse (sum) _d, by(`cellid' `cat_vars_raw')
    sort `cellid'
    local N_cells = _N

    // Cell-level categorical design: explicit dummies for non-base levels
    fvexpand `cat_terms'
    local xvars_fv_all `r(varlist)'
    local xvars_fv ""
    local wcols ""
    foreach v of local xvars_fv_all {
        _ms_parse_parts `v'
        local _omit = cond("`r(omit)'" == "", "0", "`r(omit)'")
        if "`_omit'" != "1" {
            local xvars_fv `xvars_fv' `v'
            tempvar dum
            qui gen byte `dum' = (`r(name)' == `r(level)')
            local wcols `wcols' `dum'
        }
    }
    local k : word count `wcols'

    // ── Newton iterations in Mata ──────────────────────────────────────────
    if "`nolog'" == "" di as txt "  (fast moments) Newton on `p' continuous + `k' cell parameters, `N_cells' cells..."

    tempname b_mat V_mat
    scalar _stpm_ll    = .
    scalar _stpm_ll0   = .
    scalar _stpm_iter  = .
    scalar _stpm_conv  = .
    scalar _stpm_nclu  = .
    scalar _stpm_rank  = .

    mata: _stpm_newton(         ///
        "`wcols'",              ///
        "_d",                   ///
        `hascons',              ///
        `tol',                  ///
        `maxiter',              ///
        "`ctype'",              ///
        "`b_mat'",              ///
        "`V_mat'",              ///
        "_stpm_ll",             ///
        "_stpm_ll0",            ///
        "_stpm_iter",           ///
        "_stpm_conv",           ///
        "_stpm_nclu",           ///
        "_stpm_rank"            ///
    )

    local ll      = scalar(_stpm_ll)
    local ll_0    = scalar(_stpm_ll0)
    local iter    = scalar(_stpm_iter)
    local conv    = scalar(_stpm_conv)
    local n_clust = scalar(_stpm_nclu)
    local rank    = scalar(_stpm_rank)
    scalar drop _stpm_ll _stpm_ll0 _stpm_iter _stpm_conv _stpm_nclu _stpm_rank

    restore

    // ── Names, chi2, post ──────────────────────────────────────────────────
    local bnames `cont_vars' `xvars_fv'
    if `hascons' local bnames `bnames' _cons
    matrix colnames `b_mat' = `bnames'
    matrix colnames `V_mat' = `bnames'
    matrix rownames `V_mat' = `bnames'

    local ncols = colsof(`b_mat')
    local blankeq
    forvalues i = 1/`ncols' {
        local blankeq `blankeq' ""
    }
    matrix coleq `b_mat' = `blankeq'
    matrix coleq `V_mat' = `blankeq'
    matrix roweq `V_mat' = `blankeq'

    local df_m = `rank' - `hascons'
    local chi2 = 2 * (`ll' - `ll_0')

    ereturn clear
    ereturn post `b_mat' `V_mat', esample(`touse') obs(`N_orig')

    ereturn scalar ll        = `ll'
    ereturn scalar ll_0      = `ll_0'
    ereturn scalar chi2      = `chi2'
    ereturn scalar df_m      = `df_m'
    ereturn scalar rank      = `rank'
    ereturn scalar ic        = `iter'
    ereturn scalar converged = `conv'
    ereturn scalar N_cells   = `N_cells'
    ereturn local  vce       "`ctype'"
    ereturn local  cont_vars "`cont_vars'"
    ereturn local  cat_vars  "`cat_vars_raw'"
    if "`ctype'" == "cluster" {
        ereturn scalar N_clust  = `n_clust'
        ereturn local  clustvar "`clusterv'"
    }

    if !`conv' {
        di as txt "  Warning: Newton did not converge in `maxiter' iterations"
    }
    if "`nolog'" == "" {
        di as txt "  (fast moments) `N_orig' obs, `N_cells' cells, `iter' Newton iterations, converged=" `conv'
        if "`ctype'" == "cluster" {
            di as txt "  (Std. err. adjusted for " `n_clust' " clusters in `clusterv')"
        }
    }
end


// ============================================================================
// Mata
// ============================================================================
mata:

// Microdata held across the collapse via externals
void _stpm_load(string scalar dvar, string scalar xstr,
                string scalar tvar, string scalar cellvar,
                string scalar cluvar)
{
    external real colvector STPM_d, STPM_t, STPM_cell, STPM_clu
    external real matrix    STPM_X

    STPM_d    = st_data(., dvar)
    STPM_t    = st_data(., tvar)
    STPM_cell = st_data(., cellvar)
    STPM_X    = st_data(., tokens(xstr))
    if (cluvar != "") STPM_clu = st_data(., cluvar)
    else              STPM_clu = J(0, 1, .)
}


void _stpm_newton(string scalar wstr, string scalar Dvar,
    real scalar hascons, real scalar tol, real scalar maxiter,
    string scalar ctype,
    string scalar b_ret, string scalar V_ret,
    string scalar ll_ret, string scalar ll0_ret,
    string scalar iter_ret, string scalar conv_ret,
    string scalar nclu_ret, string scalar rank_ret)
{
    external real colvector STPM_d, STPM_t, STPM_cell, STPM_clu
    external real matrix    STPM_X

    real matrix    W, info, Xj, H, Hinv, V, S, B, MX
    real colvector Dj, d, t, cell, mu, eta, r, theta, step, theta_new
    real colvector Mj, Rj, s, etac, ucids, si
    real scalar    n, p, k, q, J_, j, iter, ll, ll_new, ll0, a0
    real scalar    converged, halv, delta_c, sum_d, sum_t, rank
    real scalar    G_clust, i_c
    real colvector cids

    d    = STPM_d
    t    = STPM_t
    cell = STPM_cell
    n    = rows(d)
    p    = cols(STPM_X)

    // Cell design from the collapsed data currently in memory
    W  = st_data(., tokens(wstr))
    Dj = st_data(., Dvar)
    if (hascons) W = W, J(rows(W), 1, 1)
    J_ = rows(W)
    k  = cols(W)
    q  = p + k

    // Microdata are sorted by cell; panel info once
    info = panelsetup(cell, 1)

    // Fixed data-side score component
    // S_dx = sum_i d_i X_i ; cell sums D_j already in Dj

    // Null model (exposure + constant only): closed form
    sum_d = quadsum(d)
    sum_t = quadsum(t)
    a0    = ln(sum_d / sum_t)
    ll0   = quadsum(d :* (ln(t :+ 1e-300) :+ a0)) - sum_d

    // Initialize: gamma = 0, delta = 0 except constant
    theta = J(q, 1, 0)
    if (hascons) theta[q] = a0

    converged = 0
    ll = .

    for (iter = 1; iter <= maxiter; iter++) {
        // eta_i = X_i gamma + w_{j(i)} delta
        etac = W * theta[|p+1 \ q|]
        eta  = STPM_X * theta[|1 \ p|] + etac[cell]
        mu   = t :* exp(eta)
        mu   = mu + (mu :< 1e-300) :* 1e-300
        ll   = quadsum(d :* ln(mu) - mu)
        r    = d - mu

        // One pass: cell sums of mu and mu*X
        Mj = J(J_, 1, 0)
        MX = J(J_, p, 0)
        for (j = 1; j <= J_; j++) {
            Mj[j]    = quadsum(panelsubmatrix(mu, j, info))
            MX[j, .] = quadcolsum(panelsubmatrix(mu, j, info) :*
                                  panelsubmatrix(STPM_X, j, info))
        }
        Rj = Dj - Mj

        // Score
        s = quadcross(STPM_X, r) \ quadcross(W, Rj)

        // Hessian (negative of): blocks gg, gd, dd
        H = quadcross(STPM_X, mu, STPM_X), (MX' * W) \
            (MX' * W)', quadcross(W, Mj, W)

        Hinv = invsym(H)
        step = Hinv * s

        // Step-halving on the log likelihood
        theta_new = theta + step
        for (halv = 1; halv <= 30; halv++) {
            etac   = W * theta_new[|p+1 \ q|]
            eta    = STPM_X * theta_new[|1 \ p|] + etac[cell]
            mu     = t :* exp(eta)
            mu     = mu + (mu :< 1e-300) :* 1e-300
            ll_new = quadsum(d :* ln(mu) - mu)
            if (ll_new >= ll - 1e-12) break
            step      = step / 2
            theta_new = theta + step
        }

        delta_c = max(abs(theta_new - theta) :/ (1 :+ abs(theta_new)))
        theta   = theta_new
        ll      = ll_new
        if (delta_c < tol) {
            converged = 1
            break
        }
    }

    // Final Hessian at the optimum for the VCE
    etac = W * theta[|p+1 \ q|]
    eta  = STPM_X * theta[|1 \ p|] + etac[cell]
    mu   = t :* exp(eta)
    mu   = mu + (mu :< 1e-300) :* 1e-300
    r    = d - mu

    Mj = J(J_, 1, 0)
    MX = J(J_, p, 0)
    for (j = 1; j <= J_; j++) {
        Mj[j]    = quadsum(panelsubmatrix(mu, j, info))
        MX[j, .] = quadcolsum(panelsubmatrix(mu, j, info) :*
                              panelsubmatrix(STPM_X, j, info))
    }
    H = quadcross(STPM_X, mu, STPM_X), (MX' * W) \
        (MX' * W)', quadcross(W, Mj, W)
    Hinv = invsym(H)
    rank = sum(diagonal(Hinv) :!= 0)

    G_clust = 0
    if (ctype == "oim") {
        V = Hinv
    }
    else {
        // Per-observation scores s_i = r_i * (X_i, w_{j(i)})
        S = r :* (STPM_X, W[cell, .])
        if (ctype == "robust") {
            B = quadcross(S, S) * (n / (n - 1))
            V = Hinv * B * Hinv
        }
        else {
            cids    = STPM_clu
            ucids   = uniqrows(cids)
            G_clust = rows(ucids)
            B       = J(q, q, 0)
            for (i_c = 1; i_c <= G_clust; i_c++) {
                si = quadcolsum(S[selectindex(cids :== ucids[i_c, 1]), .])
                B  = B + si' * si
            }
            B = B * (G_clust / (G_clust - 1))
            V = Hinv * B * Hinv
        }
    }

    st_matrix(b_ret, theta')
    st_matrix(V_ret, V)
    st_numscalar(ll_ret,   ll)
    st_numscalar(ll0_ret,  ll0)
    st_numscalar(iter_ret, iter)
    st_numscalar(conv_ret, converged)
    st_numscalar(nclu_ret, G_clust)
    st_numscalar(rank_ret, rank)

    // Release microdata
    STPM_d = STPM_t = STPM_cell = STPM_clu = J(0, 1, .)
    STPM_X = J(0, 0, .)
}

end
