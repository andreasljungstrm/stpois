*! _stpois_hdfe  version 0.5.0  13jul2026
*! IRLS Poisson with alternating-projections FE absorption (ppmlhdfe-style)
*! Called by stpois when absorb() is specified. Not intended for direct use.
program _stpois_hdfe, eclass
    version 14

    syntax varlist(fv ts) [aw iw fw pw], ///
        touse(varname)                    ///
        exposure(varname)                 ///
        absorb(varlist)                   ///
        [tol(real 1e-8)                   ///
        maxiter(integer 100)              ///
        nolog                             ///
        noconstant                        ///
        vce(string)                       ///
        robust                            ///
        cluster(varname)                  ///
        *]

    // --- Parse VCE type and cluster variable ---------------------------------
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
            di as error "absorb() supports vce(oim), vce(robust), vce(cluster varname)"
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

    // --- Separation check ----------------------------------------------------
    tempvar touse2
    qui gen byte `touse2' = `touse'
    if "`clusterv'" != "" markout `touse2' `clusterv'

    foreach fevar of local absorb {
        tempvar n_ev
        qui bysort `fevar': egen `n_ev' = total(_d * `touse2')
        qui count if `touse2' & `n_ev' == 0
        if r(N) > 0 {
            di as txt "  (absorb(`fevar'): dropped " r(N) " obs in zero-event groups)"
            qui replace `touse2' = 0 if `n_ev' == 0
        }
        drop `n_ev'
    }

    // --- Expand factor variables ----------------------------------------------
    fvrevar `varlist' if `touse2'
    local xvars_expanded `r(varlist)'

    fvexpand `varlist' if `touse2'
    local xvars_fv_all `r(varlist)'
    local xvars_fv ""
    foreach v of local xvars_fv_all {
        _ms_parse_parts `v'
        local _omit = cond("`r(omit)'" == "", "0", "`r(omit)'")
        if "`_omit'" != "1" local xvars_fv `xvars_fv' `v'
    }

    local p : word count `xvars_expanded'

    // --- Count clusters if needed --------------------------------------------
    local n_clust 0
    if "`ctype'" == "cluster" {
        // count clusters in Mata; levelsof overflows the macro buffer
        // with very many levels
        mata: st_local("n_clust", ///
            strofreal(rows(uniqrows(st_data(., "`clusterv'", "`touse2'")))))
    }

    // --- Run Mata IRLS -------------------------------------------------------
    tempname b_mat V_mat
    qui count if `touse2'
    local N = r(N)

    scalar _hdfe_ll   = .
    scalar _hdfe_iter = .
    scalar _hdfe_conv = .

    mata: _stpois_hdfe_irls(   ///
        "_d",                  ///
        "`xvars_expanded'",    ///
        "`absorb'",            ///
        "`exposure'",          ///
        "`touse2'",            ///
        `tol',                 ///
        `maxiter',             ///
        "`ctype'",             ///
        "`clusterv'",          ///
        `n_clust',             ///
        "`b_mat'",             ///
        "`V_mat'",             ///
        "_hdfe_ll",            ///
        "_hdfe_iter",          ///
        "_hdfe_conv"           ///
    )

    local ll   = scalar(_hdfe_ll)
    local iter = scalar(_hdfe_iter)
    local conv = scalar(_hdfe_conv)
    scalar drop _hdfe_ll _hdfe_iter _hdfe_conv

    // --- Fix colnames on b and V ---------------------------------------------
    matrix colnames `b_mat' = `xvars_fv'
    matrix colnames `V_mat' = `xvars_fv'
    matrix rownames `V_mat' = `xvars_fv'

    local ncols = colsof(`b_mat')
    local blankeq
    forvalues i = 1/`ncols' {
        local blankeq `blankeq' ""
    }
    matrix coleq `b_mat' = `blankeq'
    matrix coleq `V_mat' = `blankeq'
    matrix roweq `V_mat' = `blankeq'

    // --- Null LL for LR chi2 -------------------------------------------------
    qui poisson _d if `touse2', exposure(`exposure') nolog
    local ll_0 = e(ll)
    local chi2 = 2 * (`ll' - `ll_0')

    // --- Post results --------------------------------------------------------
    ereturn clear
    ereturn post `b_mat' `V_mat', esample(`touse2') obs(`N')

    ereturn scalar ll        = `ll'
    ereturn scalar ll_0      = `ll_0'
    ereturn scalar chi2      = `chi2'
    ereturn scalar df_m      = `p'
    ereturn scalar rank      = `p'
    ereturn scalar ic        = `iter'
    ereturn scalar converged = `conv'
    ereturn local  vce       "`ctype'"
    if "`ctype'" == "cluster" {
        ereturn scalar N_clust  = `n_clust'
        ereturn local  clustvar "`clusterv'"
    }

    if !`conv' {
        di as txt "  Warning: HDFE IRLS did not converge in `maxiter' iterations"
    }
    if "`nolog'" == "" {
        di as txt "  (HDFE: `iter' IRLS iterations, converged=" `conv' ")"
        if "`ctype'" == "cluster" {
            di as txt "  (Std. err. adjusted for `n_clust' clusters in `clusterv')"
        }
    }
end


// ============================================================================
// Mata
// ============================================================================
mata:

// Weighted within-group demeaning using a precomputed sort permutation
// and panel boundaries: one O(n) pass per call, independent of the
// number of group levels.
real colvector _wdemean(real colvector v,
                        real colvector w,
                        real colvector perm,
                        real matrix    info)
{
    real colvector vp, wp, res, out
    real scalar    i, r1, r2, wsum

    vp  = v[perm]
    wp  = w[perm]
    res = vp

    for (i = 1; i <= rows(info); i++) {
        r1   = info[i, 1]
        r2   = info[i, 2]
        wsum = quadsum(wp[|r1 \ r2|])
        if (wsum > 0) {
            res[|r1 \ r2|] = vp[|r1 \ r2|] :-
                (quadsum(wp[|r1 \ r2|] :* vp[|r1 \ r2|]) / wsum)
        }
    }
    out       = v
    out[perm] = res
    return(out)
}


void _stpois_hdfe_irls(
    string scalar yvar,
    string scalar xvars_str,
    string scalar fevars_str,
    string scalar evar,
    string scalar tousevar,
    real   scalar tol,
    real   scalar maxiter,
    string scalar ctype,
    string scalar clustervar,
    real   scalar n_clust,
    string scalar b_ret,
    string scalar V_ret,
    string scalar ll_ret,
    string scalar iter_ret,
    string scalar conv_ret)
{
    // Declare all variables at function scope
    real colvector touse_col, idx
    real colvector y, E
    string rowvector xnames, fenames
    real scalar    p, nfe, n
    real matrix    X, G, P
    pointer(real matrix) rowvector infos
    real scalar    sum_y, sum_E
    real colvector eta, b_old, b_new
    real colvector mu, z, w, z_tilde, res, res_tilde, FE_c, z_prev, r_prev
    real matrix    X_tilde, XtW, XtWX, Xt_prev, H, Hinv, V, S, B
    real colvector XtWz, lmu, cids, ucids, si
    real scalar    iter, inner, delta, converged, j, k, ll, G_clust, i_c

    // --- Load data -----------------------------------------------------------
    touse_col = st_data(., tousevar)
    idx       = selectindex(touse_col)
    n         = rows(idx)

    y  = st_data(idx, yvar)
    E  = st_data(idx, evar)

    xnames  = tokens(xvars_str)
    fenames = tokens(fevars_str)
    p   = cols(xnames)
    nfe = cols(fenames)

    X = st_data(idx, xnames)
    G = st_data(idx, fenames)

    // Precompute, per FE, the sort permutation and panel boundaries once;
    // every demeaning pass is then O(n) regardless of the level count
    P     = J(n, nfe, .)
    infos = J(1, nfe, NULL)
    for (k = 1; k <= nfe; k++) {
        P[., k]  = order(G[., k], 1)
        infos[k] = &panelsetup(G[P[., k], k], 1)
    }

    // --- Initialize ----------------------------------------------------------
    sum_y = quadsum(y)
    sum_E = quadsum(E)
    eta   = J(n, 1, log(sum_y / sum_E))
    b_old = J(p, 1, 0)
    b_new = J(p, 1, 0)

    // --- IRLS loop -----------------------------------------------------------
    converged = 0

    for (iter = 1; iter <= maxiter; iter++) {
        mu = exp(eta) :* E
        mu = mu + (mu :< 1e-300) :* 1e-300
        z  = eta + (y - mu) :/ mu
        w  = mu

        z_tilde = z
        X_tilde = X

        // Inner alternating-projections (demean z and X by each FE)
        for (inner = 1; inner <= 200; inner++) {
            z_prev = z_tilde
            for (k = 1; k <= nfe; k++) {
                z_tilde = _wdemean(z_tilde, w, P[., k], *infos[k])
                for (j = 1; j <= p; j++) {
                    X_tilde[., j] = _wdemean(X_tilde[., j], w, P[., k], *infos[k])
                }
            }
            if (max(abs(z_tilde - z_prev)) < tol * 0.01) break
        }

        // WLS: b = (X_tilde'W X_tilde)^{-1} X_tilde'W z_tilde
        XtW   = (X_tilde :* w)'
        XtWX  = XtW * X_tilde
        XtWz  = XtW * z_tilde
        b_new = invsym(XtWX) * XtWz

        // Recover FE contributions: res - M_W(res)
        res       = z - X * b_new
        res_tilde = res
        for (inner = 1; inner <= 200; inner++) {
            r_prev = res_tilde
            for (k = 1; k <= nfe; k++) {
                res_tilde = _wdemean(res_tilde, w, P[., k], *infos[k])
            }
            if (max(abs(res_tilde - r_prev)) < tol * 0.01) break
        }
        FE_c = res - res_tilde

        eta   = X * b_new + FE_c
        delta = max(abs(b_new - b_old) :/ (1 :+ abs(b_new)))
        b_old = b_new
        if (delta < tol) {
            converged = 1
            break
        }
    }

    // --- VCE: demean X at final weights -------------------------------------
    mu      = exp(eta) :* E
    mu      = mu + (mu :< 1e-300) :* 1e-300
    w       = mu
    X_tilde = X

    for (inner = 1; inner <= 200; inner++) {
        Xt_prev = X_tilde
        for (k = 1; k <= nfe; k++) {
            for (j = 1; j <= p; j++) {
                X_tilde[., j] = _wdemean(X_tilde[., j], w, P[., k], *infos[k])
            }
        }
        if (max(abs(X_tilde - Xt_prev)) < tol * 0.01) break
    }

    H    = quadcross(X_tilde, w, X_tilde)
    Hinv = invsym(H)

    // Scores: s_i = (y_i - mu_i) * X_tilde_i  (n x p matrix)
    S = (y - mu) :* X_tilde

    if (ctype == "oim") {
        // OIM: V = H^{-1}  (exact)
        V = Hinv
    }
    else if (ctype == "robust") {
        // HC1: multiply meat by n/(n-1) to match Stata's robust
        B = quadcross(S, S) * (n / (n - 1))
        V = Hinv * B * Hinv
    }
    else {
        // Cluster sandwich: sum scores within cluster, G/(G-1) correction
        cids    = st_data(idx, clustervar)
        ucids   = uniqrows(cids)
        G_clust = rows(ucids)
        B       = J(p, p, 0)
        for (i_c = 1; i_c <= G_clust; i_c++) {
            si = quadcolsum(S[selectindex(cids :== ucids[i_c, 1]), .])
            B  = B + si' * si
        }
        B = B * (G_clust / (G_clust - 1))
        V = Hinv * B * Hinv
    }

    // --- Log-likelihood ------------------------------------------------------
    lmu = ln(mu)
    ll  = quadsum(y :* lmu - mu)

    // --- Return --------------------------------------------------------------
    st_matrix(b_ret, b_new')
    st_matrix(V_ret, V)
    st_numscalar(ll_ret,   ll)
    st_numscalar(iter_ret, iter)
    st_numscalar(conv_ret, converged)
}

end
