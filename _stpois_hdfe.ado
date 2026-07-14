*! _stpois_hdfe  version 0.6.0  14jul2026
*! IRLS Poisson with alternating-projections FE absorption (ppmlhdfe-style)
*! Called by stpois when absorb() is specified. Not intended for direct use.
program _stpois_hdfe, eclass
    version 14

    syntax varlist(fv ts) [aw iw fw pw], ///
        touse(varname)                    ///
        exposure(varname)                 ///
        absorb(string)                    ///
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

    foreach aterm of local absorb {
        local araw = subinstr("`aterm'", "#", " ", .)
        local araw = subinstr("`araw'", "i.", "", .)
        tempvar n_ev
        qui bysort `araw': egen `n_ev' = total(_d * `touse2')
        qui count if `touse2' & `n_ev' == 0
        if r(N) > 0 {
            di as txt "  (absorb(`aterm'): dropped " r(N) " obs in zero-event groups)"
            qui replace `touse2' = 0 if `n_ev' == 0
        }
        drop `n_ev'
    }

    // --- Expand factor variables ----------------------------------------------
    // fvexpand and fvrevar return parallel lists (one entry per expanded
    // term, including interactions); keep the non-omitted pairs.
    fvexpand `varlist' if `touse2'
    local xvars_fv_all `r(varlist)'
    fvrevar `varlist' if `touse2'
    local xcols_all `r(varlist)'

    local n_terms : word count `xvars_fv_all'
    local n_cols  : word count `xcols_all'
    if `n_terms' != `n_cols' {
        di as error "absorb(): could not expand factor-variable terms"
        exit 198
    }

    local xvars_fv ""
    local xvars_expanded ""
    forvalues ti = 1/`n_terms' {
        local trm : word `ti' of `xvars_fv_all'
        local cv  : word `ti' of `xcols_all'
        _ms_parse_parts `trm'
        local _omit = cond("`r(omit)'" == "", "0", "`r(omit)'")
        if "`_omit'" != "1" {
            local xvars_fv       `xvars_fv' `trm'
            local xvars_expanded `xvars_expanded' `cv'
        }
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

    // FE spec for Mata: terms separated by ";", vars within a term by space
    local fespec ""
    foreach aterm of local absorb {
        local araw = subinstr("`aterm'", "#", " ", .)
        local araw = subinstr("`araw'", "i.", "", .)
        local fespec "`fespec';`araw'"
    }
    local fespec = substr("`fespec'", 2, .)

    mata: _stpois_hdfe_irls(   ///
        "_d",                  ///
        "`xvars_expanded'",    ///
        "`fespec'",            ///
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

// Panel boundaries of a sorted key matrix: rows where any key column
// changes start a new panel. Returns (start, end) rows, the same
// format as panelsetup().
real matrix _sthdfe_panels(real matrix Ks)
{
    real colvector flag, starts, ends
    real scalar    n

    n    = rows(Ks)
    flag = J(n, 1, 1)
    if (n > 1) {
        flag[|2 \ n|] = rowmax(Ks[|2, 1 \ n, .|] :!= Ks[|1, 1 \ n - 1, .|])
    }
    starts = selectindex(flag)
    if (rows(starts) > 1) ends = starts[|2 \ rows(starts)|] :- 1 \ n
    else                  ends = J(1, 1, n)
    return((starts, ends))
}

// Segment sums over a sorted matrix via running sums: O(n) regardless
// of the number of segments.
real matrix _sthdfe_segsum(real matrix Ap, real colvector ends)
{
    real matrix    out
    real colvector cs
    real scalar    j, nseg

    nseg = rows(ends)
    out  = J(nseg, cols(Ap), .)
    for (j = 1; j <= cols(Ap); j++) {
        cs         = quadrunningsum(Ap[., j])
        out[., j]  = cs[ends]
        if (nseg > 1) {
            out[|2, j \ nseg, j|] = out[|2, j \ nseg, j|] -
                                    cs[ends[|1 \ nseg - 1|]]
        }
    }
    return(out)
}

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
    string rowvector xnames, fes
    real scalar    p, nfe, n
    real matrix    X, M, P
    real colvector prm
    pointer(real matrix) rowvector infos
    real scalar    sum_y, sum_E
    real colvector eta, b_old, b_new
    real colvector mu, z, w, z_tilde, res, res_tilde, FE_c, z_prev, r_prev
    real matrix    X_tilde, XtW, XtWX, Xt_prev, H, Hinv, V, S, B, SG, info_c
    real colvector XtWz, lmu, cids, permc
    real scalar    iter, inner, delta, converged, j, k, ll, G_clust

    // --- Load data -----------------------------------------------------------
    touse_col = st_data(., tousevar)
    idx       = selectindex(touse_col)
    n         = rows(idx)

    y  = st_data(idx, yvar)
    E  = st_data(idx, evar)

    xnames = tokens(xvars_str)
    p      = cols(xnames)
    X      = st_data(idx, xnames)

    // FE spec: terms separated by ";"; each term lists the variables whose
    // cross-classification is absorbed (single variable or interaction)
    fes = tokens(fevars_str, ";")
    fes = select(fes, fes :!= ";")
    nfe = cols(fes)

    // Precompute, per FE, the sort permutation and panel boundaries once;
    // every demeaning pass is then O(n) regardless of the level count
    P     = J(n, nfe, .)
    infos = J(1, nfe, NULL)
    for (k = 1; k <= nfe; k++) {
        M        = st_data(idx, tokens(fes[k]))
        prm      = order(M, (1..cols(M)))
        P[., k]  = prm
        infos[k] = &_sthdfe_panels(M[prm, .])
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
        // Cluster sandwich: scores summed within clusters via one sort
        // and segment sums, G/(G-1) correction
        cids    = st_data(idx, clustervar)
        permc   = order(cids, 1)
        info_c  = _sthdfe_panels(cids[permc])
        G_clust = rows(info_c)
        SG      = _sthdfe_segsum(S[permc, .], info_c[., 2])
        B       = quadcross(SG, SG) * (G_clust / (G_clust - 1))
        V       = Hinv * B * Hinv
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
