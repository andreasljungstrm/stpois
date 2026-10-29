*! _stpois_fast  version 0.8.1  15jul2026
*! Exact Poisson EHA accelerated through the cell structure of the
*! categorical covariates, with optional absorbed fixed effects and
*! weights.
*!
*! Terms in the varlist are classified as cell-level (built only from
*! factor-variable pieces, including factor#factor interactions) or
*! individual-level (involving at least one continuous piece, including
*! continuous#continuous and continuous#factor interactions). Newton
*! iterations compute the score and Hessian from exponentially tilted
*! within-cell moments of the individual-level columns: one O(n p^2)
*! pass over the microdata per iteration, all remaining algebra on the
*! J cells. When every term is cell-level, nothing is absorbed, and the
*! default observed-information VCE is requested, estimation is routed
*! to Stata's compiled collapse + poisson on the (events, exposure)
*! cell sums -- exact by sufficiency and faster than iterating in Mata.
*! Robust/cluster VCEs keep the Mata cell-sum Newton, since they need
*! the per-observation microdata scores.
*!
*! With fespec() (stpois option absorb()), the fixed effects are
*! concentrated out of the likelihood by closed-form Gauss-Seidel
*! updates -- for Poisson, the profile-maximizing effect of group g is
*! alpha_g = ln(sum_g w*d / sum_g w*mu) -- iterated to convergence
*! inside every Newton step (the approach of R's fixest). The VCE then
*! comes from the Frisch-Waugh-Lovell partitioned information: all
*! design columns are weighted-demeaned within the absorbed groups and
*! the sandwich is built on the demeaned columns, exactly as in the
*! ppmlhdfe approach.
*!
*! Weights (fw/iw/pw) multiply the likelihood contributions; pweights
*! are given robust standard errors by the caller.
*!
*! Estimates are the exact (weighted) individual-level MLE:
*! coefficients, log likelihood, and OIM/robust/cluster VCEs match
*! poisson numerically.
*!
*! Called by stpois when fast is specified. Not intended for direct use.
program _stpois_fast, eclass
    version 14

    syntax varlist(fv ts), ///
        touse(varname)                    ///
        exposure(varname)                 ///
        [tol(real 1e-8)                   ///
        maxiter(integer 100)              ///
        wvar(varname)                     ///
        wtype(string)                     ///
        fespec(string)                    ///
        dvar(string)                      ///
        nolog                             ///
        noconstant                        ///
        vce(string)                       ///
        robust                            ///
        cluster(varname)                  ///
        *]

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
            di as error "fast supports vce(oim), vce(robust), vce(cluster varname)"
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
            di as error "cluster variable must be numeric with fast"
            exit 198
        }
    }

    // Absorbed fixed effects soak up the constant
    local hasfe = (`"`fespec'"' != "")
    local hascons = cond("`constant'" == "noconstant" | `hasfe', 0, 1)

    // ── Classify expanded terms: cell-level vs individual-level ────────────
    // Each expanded term is a product of pieces (split on #). A piece that
    // is <level>.<varname> is categorical; c.<varname> or a bare name is
    // continuous. Terms whose pieces are all categorical are constant
    // within cells and enter the cell design W; the rest enter X.
    fvexpand `varlist' if `touse'
    local terms_all `r(varlist)'

    local contnames ""    // unique continuous variables (X source columns)
    local catnames ""     // unique categorical variables (C source columns)
    local xspec ""        // piece-encoded X terms
    local wspec ""        // piece-encoded W terms
    local xterms ""       // display names for X columns
    local wterms ""       // display names for W columns
    local wcat_idx ""     // C-column indices that define the cells

    foreach trm of local terms_all {
        _ms_parse_parts `trm'
        local _omit = cond("`r(omit)'" == "", "0", "`r(omit)'")
        if "`_omit'" == "1" continue

        local pieces = subinstr("`trm'", "#", " ", .)
        local spec ""
        local iscat 1
        local thiscat ""
        foreach pc of local pieces {
            if regexm("`pc'", "^([0-9]+)[bo]*\.(.+)$") {
                local lvl = regexs(1)
                local vn  = regexs(2)
                local k : list posof "`vn'" in catnames
                if `k' == 0 {
                    local catnames `catnames' `vn'
                    local k : list sizeof catnames
                }
                local spec `spec'*d`k'_`lvl'
                local thiscat `thiscat' `k'
            }
            else {
                local vn "`pc'"
                if substr("`pc'", 1, 2) == "c." local vn = substr("`pc'", 3, .)
                if strpos("`vn'", ".") > 0 {
                    di as error "operator `pc' not supported with fast"
                    exit 198
                }
                local k : list posof "`vn'" in contnames
                if `k' == 0 {
                    local contnames `contnames' `vn'
                    local k : list sizeof contnames
                }
                local spec `spec'*c`k'
                local iscat 0
            }
        }
        local spec = substr("`spec'", 2, .)
        if `iscat' {
            local wspec  `wspec' `spec'
            local wterms `wterms' `trm'
            local wcat_idx : list wcat_idx | thiscat
        }
        else {
            local xspec  `xspec' `spec'
            local xterms `xterms' `trm'
        }
    }

    if "`wterms'" == "" {
        di as error "fast requires at least one categorical term (i.varname)"
        di as error "categorical terms define the cells over which computation is accelerated"
        if `hasfe' di as error "(for models with only absorbed fixed effects, use absorb() without fast)"
        exit 198
    }

    local p : word count `xterms'
    local k : word count `wterms'

    // ── Count obs (Stata convention: sum of fweights, else # of obs) ──────
    qui count if `touse'
    local N_orig = r(N)
    local N_stata = `N_orig'
    if "`wtype'" == "fweight" {
        tempvar wsum
        qui gen double `wsum' = `wvar' if `touse'
        qui summ `wsum' if `touse', meanonly
        local N_stata = round(r(sum))
    }

    // ── Purely categorical model: native collapse + poisson ───────────────
    // When every term is cell-level (p == 0), no fixed effects are absorbed,
    // and the default observed-information VCE is requested, the sufficiency
    // of the cell sums (D_j, T_j) makes estimation on the collapsed cells
    // exact -- coefficients, standard errors, and the log likelihood all
    // reproduce the microdata poisson. Stata's compiled collapse + poisson
    // on a few hundred cells is faster than iterating in Mata over the
    // microdata, so this case is routed there. Robust and cluster-robust
    // VCEs need the per-observation microdata scores (the cell residual sum
    // squared is not the sum of squared residuals), so they fall through to
    // the Mata engine below.
    if (`p' == 0 & !`hasfe' & "`ctype'" == "oim") {
        if "`nolog'" == "" {
            di as txt "  (fast) `k' cell-level parameters, collapse + poisson..."
        }
        tempvar dltv
        qui gen double `dltv' = _d * ln(`exposure') if `touse'
        tempname b_c V_c dlogt sumD sumT llc
        preserve
            qui keep if `touse'
            local cwgt ""
            if "`wvar'" != "" local cwgt "[`wtype' = `wvar']"
            collapse (sum) _d `exposure' `dltv' `cwgt', by(`catnames')
            local N_cells = _N
            qui summ `dltv', meanonly
            scalar `dlogt' = r(sum)
            qui summ _d, meanonly
            scalar `sumD' = r(sum)
            qui summ `exposure', meanonly
            scalar `sumT' = r(sum)
            qui poisson _d `varlist', exposure(`exposure') `constant' nolog
            matrix `b_c' = e(b)
            matrix `V_c' = e(V)
            local df_m = e(df_m)
            local rank = e(rank)
            local iter = e(ic)
            local conv = e(converged)
            // individual-level log likelihood from the cell linear predictors
            // (etac_j = x_j'beta, no offset): ll = sum d_i ln t_i
            //   + sum_j D_j etac_j - sum_j T_j exp(etac_j)
            tempvar etac
            qui predict double `etac', xb nooffset
            mata: st_numscalar("`llc'",                                  ///
                st_numscalar("`dlogt'")                                  ///
              + quadsum(st_data(., "_d") :* st_data(., "`etac'"))        ///
              - quadsum(st_data(., "`exposure'") :* exp(st_data(., "`etac'"))))
        restore
        local ll   = scalar(`llc')
        local a0   = ln(scalar(`sumD') / scalar(`sumT'))
        local ll_0 = scalar(`dlogt') + `a0' * scalar(`sumD') - scalar(`sumD')

        // strip equation names, then post exactly as the Mata path does
        local ncols = colsof(`b_c')
        local blankeq
        forvalues i = 1/`ncols' {
            local blankeq `blankeq' _
        }
        matrix coleq `b_c' = `blankeq'
        matrix coleq `V_c' = `blankeq'
        matrix roweq `V_c' = `blankeq'

        ereturn clear
        ereturn post `b_c' `V_c', esample(`touse') obs(`N_stata')
        ereturn scalar ll        = `ll'
        ereturn scalar ll_0      = `ll_0'
        ereturn scalar chi2      = 2 * (`ll' - `ll_0')
        ereturn scalar df_m      = `df_m'
        ereturn scalar rank      = `rank'
        ereturn scalar ic        = `iter'
        ereturn scalar converged = `conv'
        ereturn scalar N_cells   = `N_cells'
        ereturn local  vce       "oim"
        ereturn local  cont_vars ""
        ereturn local  cat_vars  "`catnames'"
        if "`nolog'" == "" {
            di as txt "  (fast) `N_orig' obs, `N_cells' cells (collapse + poisson)"
        }
        exit
    }

    // ── Estimate in Mata ───────────────────────────────────────────────────
    if "`nolog'" == "" {
        di as txt "  (fast) `p' individual-level + `k' cell-level parameters..."
    }

    tempname b_mat V_mat
    scalar _stpf_ll    = .
    scalar _stpf_ll0   = .
    scalar _stpf_iter  = .
    scalar _stpf_conv  = .
    scalar _stpf_nclu  = .
    scalar _stpf_rank  = .
    scalar _stpf_J     = .
    scalar _stpf_wald  = .

    mata: _stpois_fast_est(     ///
        "_d",                   ///
        "`exposure'",           ///
        "`touse'",              ///
        "`wvar'",               ///
        `N_stata',              ///
        "`fespec'",             ///
        "`dvar'",               ///
        "`contnames'",          ///
        "`catnames'",           ///
        "`xspec'",              ///
        "`wspec'",              ///
        "`wcat_idx'",           ///
        `hascons',              ///
        `tol',                  ///
        `maxiter',              ///
        "`ctype'",              ///
        "`clusterv'",           ///
        "`b_mat'",              ///
        "`V_mat'",              ///
        "_stpf_ll",             ///
        "_stpf_ll0",            ///
        "_stpf_iter",           ///
        "_stpf_conv",           ///
        "_stpf_nclu",           ///
        "_stpf_rank",           ///
        "_stpf_J",              ///
        "_stpf_wald"            ///
    )

    local ll      = scalar(_stpf_ll)
    local ll_0    = scalar(_stpf_ll0)
    local iter    = scalar(_stpf_iter)
    local conv    = scalar(_stpf_conv)
    local n_clust = scalar(_stpf_nclu)
    local rank    = scalar(_stpf_rank)
    local N_cells = scalar(_stpf_J)
    local wald    = scalar(_stpf_wald)
    scalar drop _stpf_ll _stpf_ll0 _stpf_iter _stpf_conv _stpf_nclu ///
        _stpf_rank _stpf_J _stpf_wald

    // ── Names, chi2, post ──────────────────────────────────────────────────
    local bnames `xterms' `wterms'
    if `hascons' local bnames `bnames' _cons
    matrix colnames `b_mat' = `bnames'
    matrix colnames `V_mat' = `bnames'
    matrix rownames `V_mat' = `bnames'

    local ncols = colsof(`b_mat')
    local blankeq
    forvalues i = 1/`ncols' {
        local blankeq `blankeq' _
    }
    matrix coleq `b_mat' = `blankeq'
    matrix coleq `V_mat' = `blankeq'
    matrix roweq `V_mat' = `blankeq'

    if `hasfe' {
        // Wald test of all reported coefficients (an LR test against a
        // null without the fixed effects would not be valid)
        local df_m = `rank'
        local chi2 = `wald'
    }
    else {
        local df_m = `rank' - `hascons'
        local chi2 = 2 * (`ll' - `ll_0')
    }

    ereturn clear
    ereturn post `b_mat' `V_mat', esample(`touse') obs(`N_stata')

    ereturn scalar ll        = `ll'
    if !`hasfe' {
        ereturn scalar ll_0  = `ll_0'
    }
    ereturn scalar chi2      = `chi2'
    ereturn scalar df_m      = `df_m'
    ereturn scalar rank      = `rank'
    ereturn scalar ic        = `iter'
    ereturn scalar converged = `conv'
    ereturn scalar N_cells   = `N_cells'
    ereturn local  vce       "`ctype'"
    ereturn local  cont_vars "`contnames'"
    ereturn local  cat_vars  "`catnames'"
    if "`ctype'" == "cluster" {
        ereturn scalar N_clust  = `n_clust'
        ereturn local  clustvar "`clusterv'"
    }

    if !`conv' {
        di as txt "  Warning: Newton did not converge in `maxiter' iterations"
    }
    if "`nolog'" == "" {
        di as txt "  (fast) `N_orig' obs, " `N_cells' " cells, `iter' Newton iterations, converged=" `conv'
        if "`ctype'" == "cluster" {
            di as txt "  (Std. err. adjusted for " `n_clust' " clusters in `clusterv')"
        }
    }
end


// ============================================================================
// Mata
// ============================================================================
mata:

// Segment sums over a sorted matrix: column-wise running sum, differenced
// at the segment ends. O(n) regardless of the number of segments.
real matrix _stpf_segsum(real matrix Ap, real colvector ends)
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

// Segment boundaries of a sorted key matrix as an (starts, ends) matrix,
// the same format as panelsetup(). Used where the result must live in a
// pointer (a function result gets its own storage, unlike a reused local).
real matrix _stpf_panels(real matrix Ks)
{
    real colvector starts, ends

    _stpf_segments(Ks, starts, ends)
    return((starts, ends))
}

// Segment boundaries of a sorted key matrix: rows where any key column
// changes start a new segment.
void _stpf_segments(real matrix Ks, real colvector starts, real colvector ends)
{
    real colvector flag
    real scalar    n

    n       = rows(Ks)
    flag    = J(n, 1, 1)
    if (n > 1) {
        flag[|2 \ n|] = rowmax(Ks[|2, 1 \ n, .|] :!= Ks[|1, 1 \ n - 1, .|])
    }
    starts = selectindex(flag)
    if (rows(starts) > 1) ends = starts[|2 \ rows(starts)|] :- 1 \ n
    else                  ends = J(1, 1, n)
}

// Build design columns from a piece-encoded spec ("c<k>" continuous,
// "d<k>_<lvl>" indicator, pieces joined by "*") on source matrices V, C.
real matrix _stpf_build(string scalar spec_s,
                        real matrix V, real matrix C)
{
    string rowvector terms, pieces
    string scalar    pc
    real matrix      OUT
    real colvector   col
    real scalar      nt, n, j, q, kk, lvl

    terms = tokens(spec_s)
    nt    = cols(terms)
    n     = max((rows(V), rows(C)))
    OUT   = J(n, nt, .)

    for (j = 1; j <= nt; j++) {
        pieces = tokens(subinstr(terms[j], "*", " "))
        col    = J(n, 1, 1)
        for (q = 1; q <= cols(pieces); q++) {
            pc = pieces[q]
            if (substr(pc, 1, 1) == "c") {
                kk  = strtoreal(substr(pc, 2, .))
                col = col :* V[., kk]
            }
            else {
                kk  = strtoreal(substr(pc, 2, strpos(pc, "_") - 2))
                lvl = strtoreal(substr(pc, strpos(pc, "_") + 1, .))
                col = col :* (C[., kk] :== lvl)
            }
        }
        OUT[., j] = col
    }
    return(OUT)
}

// 0/1 start flags of length n given segment start positions
real colvector _stpf_startflags(real scalar n, real colvector starts)
{
    real colvector f

    f         = J(n, 1, 0)
    f[starts] = J(rows(starts), 1, 1)
    return(f)
}

// Weighted within-group demeaning using a precomputed sort permutation
// and segment ends: one O(n) pass per call.
real colvector _stpf_wdemean(real colvector x,
                             real colvector w,
                             real colvector perm,
                             real colvector ends)
{
    real colvector xp, wp, wsum, wxsum, gmean, res, out, gidx, starts
    real scalar    nseg, n

    n     = rows(x)
    xp    = x[perm]
    wp    = w[perm]
    nseg  = rows(ends)
    wsum  = _stpf_segsum(wp, ends)
    wxsum = _stpf_segsum(wp :* xp, ends)
    gmean = wxsum :/ (wsum + (wsum :== 0))
    // group index of each sorted row
    starts = J(nseg, 1, 1)
    if (nseg > 1) starts[|2 \ nseg|] = ends[|1 \ nseg - 1|] :+ 1
    gidx = runningsum(_stpf_startflags(n, starts))
    res  = xp - gmean[gidx]
    out       = x
    out[perm] = res
    return(out)
}



// Closed-form Gauss-Seidel update of the absorbed fixed effects given the
// covariate part of the linear predictor: for Poisson the profile optimum
// of group g is alpha_g = ln(sum_g v*d / sum_g v*mu), iterated over the
// absorbed terms until stable.
real colvector _stpf_feupdate(real colvector fe0, real colvector base,
                              real colvector v,   real colvector t,
                              real matrix Pf,
                              pointer(real matrix)    rowvector infosf,
                              pointer(real colvector) rowvector VDf,
                              real matrix Gidx, real scalar nfe,
                              real scalar tol)
{
    real colvector fe, vmu, adj
    real scalar    sweep, kk, maxadj

    fe = fe0
    for (sweep = 1; sweep <= 200; sweep++) {
        maxadj = 0
        for (kk = 1; kk <= nfe; kk++) {
            vmu = v :* t :* exp(base + fe)
            adj = ln(*VDf[kk] :/ _stpf_segsum(
                vmu[Pf[., kk]], (*infosf[kk])[., 2]))
            fe  = fe + adj[Gidx[., kk]]
            if (max(abs(adj)) > maxadj) maxadj = max(abs(adj))
        }
        if (maxadj < tol) break
    }
    return(fe)
}

void _stpois_fast_est(
    string scalar dvar,       string scalar evar,
    string scalar tousevar,
    string scalar wvarname,   real scalar n_stata,
    string scalar fespec_s,   string scalar dfevar,
    string scalar contnames_s, string scalar catnames_s,
    string scalar xspec_s,     string scalar wspec_s,
    string scalar cellcols_s,
    real scalar hascons,
    real scalar tol, real scalar maxiter,
    string scalar ctype, string scalar cluvar,
    string scalar b_ret,   string scalar V_ret,
    string scalar ll_ret,  string scalar ll0_ret,
    string scalar iter_ret, string scalar conv_ret,
    string scalar nclu_ret, string scalar rank_ret,
    string scalar J_ret,    string scalar wald_ret)
{
    real colvector idx, d, t, v, vd, vt, cell, perm, starts, ends, reps
    real matrix    V, C, X, W, Xp, Crep
    real colvector Dj, Tj, dlt
    real scalar    n, p, k, q, J_, iter, halv
    real colvector theta, step, theta_new, etac, eta, mu, mup, r, Mj, Rj, s
    real matrix    MX, H, Hinv, VV, S, B, SG, Xfull, Xt, Xt_prev
    real scalar    ll, ll_new, ll0, a0, sum_d, sum_t, converged, delta_c
    real scalar    rank, G_clust, dlogt, wald
    real colvector cids, permc, starts_c, ends_c
    real colvector cellcols
    string rowvector fes
    real scalar    nfe, kk, sweep, maxadj, inner
    real matrix    Pf, M, Gidx
    real colvector prm, fe, fe_try, base, vmu, adj
    pointer(real matrix)    rowvector infosf
    pointer(real colvector) rowvector VDf
    real scalar    jj

    idx = selectindex(st_data(., tousevar))
    n   = rows(idx)
    d   = st_data(idx, dvar)
    t   = st_data(idx, evar)
    v   = (wvarname != "" ? st_data(idx, wvarname) : J(n, 1, 1))
    vd  = v :* d
    vt  = v :* t

    V = (contnames_s != "" ? st_data(idx, tokens(contnames_s)) : J(n, 0, .))
    C = (catnames_s  != "" ? st_data(idx, tokens(catnames_s))  : J(n, 0, .))

    // ── Cells: one sort on the cell-defining categorical columns ──────────
    cellcols = strtoreal(tokens(cellcols_s))'
    perm     = order(C[., cellcols], (1..rows(cellcols)))
    _stpf_segments(C[perm, cellcols], starts, ends)
    J_       = rows(starts)
    // cell index: cumulative count of segment starts, scattered back
    cell       = J(n, 1, .)
    cell[perm] = runningsum(_stpf_startflags(n, starts))
    reps     = perm[starts]
    Crep     = C[reps, .]

    Dj = _stpf_segsum(vd[perm], ends)
    Tj = _stpf_segsum(vt[perm], ends)

    // ── Absorbed fixed effects: per-term sort, segments, group index,
    //    and (constant) weighted event sums ────────────────────────────────
    fes = tokens(fespec_s, ";")
    fes = select(fes, fes :!= ";")
    nfe = cols(fes)
    fe  = J(n, 1, 0)
    if (nfe > 0) {
        Pf     = J(n, nfe, .)
        Gidx   = J(n, nfe, .)
        infosf = J(1, nfe, NULL)
        VDf    = J(1, nfe, NULL)
        for (kk = 1; kk <= nfe; kk++) {
            M          = st_data(idx, tokens(fes[kk]))
            prm        = order(M, (1..cols(M)))
            Pf[., kk]  = prm
            infosf[kk] = &_stpf_panels(M[prm, .])
            Gidx[prm, kk] = runningsum(
                _stpf_startflags(n, (*infosf[kk])[., 1]))
            VDf[kk]    = &_stpf_segsum(vd[prm], (*infosf[kk])[., 2])
        }
    }

    // ── Designs ────────────────────────────────────────────────────────────
    X = _stpf_build(xspec_s, V, C)               // n × p (p may be 0)
    W = _stpf_build(wspec_s, J(J_, 0, .), Crep)  // J × k
    if (hascons) W = W, J(J_, 1, 1)
    p = cols(X)
    k = cols(W)
    q = p + k

    dlt   = ln(t :+ 1e-300)
    dlogt = quadsum(vd :* dlt)

    // ── Null model (exposure + constant) ───────────────────────────────────
    sum_d = quadsum(vd)
    sum_t = quadsum(vt)
    a0    = ln(sum_d / sum_t)
    ll0   = dlogt + a0 * sum_d - sum_d

    // ── Newton, with FEs concentrated out by Gauss-Seidel ─────────────────
    theta = J(q, 1, 0)
    if (hascons) theta[q] = a0
    converged = 0
    ll = .

    if (p > 0) Xp = X[perm, .]
    // With absorbed FEs the theta update is Newton on the profile
    // likelihood: design columns are weighted-demeaned within the FE
    // groups (Frisch-Waugh-Lovell) and the step uses the demeaned score
    // and information. A conditional-Hessian step would zigzag whenever
    // the covariates are correlated with the FE space.
    if (nfe > 0) Xfull = (p > 0 ? (X, W[cell, .]) : W[cell, .])

    for (iter = 1; iter <= maxiter; iter++) {
        etac = W * theta[|p + 1 \ q|]

        // Concentrate the FEs out at the current theta
        if (nfe > 0) {
            base = (p > 0 ? X * theta[|1 \ p|] : J(n, 1, 0)) + etac[cell]
            fe   = _stpf_feupdate(fe, base, v, t, Pf, infosf, VDf,
                                  Gidx, nfe, tol)
        }

        if (nfe > 0) {
            // Profile Newton: demean the design at the current weights,
            // then score and information on the demeaned columns
            eta = base + fe
            mu  = t :* exp(eta)
            mu  = mu + (mu :< 1e-300) :* 1e-300
            ll  = quadsum(v :* (d :* ln(mu) - mu))
            vmu = v :* mu
            Xt  = Xfull
            for (inner = 1; inner <= 200; inner++) {
                Xt_prev = Xt
                for (kk = 1; kk <= nfe; kk++) {
                    for (jj = 1; jj <= q; jj++) {
                        Xt[., jj] = _stpf_wdemean(Xt[., jj], vmu,
                            Pf[., kk], (*infosf[kk])[., 2])
                    }
                }
                if (max(abs(Xt - Xt_prev)) < tol * 0.01) break
            }
            s = quadcross(Xt, v :* (d - mu))
            H = quadcross(Xt, vmu, Xt)
        }
        else if (p == 0) {
            // Pure cell-level model: iterate on (D_j, T_j) only
            Mj = Tj :* exp(etac)
            ll = dlogt + quadsum(Dj :* etac) - quadsum(Mj)
            Rj = Dj - Mj
            s  = quadcross(W, Rj)
            H  = quadcross(W, Mj, W)
        }
        else {
            eta = X * theta[|1 \ p|] + etac[cell]
            mu  = t :* exp(eta)
            mu  = mu + (mu :< 1e-300) :* 1e-300
            ll  = quadsum(v :* (d :* ln(mu) - mu))
            vmu = v :* mu
            mup = vmu[perm]
            Mj  = _stpf_segsum(mup, ends)
            Rj  = Dj - Mj
            MX  = _stpf_segsum(mup :* Xp, ends)
            s   = quadcross(X, v :* (d - mu)) \ quadcross(W, Rj)
            H   = quadcross(X, vmu, X), (MX' * W) \
                  (MX' * W)', quadcross(W, Mj, W)
        }

        Hinv = invsym(H)
        step = Hinv * s

        theta_new = theta + step
        for (halv = 1; halv <= 30; halv++) {
            etac = W * theta_new[|p + 1 \ q|]
            if (nfe > 0) {
                // Profile likelihood: re-solve the FEs at the trial theta,
                // otherwise good profile steps would be rejected whenever
                // the covariates are correlated with the FE space
                base   = (p > 0 ? X * theta_new[|1 \ p|] : J(n, 1, 0)) +
                         etac[cell]
                fe_try = _stpf_feupdate(fe, base, v, t, Pf, infosf, VDf,
                                        Gidx, nfe, tol)
                mu     = t :* exp(base + fe_try)
                mu     = mu + (mu :< 1e-300) :* 1e-300
                ll_new = quadsum(v :* (d :* ln(mu) - mu))
            }
            else if (p == 0) {
                ll_new = dlogt + quadsum(Dj :* etac) -
                         quadsum(Tj :* exp(etac))
            }
            else {
                eta    = X * theta_new[|1 \ p|] + etac[cell]
                mu     = t :* exp(eta)
                mu     = mu + (mu :< 1e-300) :* 1e-300
                ll_new = quadsum(v :* (d :* ln(mu) - mu))
            }
            if (ll_new >= ll - 1e-12) break
            step      = step / 2
            theta_new = theta + step
        }
        if (nfe > 0) fe = fe_try

        delta_c = max(abs(theta_new - theta) :/ (1 :+ abs(theta_new)))
        theta   = theta_new
        ll      = ll_new
        if (delta_c < tol) {
            converged = 1
            break
        }
    }

    // ── Final FE refresh, mu, and log likelihood ───────────────────────────
    etac = W * theta[|p + 1 \ q|]
    if (nfe > 0) {
        base = (p > 0 ? X * theta[|1 \ p|] : J(n, 1, 0)) + etac[cell]
        fe   = _stpf_feupdate(fe, base, v, t, Pf, infosf, VDf,
                              Gidx, nfe, tol)
    }
    eta = (p > 0 ? X * theta[|1 \ p|] : J(n, 1, 0)) + etac[cell] + fe
    mu  = t :* exp(eta)
    mu  = mu + (mu :< 1e-300) :* 1e-300
    r   = d - mu
    ll  = quadsum(v :* (d :* ln(mu) - mu))
    vmu = v :* mu

    // Store the FE contribution for predict
    if (nfe > 0 & dfevar != "") st_store(idx, dfevar, fe)

    // ── Final Hessian and VCE ──────────────────────────────────────────────
    if (nfe == 0) {
        mup = vmu[perm]
        Mj  = _stpf_segsum(mup, ends)
        if (p > 0) {
            MX = _stpf_segsum(mup :* Xp, ends)
            H  = quadcross(X, vmu, X), (MX' * W) \
                 (MX' * W)', quadcross(W, Mj, W)
        }
        else {
            H = quadcross(W, Mj, W)
        }
        Hinv = invsym(H)
        rank = sum(diagonal(Hinv) :!= 0)
        S    = (ctype == "oim" ? J(0, 0, .) :
               (p > 0 ? (v :* r) :* (X, W[cell, .]) : (v :* r) :* W[cell, .]))
    }
    else {
        // Frisch-Waugh-Lovell: weighted-demean every design column within
        // the absorbed groups; the partitioned information and the scores
        // on the demeaned columns give the exact VCE (as in ppmlhdfe)
        Xfull = (p > 0 ? (X, W[cell, .]) : W[cell, .])
        Xt    = Xfull
        for (inner = 1; inner <= 200; inner++) {
            Xt_prev = Xt
            for (kk = 1; kk <= nfe; kk++) {
                for (jj = 1; jj <= q; jj++) {
                    Xt[., jj] = _stpf_wdemean(Xt[., jj], vmu,
                        Pf[., kk], (*infosf[kk])[., 2])
                }
            }
            if (max(abs(Xt - Xt_prev)) < tol * 0.01) break
        }
        H    = quadcross(Xt, vmu, Xt)
        Hinv = invsym(H)
        rank = sum(diagonal(Hinv) :!= 0)
        S    = (ctype == "oim" ? J(0, 0, .) : (v :* r) :* Xt)
    }

    G_clust = 0
    if (ctype == "oim") {
        VV = Hinv
    }
    else if (ctype == "robust") {
        // HC1 with Stata's N convention (sum of fweights, else # of obs)
        B  = quadcross(S, S) * (n_stata / (n_stata - 1))
        VV = Hinv * B * Hinv
    }
    else {
        // Cluster sandwich via one sort and segment sums
        cids  = st_data(idx, cluvar)
        permc = order(cids, 1)
        _stpf_segments(cids[permc], starts_c, ends_c)
        G_clust = rows(starts_c)
        SG = _stpf_segsum(S[permc, .], ends_c)
        B  = quadcross(SG, SG) * (G_clust / (G_clust - 1))
        VV = Hinv * B * Hinv
    }

    // Wald chi2 of all reported coefficients (used when FEs are absorbed)
    wald = (nfe > 0 ? theta' * invsym(VV) * theta : 0)

    st_matrix(b_ret, theta')
    st_matrix(V_ret, VV)
    st_numscalar(ll_ret,   ll)
    st_numscalar(ll0_ret,  ll0)
    st_numscalar(iter_ret, iter)
    st_numscalar(conv_ret, converged)
    st_numscalar(nclu_ret, G_clust)
    st_numscalar(rank_ret, rank)
    st_numscalar(J_ret,    J_)
    st_numscalar(wald_ret, wald)
}

end
