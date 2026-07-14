*! _stpois_fast  version 0.6.0  14jul2026
*! Exact Poisson EHA accelerated through the cell structure of the
*! categorical covariates.
*!
*! Terms in the varlist are classified as cell-level (built only from
*! factor-variable pieces, including factor#factor interactions) or
*! individual-level (involving at least one continuous piece, including
*! continuous#continuous and continuous#factor interactions). Newton
*! iterations compute the score and Hessian from exponentially tilted
*! within-cell moments of the individual-level columns: one O(n p^2)
*! pass over the microdata per iteration, all remaining algebra on the
*! J cells. When every term is cell-level, iterations run entirely on
*! the collapsed (events, exposure) cell sums at O(J) per iteration.
*!
*! Data handling is done in a single Mata pass with one sort:
*! group structures come from order() plus running-sum segment
*! aggregation, with no preserve/collapse round-trips.
*!
*! Estimates are the exact individual-level MLE: coefficients, log
*! likelihood, and OIM/robust/cluster VCEs match poisson numerically.
*!
*! Called by stpois when fast is specified. Not intended for direct use.
program _stpois_fast, eclass
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

    local hascons = cond("`constant'" == "noconstant", 0, 1)

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
        exit 198
    }

    local p : word count `xterms'
    local k : word count `wterms'

    // ── Count obs ──────────────────────────────────────────────────────────
    qui count if `touse'
    local N_orig = r(N)

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

    mata: _stpois_fast_est(     ///
        "_d",                   ///
        "`exposure'",           ///
        "`touse'",              ///
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
        "_stpf_J"               ///
    )

    local ll      = scalar(_stpf_ll)
    local ll_0    = scalar(_stpf_ll0)
    local iter    = scalar(_stpf_iter)
    local conv    = scalar(_stpf_conv)
    local n_clust = scalar(_stpf_nclu)
    local rank    = scalar(_stpf_rank)
    local N_cells = scalar(_stpf_J)
    scalar drop _stpf_ll _stpf_ll0 _stpf_iter _stpf_conv _stpf_nclu _stpf_rank _stpf_J

    // ── Names, chi2, post ──────────────────────────────────────────────────
    local bnames `xterms' `wterms'
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


void _stpois_fast_est(
    string scalar dvar,       string scalar evar,
    string scalar tousevar,
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
    string scalar J_ret)
{
    real colvector idx, d, t, cell, perm, starts, ends, reps
    real matrix    V, C, X, W, Xp, Crep
    real colvector Dj, Tj, dlt
    real scalar    n, p, k, q, J_, iter, halv
    real colvector theta, step, theta_new, etac, eta, mu, mup, r, Mj, Rj, s
    real matrix    MX, H, Hinv, VV, S, B, SG
    real scalar    ll, ll_new, ll0, a0, sum_d, sum_t, converged, delta_c
    real scalar    rank, G_clust, dlogt
    real colvector cids, permc, starts_c, ends_c
    real colvector cellcols

    idx = selectindex(st_data(., tousevar))
    n   = rows(idx)
    d   = st_data(idx, dvar)
    t   = st_data(idx, evar)

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

    Dj = _stpf_segsum(d[perm], ends)
    Tj = _stpf_segsum(t[perm], ends)

    // ── Designs ────────────────────────────────────────────────────────────
    X = _stpf_build(xspec_s, V, C)               // n × p (p may be 0)
    W = _stpf_build(wspec_s, J(J_, 0, .), Crep)  // J × k
    if (hascons) W = W, J(J_, 1, 1)
    p = cols(X)
    k = cols(W)
    q = p + k

    dlt   = ln(t :+ 1e-300)
    dlogt = quadsum(d :* dlt)

    // ── Null model (exposure + constant) ───────────────────────────────────
    sum_d = quadsum(d)
    sum_t = quadsum(t)
    a0    = ln(sum_d / sum_t)
    ll0   = dlogt + a0 * sum_d - sum_d

    // ── Newton ─────────────────────────────────────────────────────────────
    theta = J(q, 1, 0)
    if (hascons) theta[q] = a0
    converged = 0
    ll = .

    if (p > 0) Xp = X[perm, .]

    for (iter = 1; iter <= maxiter; iter++) {
        etac = W * theta[|p + 1 \ q|]
        if (p > 0) {
            eta = X * theta[|1 \ p|] + etac[cell]
            mu  = t :* exp(eta)
            mu  = mu + (mu :< 1e-300) :* 1e-300
            ll  = quadsum(d :* ln(mu) - mu)
            mup = mu[perm]
            Mj  = _stpf_segsum(mup, ends)
            MX  = _stpf_segsum(mup :* Xp, ends)
            Rj  = Dj - Mj
            s   = quadcross(X, d - mu) \ quadcross(W, Rj)
            H   = quadcross(X, mu, X), (MX' * W) \
                  (MX' * W)', quadcross(W, Mj, W)
        }
        else {
            // Pure cell-level model: iterate on (D_j, T_j) only
            Mj = Tj :* exp(etac)
            ll = dlogt + quadsum(Dj :* etac) - quadsum(Mj)
            Rj = Dj - Mj
            s  = quadcross(W, Rj)
            H  = quadcross(W, Mj, W)
        }

        Hinv = invsym(H)
        step = Hinv * s

        theta_new = theta + step
        for (halv = 1; halv <= 30; halv++) {
            etac = W * theta_new[|p + 1 \ q|]
            if (p > 0) {
                eta    = X * theta_new[|1 \ p|] + etac[cell]
                mu     = t :* exp(eta)
                mu     = mu + (mu :< 1e-300) :* 1e-300
                ll_new = quadsum(d :* ln(mu) - mu)
            }
            else {
                ll_new = dlogt + quadsum(Dj :* etac) - quadsum(Tj :* exp(etac))
            }
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

    // ── Final Hessian and VCE ──────────────────────────────────────────────
    etac = W * theta[|p + 1 \ q|]
    if (p > 0) {
        eta = X * theta[|1 \ p|] + etac[cell]
    }
    else {
        eta = etac[cell]
    }
    mu = t :* exp(eta)
    mu = mu + (mu :< 1e-300) :* 1e-300
    r  = d - mu
    ll = quadsum(d :* ln(mu) - mu)

    mup = mu[perm]
    Mj  = _stpf_segsum(mup, ends)
    if (p > 0) {
        MX = _stpf_segsum(mup :* Xp, ends)
        H  = quadcross(X, mu, X), (MX' * W) \
             (MX' * W)', quadcross(W, Mj, W)
    }
    else {
        H = quadcross(W, Mj, W)
    }
    Hinv = invsym(H)
    rank = sum(diagonal(Hinv) :!= 0)

    G_clust = 0
    if (ctype == "oim") {
        VV = Hinv
    }
    else {
        S = (p > 0 ? (r :* X, r :* W[cell, .]) : r :* W[cell, .])
        if (ctype == "robust") {
            B  = quadcross(S, S) * (n / (n - 1))
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
    }

    st_matrix(b_ret, theta')
    st_matrix(V_ret, VV)
    st_numscalar(ll_ret,   ll)
    st_numscalar(ll0_ret,  ll0)
    st_numscalar(iter_ret, iter)
    st_numscalar(conv_ret, converged)
    st_numscalar(nclu_ret, G_clust)
    st_numscalar(rank_ret, rank)
    st_numscalar(J_ret,    J_)
}

// 0/1 start flags of length n given segment start positions
real colvector _stpf_startflags(real scalar n, real colvector starts)
{
    real colvector f

    f         = J(n, 1, 0)
    f[starts] = J(rows(starts), 1, 1)
    return(f)
}

end
