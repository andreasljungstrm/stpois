*! _stpois_fast_moments  version 0.2.0  13jul2026
*! CGF/Jensen moment-correction method for fast Poisson EHA
*!
*! Collapses individual-level data to categorical cells, computes within-cell
*! means, variances, and covariances of continuous covariates, then estimates
*! Poisson including these moment terms as aggregation-bias corrections.
*!
*! Basis: 2nd-order CGF expansion  log E[exp(γ'X)] ≈ γ'μ + ½γ'Σγ
*!
*! APPROXIMATION WARNINGS:
*!  · SEs are conditional on moment estimates; first-stage uncertainty not propagated
*!  · Stable when continuous X is approximately normal within cells
*!  · If coef(var_xi) ≈ γi²/2 the 2nd-order expansion holds well
*!
*! Called by stpois when fast(moments) is specified. Not intended for direct use.
program _stpois_fast_moments, eclass
    version 14

    syntax varlist(fv ts) [aw iw fw pw], ///
        touse(varname)                    ///
        exposure(varname)                 ///
        [nolog                            ///
        noconstant                        ///
        vce(string)                       ///
        robust                            ///
        cluster(varname)                  ///
        SKEWness                          ///
        *]

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

    local nc : word count `cont_vars'
    tokenize `cont_vars'   // now ``1''=first cont var, ``2''=second, etc.

    // ── VCE pass-through ──────────────────────────────────────────────────
    local vceopts
    if "`robust'"  != "" local vceopts `vceopts' robust
    if "`vce'"     != "" local vceopts `vceopts' vce(`vce')
    if "`cluster'" != "" local vceopts `vceopts' cluster(`cluster')

    // ── Count obs before preserve ──────────────────────────────────────────
    qui count if `touse'
    local N_orig = r(N)

    // ── Generate moment ingredients on individual data ─────────────────────
    preserve
    qui keep if `touse'

    // Squared terms (for variances)
    local sq_names ""
    forvalues ci = 1/`nc' {
        tempvar sq_`ci'
        qui gen double `sq_`ci'' = ``ci''^2
        local sq_names `sq_names' msq_x`ci'=`sq_`ci''
    }

    // Cubed terms (for skewness correction, optional)
    local cube_names ""
    if "`skewness'" != "" {
        forvalues ci = 1/`nc' {
            tempvar cu_`ci'
            qui gen double `cu_`ci'' = ``ci''^3
            local cube_names `cube_names' mcu_x`ci'=`cu_`ci''
        }
    }

    // Cross-product terms (for covariances)
    local cross_names ""
    if `nc' >= 2 {
        forvalues ci = 1/`nc' {
            local cj_start = `ci' + 1
            forvalues cj = `cj_start'/`nc' {
                tempvar cr_`ci'_`cj'
                qui gen double `cr_`ci'_`cj'' = ``ci'' * ``cj''
                local cross_names `cross_names' mcr_x`ci'x`cj'=`cr_`ci'_`cj''
            }
        }
    }

    // Mean names for collapse
    local mean_names ""
    forvalues ci = 1/`nc' {
        local mean_names `mean_names' mx`ci'=``ci''
    }

    // ── Person-time variable (real column required for collapse) ───────────
    tempvar ptime
    qui gen double `ptime' = `exposure'

    // ── One-pass collapse ──────────────────────────────────────────────────
    if "`nolog'" == "" di as txt "  (fast moments) collapsing to categorical cells..."

    // Build collapse command piece by piece
    local collapse_cmd "collapse (sum) _d (sum) `ptime' (mean) `mean_names' (mean) `sq_names'"
    if "`skewness'" != "" & "`cube_names'" != "" {
        local collapse_cmd "`collapse_cmd' (mean) `cube_names'"
    }
    if "`cross_names'" != "" {
        local collapse_cmd "`collapse_cmd' (mean) `cross_names'"
    }
    local collapse_cmd "`collapse_cmd', by(`cat_vars_raw')"
    `collapse_cmd'

    local N_cells = _N

    // ── Reconstruct within-cell moments ───────────────────────────────────
    // Var(X_i)   = E[X²] - (E[X])²
    // Cov(X_i,X_j) = E[X_i*X_j] - E[X_i]*E[X_j]
    // Skew_3(X_i) = E[X³] - 3*E[X]*Var(X) - (E[X])³  (central 3rd moment)

    local var_list  ""
    local cov_list  ""
    local skew_list ""

    forvalues ci = 1/`nc' {
        qui gen double var_x`ci' = msq_x`ci' - (mx`ci')^2
        local var_list `var_list' var_x`ci'

        if "`skewness'" != "" {
            qui gen double skew_x`ci' = mcu_x`ci' - 3*mx`ci'*var_x`ci' - (mx`ci')^3
            local skew_list `skew_list' skew_x`ci'
        }
    }

    if `nc' >= 2 {
        forvalues ci = 1/`nc' {
            local cj_start = `ci' + 1
            forvalues cj = `cj_start'/`nc' {
                qui gen double cov_x`ci'x`cj' = mcr_x`ci'x`cj' - mx`ci'*mx`cj'
                local cov_list `cov_list' cov_x`ci'x`cj'
            }
        }
    }

    // Rename mx* to m_originalname for cleaner output
    local mean_model_vars ""
    forvalues ci = 1/`nc' {
        rename mx`ci' m_``ci''
        local mean_model_vars `mean_model_vars' m_``ci''
    }

    // ── Collapsed Poisson with moment corrections ──────────────────────────
    local moment_terms `mean_model_vars' `var_list' `cov_list' `skew_list'

    if "`nolog'" == "" di as txt "  (fast moments) Poisson on `N_cells' cells..."

    qui poisson _d `cat_terms' `moment_terms', ///
        exposure(`ptime') nolog `vceopts' `options'

    local ll    = e(ll)
    local ll_0  = e(ll_0)
    local chi2  = e(chi2)
    local df_m  = e(df_m)
    local rank  = e(rank)
    local vce_r = e(vce)
    local conv  = e(converged)

    tempname b V
    matrix `b' = e(b)
    matrix `V' = e(V)

    local ncols = colsof(`b')
    local blankeq
    forvalues i = 1/`ncols' {
        local blankeq `blankeq' ""
    }
    matrix coleq `b' = `blankeq'
    matrix coleq `V' = `blankeq'
    matrix roweq `V' = `blankeq'

    restore

    // ── Post (original touse for esample) ─────────────────────────────────
    ereturn clear
    ereturn post `b' `V', esample(`touse') obs(`N_orig')

    ereturn scalar ll        = `ll'
    ereturn scalar ll_0      = `ll_0'
    ereturn scalar chi2      = `chi2'
    ereturn scalar df_m      = `df_m'
    ereturn scalar rank      = `rank'
    ereturn scalar converged = `conv'
    ereturn scalar N_cells   = `N_cells'
    ereturn local  vce       "`vce_r'"
    ereturn local  cont_vars "`cont_vars'"
    ereturn local  cat_vars  "`cat_vars_raw'"

    if "`nolog'" == "" {
        di as txt "  (fast moments) `N_orig' obs → `N_cells' cells"
        di as txt "  Structural γ: coefficients on m_* variables"
        di as txt "  Stability check: coef(var_xi) ≈ γi²/2 → 2nd-order approx holds"
    }
end
