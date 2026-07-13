*! _stpois_fast_offset  version 0.2.0  13jul2026
*! Two-stage multiplicative offset method for fast Poisson EHA
*! Stage 1: estimate continuous-covariate effects on individual data
*! Stage 2: mutate exposure, collapse to categorical cells, re-estimate
*! Called by stpois when fast(offset) is specified. Not intended for direct use.
*!
*! APPROXIMATION WARNING: SEs do not account for first-stage estimation
*! uncertainty. Bias grows when continuous covariates correlate with
*! categorical grouping variables. See stpois help for details.
program _stpois_fast_offset, eclass
    version 14

    syntax varlist(fv ts) [aw iw fw pw], ///
        touse(varname)                    ///
        exposure(varname)                 ///
        [nolog                            ///
        noconstant                        ///
        vce(string)                       ///
        robust                            ///
        cluster(varname)                  ///
        MTopel                            ///
        *]

    // ── Parse varlist: categorical (i.) vs continuous (bare) ───────────────
    local cont_vars ""     // raw variable names, continuous
    local cat_vars_raw ""  // raw variable names, categorical (for by())
    local cat_terms ""     // factor-variable terms (for model)

    foreach tok of local varlist {
        if strpos("`tok'", "#") > 0 {
            di as error "Interactions not supported with fast(offset)"
            exit 198
        }
        // Factor variable: starts with i. or matches n. or nb. level notation
        if regexm("`tok'", "^[0-9]*b?i\.") | substr("`tok'", 1, 2) == "i." {
            // Strip factor notation to get raw variable name
            local rawname = regexr("`tok'", "^[0-9]*b?i\.", "")
            local cat_vars_raw `cat_vars_raw' `rawname'
            local cat_terms    `cat_terms' `tok'
        }
        else {
            // Bare name, c., L., D., F. — treat as continuous
            // Strip c. prefix if present
            local rawname = regexr("`tok'", "^[a-zA-Z]+\.", "")
            if "`rawname'" == "" local rawname `tok'
            // For bare names, no stripping needed
            if !regexm("`tok'", "\.") local rawname `tok'
            local cont_vars `cont_vars' `rawname'
        }
    }

    if "`cont_vars'" == "" {
        di as error "fast(offset) requires at least one continuous covariate (bare variable name)"
        di as error "Example: stpois age income i.edu i.region, fast(offset)"
        exit 198
    }
    if "`cat_vars_raw'" == "" {
        di as error "fast(offset) requires at least one categorical covariate (i.varname)"
        di as error "Without categorical variables there is nothing to collapse by"
        exit 198
    }

    // Build VCE pass-through
    local vceopts
    if "`robust'"  != "" local vceopts `vceopts' robust
    if "`vce'"     != "" local vceopts `vceopts' vce(`vce')
    if "`cluster'" != "" local vceopts `vceopts' cluster(`cluster')

    // ── Count observations before preserve ────────────────────────────────
    qui count if `touse'
    local N_orig = r(N)

    // ── Stage 1: estimate continuous effects on individual data ────────────
    // Include categorical indicators to avoid OVB from correlation.
    // With mtopel, run stage 1 under the requested VCE so its joint
    // covariance matrix can supply the corrected second-stage VCE block.
    if "`nolog'" == "" di as txt "  (fast offset) Stage 1: individual-level Poisson..."
    local s1_vceopts
    if "`mtopel'" != "" local s1_vceopts `vceopts'
    qui poisson _d `cont_vars' `cat_terms' if `touse', ///
        exposure(`exposure') nolog `s1_vceopts' `options'

    // Extract continuous-covariate coefficients (and, for mtopel, the
    // full joint VCE)
    tempname b_stage1 V_stage1
    matrix `b_stage1' = e(b)
    if "`mtopel'" != "" {
        matrix `V_stage1' = e(V)
        local ncols1 = colsof(`V_stage1')
        local blankeq1
        forvalues i = 1/`ncols1' {
            local blankeq1 `blankeq1' ""
        }
        matrix coleq `V_stage1' = `blankeq1'
        matrix roweq `V_stage1' = `blankeq1'
    }
    local nc : word count `cont_vars'
    // Coefficients for cont_vars are the first nc columns of b
    // (continuous vars were listed first in the stage-1 model)
    local s1_xb ""
    forvalues ci = 1/`nc' {
        local vn : word `ci' of `cont_vars'
        local coef = `b_stage1'[1, `ci']
        if `ci' == 1 local s1_xb "`coef' * `vn'"
        else         local s1_xb "`s1_xb' + `coef' * `vn'"
    }

    // ── Preserve, collapse, Stage 2 ───────────────────────────────────────
    preserve

    // Generate structural weight and mutated exposure
    tempvar struct_w mutated_exp
    qui gen double `struct_w'    = exp(`s1_xb') if `touse'
    qui gen double `mutated_exp' = `exposure' * `struct_w' if `touse'

    // Keep only touse obs before collapse
    qui keep if `touse'

    // Collapse: sum events and mutated exposure within categorical cells
    if "`nolog'" == "" di as txt "  (fast offset) Stage 2: collapsing to categorical cells..."
    collapse (sum) _d (sum) `mutated_exp', by(`cat_vars_raw')

    local N_cells = _N

    // Re-estimate on collapsed data
    if "`nolog'" == "" di as txt "  (fast offset) Stage 2: Poisson on collapsed data..."
    qui poisson _d `cat_terms' if _d < ., ///
        exposure(`mutated_exp') nolog `vceopts' `options'

    local N2    = e(N)
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

    // ── Murphy–Topel correction ────────────────────────────────────────────
    // Stage 1 estimates the full model (continuous + categorical) jointly,
    // so the MT-corrected covariance of the second-stage parameters reduces
    // to the corresponding block of the stage-1 joint VCE. Extract it by name.
    if "`mtopel'" != "" {
        // Overwrite V element-wise from the stage-1 joint VCE, matching
        // parameters by name; the stripe of V (incl. omit flags) is untouched.
        local s2names : colnames `b'
        local n2 : word count `s2names'
        local ok 1
        forvalues i = 1/`n2' {
            local ni : word `i' of `s2names'
            local pi = colnumb(`V_stage1', "`ni'")
            if `pi' >= . {
                local ok 0
                continue, break
            }
            forvalues j = `i'/`n2' {
                local nj : word `j' of `s2names'
                local pj = colnumb(`V_stage1', "`nj'")
                if `pj' >= . {
                    local ok 0
                    continue, break
                }
                matrix `V'[`i', `j'] = `V_stage1'[`pi', `pj']
                matrix `V'[`j', `i'] = `V_stage1'[`pj', `pi']
            }
            if !`ok' continue, break
        }
        if !`ok' {
            di as txt "  (mtopel: could not match stage-1 and stage-2 parameter names; " ///
                      "conditional SEs reported)"
            local mtopel
        }
    }

    restore

    // ── Post (using original touse for esample) ────────────────────────────
    ereturn clear
    ereturn post `b' `V', esample(`touse') obs(`N_orig')

    ereturn scalar ll        = `ll'
    ereturn scalar ll_0      = `ll_0'
    ereturn scalar chi2      = `chi2'
    ereturn scalar df_m      = `df_m'
    ereturn scalar rank      = `rank'
    ereturn scalar converged = `conv'
    ereturn scalar N_cells   = `N_cells'
    ereturn scalar mtopel    = ("`mtopel'" != "")
    ereturn local  vce       "`vce_r'"
    ereturn local  cont_vars "`cont_vars'"
    ereturn local  cat_vars  "`cat_vars_raw'"

    if "`nolog'" == "" {
        di as txt "  (fast offset) `N_orig' individual obs → `N_cells' cells"
    }
end
