*! version 0.1.0  13jul2026  Andreas Ljungström, SOFI Stockholm University
program stpois_p
    version 14

    syntax newvarname [if] [in] [,   ///
        XB                           ///
        N                            ///
        Hazard                       ///
        HR                           ///
        SURVival                     ///
        noOFFset                     ///
        SCores                       ///
        ]

    // Exactly one statistic
    local nopt : word count `xb' `n' `hazard' `hr' `survival' `scores'
    if `nopt' > 1 {
        di as error "only one statistic may be specified"
        exit 198
    }
    if `nopt' == 0 local xb xb

    // scores: one score variable per parameter
    if "`scores'" != "" {
        di as error "scores not supported; use predict with xb, n, hazard, or surv"
        exit 198
    }

    local newv `varlist'

    // Copy e(b) to a named tempmatrix before scoring (avoids r(111) inside margins)
    tempname bmat
    matrix `bmat' = e(b)

    // Linear predictor Xβ (no offset)
    tempvar xbhat
    matrix score double `xbhat' = `bmat' `if' `in'

    if "`xb'" != "" {
        // `offset' = "nooffset" when user specifies nooffset; "" otherwise
        if "`offset'" != "" {
            gen double `newv' = `xbhat' `if' `in'
        }
        else {
            // xb + ln(exposure)
            gen double `newv' = `xbhat' + ln(_t - _t0) `if' `in'
        }
    }
    else if "`n'" != "" {
        // Expected events = exp(Xβ) * exposure
        gen double `newv' = exp(`xbhat') * (_t - _t0) `if' `in'
    }
    else if "`hazard'" != "" | "`hr'" != "" {
        // Hazard rate = exp(Xβ)
        gen double `newv' = exp(`xbhat') `if' `in'
    }
    else if "`survival'" != "" {
        // S(t) = exp(-exp(Xβ) * (t - t0))
        gen double `newv' = exp(-exp(`xbhat') * (_t - _t0)) `if' `in'
    }
end
