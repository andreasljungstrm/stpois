{smcl}
{* *! version 0.6.0  14jul2026}{...}
{viewerjumpto "Syntax" "stpois##syntax"}{...}
{viewerjumpto "Description" "stpois##description"}{...}
{viewerjumpto "Options" "stpois##options"}{...}
{viewerjumpto "Fast estimation" "stpois##fast"}{...}
{viewerjumpto "HDFE" "stpois##hdfe"}{...}
{viewerjumpto "Postestimation" "stpois##postestimation"}{...}
{viewerjumpto "Examples" "stpois##examples"}{...}
{viewerjumpto "References" "stpois##references"}{...}
{title:Title}

{phang}
{bf:stpois} {hline 2} Poisson event-history regression for {cmd:stset} data


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:stpois}
{varlist}
{ifin}
{weight}
[{cmd:,} {it:options}]

{synoptset 26 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Model}
{synopt:{opt irr}}display incidence-rate ratios (= hazard ratios){p_end}
{synopt:{opt nolog}}suppress iteration log{p_end}
{synopt:{opt level(#)}}set confidence level; default is {cmd:level(95)}{p_end}
{synopt:{opt nocons:tant}}suppress constant term{p_end}
{synopt:{opth vce(vcetype)}}{it:vcetype}: {opt oim}, {opt r:obust},
    {opt cl:uster} {it:clustvar}; all exact on every path{p_end}
{synopt:{opt robust}}synonym for {cmd:vce(robust)}{p_end}
{synopt:{opth cluster(varname)}}synonym for {cmd:vce(cluster} {it:varname}{cmd:)}{p_end}

{syntab:Fast estimation on collapsed cells}
{synopt:{opt fast}}exact MLE accelerated through the cell structure of the
    categorical covariates; see {help stpois##fast:Fast estimation}{p_end}

{syntab:High-dimensional fixed effects}
{synopt:{opt absorb(terms)}}absorb fixed effects; {it:terms} are variable
    names or {it:var1}{cmd:#}{it:var2} interactions;
    see {help stpois##hdfe:HDFE}{p_end}
{synopt:{opt tol:erance(#)}}iteration convergence tolerance; default {cmd:1e-8}{p_end}
{synopt:{opt maxiter(#)}}maximum iterations; default {cmd:100}{p_end}

{synopt:{it:poisson_options}}other options passed to {helpb poisson} (standard path){p_end}
{synoptline}
{p2colreset}{...}

{p 4 6 2}
{opt fweight}s, {opt iweight}s, and {opt pweight}s are allowed;
see {help weight}.{p_end}

{p 4 6 2}
Data must be {cmd:stset} before using {cmd:stpois}; see {helpb stset}.{p_end}

{p 4 6 2}
{cmd:fast} and {cmd:absorb()} may not be combined.{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:stpois} fits a Poisson regression model to {cmd:stset} survival data.
Each record contributes its event indicator ({cmd:_d}) as outcome and
time at risk ({cmd:_t} - {cmd:_t0}) as exposure:

{pmore}
log E[_d{sub:i}] = x{sub:i}β + log(t{sub:i} − t{sub:0i})

{pstd}
This corresponds to a constant-hazard (exponential) model with hazard rate
exp(x{sub:i}β); duration dependence and period effects are accommodated by
episode splitting ({helpb stsplit}) followed by indicator variables for the
time bands.

{pstd}
{bf:Three estimation paths} are available, all producing the exact
maximum-likelihood estimates:

{phang2}
1. {bf:Standard} (default). Wraps Stata's {cmd:poisson} with automatic
exposure. Identical coefficients and standard errors to
{cmd:streg, distribution(exponential)}.

{phang2}
2. {bf:fast}. Exact MLE accelerated through the cell structure of the
categorical covariates. See {help stpois##fast:Fast estimation} below.

{phang2}
3. {bf:absorb()}. High-dimensional fixed-effect absorption via iteratively
reweighted alternating projections. See {help stpois##hdfe:HDFE} below.

{pstd}
Factor-variable interactions are supported on all paths, including
continuous#continuous, continuous#categorical, and
categorical#categorical terms, and interacted fixed effects in
{cmd:absorb()}.


{marker options}{...}
{title:Options}

{phang}
{opt irr} reports coefficients as incidence-rate ratios (hazard ratios).
Affects display only.

{phang}
{opt nolog} suppresses iteration logs.

{phang}
{opt level(#)} sets the confidence level as a percentage. Default 95.

{phang}
{opt noconstant} suppresses the constant term.

{phang}
{opth vce(vcetype)} specifies the standard-error type. {opt oim} (the
default), {opt robust}, and {opt cluster} {it:clustvar} are supported on
all paths and are numerically exact relative to the corresponding
{cmd:poisson} model. On the standard path, any {it:vcetype} accepted by
{cmd:poisson} may be used.

{phang}
{opt robust}, {opth cluster(varname)} are convenient synonyms.

{phang}
{opt fast} selects cell-accelerated exact estimation;
see {help stpois##fast:Fast estimation}.

{phang}
{opt absorb(terms)} absorbs fixed effects. Each term is a variable name or
an interaction {it:var1}{cmd:#}{it:var2}, absorbed as the
cross-classification of its variables.

{phang}
{opt tolerance(#)} and {opt maxiter(#)} control the iterative convergence
criterion for {cmd:absorb()} and {cmd:fast}. Defaults: {cmd:tol(1e-8)},
{cmd:maxiter(100)}.


{marker fast}{...}
{title:Fast estimation: fast}

{pstd}
When datasets contain many individual spell records but the categorical
covariates define a limited number of risk cells (time period × education
× region, etc.), {cmd:fast} moves the heavy computation from the microdata
to the cells while returning the {bf:exact} individual-level MLE.

{pstd}
{bf:Term classification.} Each term in the varlist is classified by its
pieces:

{phang2}
• Terms built only from factor-variable pieces — {cmd:i.edu},
{cmd:i.edu#i.region} — are constant within cells and enter the model at
the cell level.{p_end}
{phang2}
• Terms involving at least one continuous piece — {cmd:age},
{cmd:c.age#c.income}, {cmd:c.age#i.edu} — enter at the individual
level.{p_end}

{pstd}
The cells are the cross-classification of the categorical variables in
the cell-level terms. At least one cell-level term is required.

{pstd}
{bf:Algorithm.} Newton iterations on the full likelihood. Each iteration
accumulates exponentially tilted within-cell moments of the
individual-level columns —

{p 12 12 2}
Σ{sub:i∈j} μ{sub:i},{space 3}Σ{sub:i∈j} μ{sub:i}X{sub:i},{space 3}Σ{sub:i} μ{sub:i}X{sub:i}X{sub:i}'

{pstd}
— in one O(n·p²) pass over the microdata, with p the number of
{it:individual-level} columns only; all remaining algebra is at the cell
level. When {it:every} term is cell-level, the iterations run entirely
on the collapsed (events, exposure) cell sums at O(J) per iteration.
Group structures are computed with a single sort and running-sum
segment aggregation; there are no data copies or collapse round-trips.

{pstd}
{bf:Exactness.} Coefficients, the log likelihood, and the OIM, robust,
and cluster–robust standard errors are numerically identical to the full
{cmd:poisson} model (verified to 1e-6 in the test suite).


{marker hdfe}{...}
{title:HDFE: absorb(terms)}

{pstd}
{cmd:absorb()} absorbs one or more sets of fixed effects (e.g., firm,
individual, region, or interactions such as {cmd:region#period}) using
iteratively reweighted alternating projections — the same mathematical
approach as
{browse "https://github.com/sergiocorreia/ppmlhdfe":ppmlhdfe}, implemented
in pure Mata with no external dependency.

{pstd}
{bf:Algorithm.} At each IRLS iteration:

{phang2}
1. Compute the linearized working variable {it:z_i} and weights {it:w_i = μ_i}.

{phang2}
2. Weighted-demean {it:z} and all covariates X by each absorbed term in
turn (Gauss–Seidel alternating projections), repeating until the inner
loop converges. Each demeaning pass is a single O(n) sweep using
precomputed sort permutations and group boundaries, independent of the
number of levels.

{phang2}
3. Run weighted least squares of the demeaned {it:z̃} on demeaned X̃ to
update β.

{phang2}
4. Recover the FE contributions and update {it:η}; repeat until
max|Δβ| / (1 + |β|) < {cmd:tolerance()}.

{pstd}
{bf:Standard errors.} {cmd:vce(oim)} (default), {cmd:vce(robust)}, and
{cmd:vce(cluster} {it:varname}{cmd:)} are all numerically exact relative
to the model with explicit dummy variables: OIM uses the partitioned
Hessian inverse; robust the HC1 sandwich; cluster the G/(G−1)-corrected
clustered sandwich with scores summed within clusters.

{pstd}
{bf:Separation.} Observations in FE groups with no events (Σd = 0) are
dropped before estimation with a message (those FE effects diverge to
−∞).


{marker postestimation}{...}
{title:Postestimation}

{pstd}
Standard postestimation commands available after {cmd:stpois}:
{helpb contrast}, {helpb estat vce}, {helpb estimates}, {helpb lincom},
{helpb lrtest}, {helpb margins}, {helpb nlcom}, {helpb predict},
{helpb predictnl}, {helpb test}, {helpb testparm}.

{pstd}
{ul:predict} after {cmd:stpois}

{pstd}
{cmd:predict} [{it:type}] {it:newvar} [{it:if}] [{it:in}] [, {it:statistic}]

{synoptset 16 tabbed}{...}
{synopthdr:statistic}
{synoptline}
{synopt:{opt xb}}linear predictor; the default{p_end}
{synopt:{opt n}}predicted number of events{p_end}
{synopt:{opt hazard}}predicted hazard rate: exp(x{sub:i}β){p_end}
{synopt:{opt hr}}synonym for {opt hazard}{p_end}
{synopt:{opt surv:ival}}predicted survival: exp(−exp(x{sub:i}β)×(t−t0)){p_end}
{synoptline}
{p2colreset}{...}

{pstd}
After {cmd:absorb()}, {helpb margins} operates on the non-absorbed
regressors.


{marker examples}{...}
{title:Examples}

{pstd}Setup{p_end}
{phang2}{cmd:. webuse stan3}{p_end}
{phang2}{cmd:. stset t1, failure(died) id(id)}{p_end}

{pstd}Standard estimation{p_end}
{phang2}{cmd:. stpois age posttran}{p_end}
{phang2}{cmd:. stpois age posttran, irr}{p_end}
{phang2}{cmd:. stpois age posttran, vce(robust)}{p_end}

{pstd}HDFE: absorb fixed effects{p_end}
{phang2}{cmd:. stpois age, absorb(surgery)}{p_end}
{phang2}{cmd:. stpois age, absorb(surgery) vce(robust)}{p_end}

{pstd}Fast: cell-accelerated exact MLE{p_end}
{phang2}{cmd:. stpois age i.surgery, fast}{p_end}

{pstd}Interactions{p_end}
{phang2}{cmd:. stpois c.age##i.posttran i.surgery, fast}{p_end}
{phang2}{cmd:. stpois age, absorb(posttran#surgery)}{p_end}

{pstd}With episode-split data{p_end}
{phang2}{cmd:. stsplit caltime, at(1970 1980 1990 2000) after(time=dob)}{p_end}
{phang2}{cmd:. stpois i.caltime age i.edu, irr}{p_end}
{phang2}{cmd:. stpois i.caltime age i.edu, fast}{p_end}


{marker references}{...}
{title:References}

{phang}
Correia, S., P. Guimarães, and T. Zylkin. 2020.
Fast Poisson estimation with high-dimensional fixed effects.
{it:Stata Journal} 20: 95–115.

{phang}
Holford, T. R. 1980. The analysis of rates and survivorship using log-linear
models. {it:Biometrics} 36: 299–305.

{phang}
Laird, N., and D. Olivier. 1981. Covariance analysis of censored survival
data using log-linear analysis techniques. {it:Journal of the American}
{it:Statistical Association} 76: 231–240.

{phang}
Rabe-Hesketh, S., and A. Skrondal. 2022. {it:Multilevel and Longitudinal}
{it:Modeling Using Stata}, 4th ed. College Station, TX: Stata Press.


{title:Author}

{pstd}
Andreas Ljungström{break}
Swedish Institute for Social Research (SOFI), Stockholm University{break}
andreas.ljungstrom@sofi.su.se


{title:Also see}

{psee}
{helpb streg}, {helpb stcox}, {helpb stsplit}, {helpb poisson},
{browse "https://github.com/sergiocorreia/ppmlhdfe":ppmlhdfe}
{p_end}
