{smcl}
{* *! version 0.2.0  13jul2026}{...}
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
    {opt cl:uster} {it:clustvar}, {opt boot:strap}, {opt jack:knife};
    all three exact with {cmd:absorb()}{p_end}
{synopt:{opt robust}}synonym for {cmd:vce(robust)}{p_end}
{synopt:{opth cluster(varname)}}synonym for {cmd:vce(cluster} {it:varname}{cmd:)}; supported with {cmd:absorb()}{p_end}

{syntab:Fast estimation on collapsed cells}
{synopt:{opt fast(method)}}{it:method} is {opt offset} or {opt moments};
    see {help stpois##fast:Fast estimation}{p_end}
{synopt:{opt mt:opel}}with {cmd:fast(offset)}: Murphy–Topel-corrected SEs
    that propagate first-stage uncertainty{p_end}

{syntab:High-dimensional fixed effects}
{synopt:{opth absorb(varlist)}}absorb fixed effects via IRLS demeaning;
    see {help stpois##hdfe:HDFE}{p_end}
{synopt:{opt tol:erance(#)}}IRLS convergence tolerance; default {cmd:1e-8}{p_end}
{synopt:{opt maxiter(#)}}maximum IRLS iterations; default {cmd:100}{p_end}

{synopt:{it:poisson_options}}other options passed to {helpb poisson} (standard path){p_end}
{synoptline}
{p2colreset}{...}

{p 4 6 2}
{opt fweight}s, {opt iweight}s, and {opt pweight}s are allowed;
see {help weight}.{p_end}

{p 4 6 2}
Data must be {cmd:stset} before using {cmd:stpois}; see {helpb stset}.{p_end}

{p 4 6 2}
{cmd:fast()} and {cmd:absorb()} may not be combined.{p_end}


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
exp(x{sub:i}β).

{pstd}
{bf:Three estimation paths} are available:

{phang2}
1. {bf:Standard} (default). Wraps Stata's {cmd:poisson} with automatic exposure.
Identical coefficients to {cmd:streg, distribution(exponential)}.

{phang2}
2. {bf:fast(offset) and fast(moments)}. Fast paths that exploit the cell
structure of categorical covariates in large individual-level datasets.
{cmd:fast(moments)} is {bf:exact} (identical coefficients and SEs to the
full model); {cmd:fast(offset)} is exact for coefficients but conditions
SEs on a first stage. See {help stpois##fast:Fast estimation} below.

{phang2}
3. {bf:absorb()}. High-dimensional fixed-effect absorption via iteratively
reweighted alternating projections.
See {help stpois##hdfe:HDFE} below.


{marker options}{...}
{title:Options}

{phang}
{opt irr} reports coefficients as incidence-rate ratios (hazard ratios).
Affects display only.

{phang}
{opt nolog} suppresses iteration logs from sub-routines.

{phang}
{opt level(#)} sets the confidence level as a percentage. Default 95.

{phang}
{opt noconstant} suppresses the constant term.

{phang}
{opth vce(vcetype)} specifies standard error type. For {cmd:absorb()},
only {opt oim} and {opt robust} are supported; other types require
{browse "https://github.com/sergiocorreia/ppmlhdfe":ppmlhdfe}.

{phang}
{opt robust}, {opth cluster(varname)} are convenient synonyms.

{phang}
{opt tolerance(#)} and {opt maxiter(#)} control IRLS convergence for
{cmd:absorb()}. Defaults: {cmd:tol(1e-8)}, {cmd:maxiter(100)}.


{marker fast}{...}
{title:Fast estimation: fast(offset) and fast(moments)}

{pstd}
When datasets contain millions of individual spell records but a small
number of categorical risk cells (time period × education × region, etc.),
the fast paths move the heavy computation from the microdata to the
cells, achieving large speed gains. Both methods handle continuous
covariates.

{pstd}
{bf:Syntax convention.} {cmd:stpois} uses factor-variable notation to
distinguish covariate types in the fast paths:

{phang2}
• Bare variable names (e.g., {cmd:age}, {cmd:income}) → {bf:continuous}
  (moments computed or absorbed into exposure){p_end}
{phang2}
• Factor notation ({cmd:i.edu}, {cmd:i.region}) → {bf:categorical}
  (cells collapsed by these grouping variables){p_end}

{pstd}
{bf:fast(offset) — two-stage multiplicative offset}

{phang2}
{ul:Stage 1.} Estimate continuous-covariate coefficients on the full
individual-level dataset (categorical indicators included to avoid OVB):

{pmore2}
poisson _d {it:cont_vars cat_vars}, exposure(ptime)

{phang2}
{ul:Stage 2.} Generate a personalized exposure multiplier, collapse
to cells by the categorical variables, and re-estimate:

{pmore2}
structural_weight = exp(γ̂₁x₁ + γ̂₂x₂ + ...){break}
mutated_ptime = ptime × structural_weight{break}
collapse (sum) _d (sum) mutated_ptime, by({it:cat_vars}){break}
poisson _d {it:cat_terms}, exposure(mutated_ptime)

{phang2}
{ul:Mathematical basis.} The multiplicative property of the log link,
exp(γX + αD) = exp(γX) × exp(αD), allows the continuous effect to be
folded into the exposure. Within a categorical cell j:

{p 12 12 2}
Ẽ_j = Σ_{i∈j} ptime_i × exp(γ'X_i)

{phang2}
When this aggregated mutated exposure enters the Poisson model as an
offset, it correctly tracks the nonlinear shift from the within-cell
distribution of the continuous covariates.

{pstd}
{bf:fast(moments) — iterated cell moments (exact MLE)}

{phang2}
Fits the {bf:exact} individual-level maximum-likelihood model by Newton
iterations that separate the computation into a cheap microdata pass and
cell-level algebra. Each iteration accumulates exponentially tilted
within-cell moments of the continuous covariates —

{p 12 12 2}
Σ_{i∈j} μ_i,{space 3}Σ_{i∈j} μ_i X_i,{space 3}Σ_{i∈j} μ_i X_i X_i'
{space 3}with μ_i = t_i exp(γ'X_i + δ'w_j)

{phang2}
— and updates all coefficients (continuous γ and categorical δ jointly)
from the resulting score and Hessian. The per-iteration cost is
O(n·p²) with p the number of {it:continuous} covariates, versus
O(n·(p+k)²) for {cmd:poisson} with k categorical parameters, so the
speed advantage grows with the number of categorical terms.

{phang2}
Coefficients, log likelihood, and standard errors are numerically
identical to the full {cmd:poisson} model (verified to 1e-6 in the test
suite). {cmd:vce(oim)} (default), {cmd:vce(robust)}, and
{cmd:vce(cluster} {it:varname}{cmd:)} are all exact.
{opt tolerance(#)} and {opt maxiter(#)} control the Newton iterations.

{phang2}
{opt skewness} is obsolete (retained for backward compatibility; ignored
with a note): the estimator is exact, so no higher-order correction is
needed.

{pstd}
{bf:Standard errors under fast(offset)}

{phang2}
• By default, {cmd:fast(offset)} standard errors are conditional on the
  stage-1 estimates of the continuous coefficients and do not account
  for first-stage estimation uncertainty. The constant is most affected;
  cell-contrast SEs much less so.

{phang2}
• {opt mtopel} requests Murphy–Topel-corrected standard errors
  (Murphy and Topel 1985). Because stage 1 estimates the full model
  (continuous and categorical terms jointly), the corrected covariance
  of the second-stage parameters reduces to the corresponding block of
  the stage-1 joint VCE, which {cmd:stpois} extracts exactly: with
  {opt mtopel}, {cmd:fast(offset)} SEs equal the full joint-model SEs
  to numerical precision.

{phang2}
• {cmd:fast(moments)} carries no such caveat — it is exact — and is the
  recommended fast path. {cmd:fast(offset)} remains useful when the
  collapsed cell file itself is the object of interest (for example,
  when only aggregated data may leave a secure server).


{marker hdfe}{...}
{title:HDFE: absorb(varlist)}

{pstd}
{cmd:absorb(varlist)} absorbs one or more sets of fixed effects (e.g., firm,
individual, region) using iteratively reweighted alternating projections —
the same mathematical approach as
{browse "https://github.com/sergiocorreia/ppmlhdfe":ppmlhdfe} but
implemented in pure Mata with no external dependency.

{pstd}
{bf:Algorithm.} At each IRLS iteration:

{phang2}
1. Compute the linearized working variable {it:z_i} and weights {it:w_i = μ_i}.

{phang2}
2. Weighted-demean {it:z} and all covariates X by each absorbed variable
in turn (Gauss–Seidel alternating projections), repeating until the inner
loop converges.

{phang2}
3. Run weighted least squares of the demeaned {it:z̃} on demeaned X̃ to
update the coefficient vector β.

{phang2}
4. Recover the FE contributions via the complement of the FE projector
and update {it:η = Xβ + FE_contribution}.

{phang2}
5. Repeat until max|Δβ| / (1 + |β|) < {cmd:tolerance()}.

{pstd}
{bf:Standard errors.} The VCE is computed from the partitioned Hessian at
convergence: {it:V = (X̃'WX̃)^{−1}} (OIM) or the sandwich form (robust).
This is asymptotically equivalent to the VCE from the full Poisson model
with FE dummies. A single-FE acid test: {cmd:stpois x, absorb(g)} should
match {cmd:poisson _d x i.g, exposure(e)} on both coefficients and SEs.

{pstd}
{bf:Separation.} Observations in FE groups with no events (Σd = 0) are
automatically dropped before estimation (those FE effects diverge to
−∞). A warning is printed.

{pstd}
{bf:Performance.} The Mata implementation uses {it:O(n × n_groups)} per
inner demeaning pass. This is efficient for typical EHA use cases with
moderate numbers of FE levels (period, education, region). For
{bf:very} high-dimensional FEs (e.g., individual fixed effects with millions
of levels), {browse "https://github.com/sergiocorreia/ppmlhdfe":ppmlhdfe}
will be substantially faster due to its optimized HDFE solver.

{pstd}
{bf:VCE support.} {cmd:vce(oim)} (default), {cmd:vce(robust)}, and
{cmd:vce(cluster} {it:varname}{cmd:)} are all supported with {cmd:absorb()}.
All three are {bf:numerically exact} relative to the equivalent model
with explicit dummy variables: OIM uses the partitioned Hessian inverse;
robust uses the HC1 sandwich (n/(n-1) corrected meat); cluster uses the
G/(G-1) corrected clustered sandwich with scores summed within clusters.
Differences from {cmd:poisson ... i.group} are at the numerical precision
level (~1e-9).


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
After {cmd:absorb()} or {cmd:fast()}, postestimation commands that require
the full likelihood (e.g., {helpb lrtest}) are less meaningful.
{helpb margins} operates on the non-absorbed regressors.


{marker examples}{...}
{title:Examples}

{pstd}Setup{p_end}
{phang2}{cmd:. webuse stan3}{p_end}
{phang2}{cmd:. stset t1, failure(died) id(id)}{p_end}

{pstd}Standard estimation{p_end}
{phang2}{cmd:. stpois age posttran}{p_end}
{phang2}{cmd:. stpois age posttran, irr}{p_end}
{phang2}{cmd:. stpois age posttran, vce(robust)}{p_end}

{pstd}HDFE: absorb transplant center{p_end}
{phang2}{cmd:. stpois age, absorb(surgery)}{p_end}
{phang2}{cmd:. stpois age, absorb(surgery) vce(robust)}{p_end}

{pstd}Fast: two-stage offset (age is continuous; surgery is categorical){p_end}
{phang2}{cmd:. stpois age i.surgery, fast(offset)}{p_end}

{pstd}Fast: moment correction{p_end}
{phang2}{cmd:. stpois age i.surgery, fast(moments)}{p_end}
{phang2}{cmd:. stpois age i.surgery, fast(moments) skewness}{p_end}

{pstd}With episode-split data{p_end}
{phang2}{cmd:. stsplit caltime, at(1970 1980 1990 2000) after(time=dob)}{p_end}
{phang2}{cmd:. stpois i.caltime age i.edu, irr}{p_end}
{phang2}{cmd:. stpois i.caltime age i.edu, fast(moments)}{p_end}

{pstd}Compare standard vs HDFE vs fast{p_end}
{phang2}{cmd:. stpois age posttran}{p_end}
{phang2}{cmd:. estimates store full}{p_end}
{phang2}{cmd:. stpois age, absorb(posttran)}{p_end}
{phang2}{cmd:. estimates store hdfe}{p_end}
{phang2}{cmd:. stpois age i.posttran, fast(moments)}{p_end}
{phang2}{cmd:. estimates store fast_moments}{p_end}
{phang2}{cmd:. estimates table full hdfe fast_moments}{p_end}


{marker references}{...}
{title:References}

{phang}
Murphy, K. M., and R. H. Topel. 1985. Estimation and inference in two-step
econometric models. {it:Journal of Business & Economic Statistics} 3: 88–97.

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
{browse "https://github.com/sergiocorreia/ppmlhdfe":ppmlhdfe (for very high-dimensional FEs)}
{p_end}
