---
artifact: local_check
scope: local
theorem_ids: [thm:tilted, cor:complexity, prop:noncanonical, lem:levelsets, thm:global, thm:rate, cor:bipartite, thm:multiG, thm:classical, thm:diverging, lem:conditional, thm:proportional, prop:penalized, prop:firth]
generated: 2026-07-18
---

# Critical-path and support-lemma checks (units that passed)

Each unit below was checked with the full per-step template (claim normalization,
dependencies, step-by-step verification, step-completeness/skip audit, edge cases,
negligibility ledger where relevant). Only the verdicts and any residual notes are
recorded here; all are COMPLETE proofs.

---

## T1 `thm:tilted` — Exact tilted-moment representation
- Strategy: direct differentiation. Provability: PROVABLE AS STATED. **Verified.**
- Score: $\phi\partial\ell/\partial\gamma=\sum_i(y_i-b'(\eta_i))\bx_i$,
  $\phi\partial\ell/\partial\delta=\sum_i(y_i-b'(\eta_i))\bw_{j(i)}=\sum_j\bw_jR_j$ ✓
  (grouping uses $\bw_j$ constant within cell). Info $\phi\bF=\sum_i b''(\eta_i)\ba_i\ba_i'$;
  block expansion + Def-tilt moment identities ✓.
- (iii) Poisson: $b'=b''=\exp$ ⇒ $\omega_i=e^{\delta'\bw_j}t_ie^{\gamma'\bx_i}$, cell factor
  cancels in normalization, $M_j^{(0)}=e^{\delta'\bw_j}\mathcal L_{\nu_j}(\gamma)$;
  $\E_{P_{j,\gamma}}[\bx]=\nabla K_j$, $\Var=\nabla^2K_j$, $\E[\bx\bx']=\Var+\E\E'$ ✓.
- (iv) exactness: same matrices entry-for-entry ⇒ same Newton step ✓ (immediate).
- Note (S3): the "$b'$-tilt" $P^{b'}_{j,\theta}$ in (i) is a probability measure only when
  $b'>0$; true for the in-scope families ($b'=\mu>0$) but the positivity is left implicit.

## C1 `cor:complexity` — Scaling
- Strategy: operation counting. **Verified.** Per-obs: $\eta_i$ ($O(p)$), $b',b''$ ($O(1)$),
  binned $\omega_i,\omega_i\bx_i$ ($O(p)$), rank-one $\omega_i\bx_i\bx_i'$ ($O(p^2)$) ⇒
  $O(N(p^2{+}p{+}1))$; assembly $O(J(k^2{+}kp))$; solve $O((p{+}k)^3)$; dense
  $O(N(p{+}k)^2)$ ✓. Ratios follow. $p=0$ degeneracy ✓.
- Note (S3): supplement "Notation recalled" writes $O(N(p^2{+}p))$ (drops the $+1$); cosmetic.

## P2 `prop:noncanonical` — Tilted Fisher scoring
- Strategy: reduction to T1 with IRLS weight. **Verified.** ML score
  $\sum_i(y_i-\mu_i)\frac{\partial\mu/\partial\eta}{V(\mu)}\ba_i=\sum_i r_i\ba_i$ using
  $\partial\mu/\partial\eta=1/g'(\mu)$; expected info $\sum_i\omega_i\ba_i\ba_i'$ ✓.
  Fixed point $\sum_i r_i\ba_i=0$ = ML score root ✓. Gamma–log:
  $\omega_i=1/(V(\mu)g'(\mu)^2)=1/(\mu^2\cdot\mu^{-2})=1$ ✓.

## LS1 `lem:levelsets` — Compact superlevel sets
- Strategy: contradiction via unbounded ray (Rockafellar Thm 8.4). **Verified.**
  $\varphi(t)=f(x^*+t\bv)$ strictly concave, $\ge a$, max at $0$ ⇒ right-derivative
  $\varphi'(t^+)\ge0$ ∀t>0 (else $\to-\infty$); but strict concavity ⇒ $\varphi'(t^+)<\varphi'(0^+)\le0$
  for $t>0$ — contradiction ✓.

## S1 `thm:global` — Global convergence (= main Thm 2)
- Strategy: block-coordinate ascent on concave (strictly, on quotient $E=\mathcal K^\perp$)
  profile objective. Provability: PROVABLE AS STATED. **Verified.**
- Strict concavity on $E$: $u\in E\setminus0\Rightarrow$ some $\ba_i'u\ne0$ (else
  $u\in\mathcal K\cap\mathcal K^\perp=0$), $b''>0$ ✓.
- Monotonicity + confinement to compact $\bar S$ via LS1 ✓. Limit-point-stationary uses
  Berge maximum theorem (unique block maximizers) + Armijo sufficient ascent
  (Bertsekas Prop 1.2.1) for the β-block ✓; cyclic BCD is legitimate here because block
  maximizers are unique (Tseng 2001 conditions met — avoids the Powell non-convergence
  pathology). Stationary ⇒ global by concavity ✓.
- (iv) log-concave links: re-inspection of where A-canonical entered (concavity, finite
  block maximizers, sufficient ascent) confirms each survives under a strictly-concave
  per-obs contribution + bounded-eigenvalue symmetric preconditioner (Fisher scoring). ✓
- Note (S3): Algorithm S1 runs the inner α-loop **to tolerance** then a single damped
  Newton β-step; the proof analyzes essentially-cyclic single block updates. The equivalence
  is asserted via "every block updated at least once every B iterations" — correct but the
  inner-to-convergence vs one-sweep distinction could be stated more explicitly.

## S2 `thm:rate` — Local linear convergence
- Strategy: implicit-function-theorem Jacobian = block Gauss–Seidel operator. **Verified.**
- (i) $DT=-(\bD+\mathbf L)^{-1}\mathbf L'$, SPD $F_{\alpha\alpha}$ on quotient ⇒ $\rho_{GS}<1$
  (Golub–Van Loan 11.2) + Ostrowski (Ortega 10.1) ✓.
- (ii) two-block GS eigenvalues $\{0\}\cup\mathrm{eig}(F_{22}^{-1}F_{21}F_{11}^{-1}F_{12})$;
  in $\Omega^*$ inner product this is $P_2P_1|_{V_2}$ with spectrum $\{\cos^2\theta_i\}$ ⇒
  $\rho=c_F^2$ on quotient ✓. Demeaning-residual rate $c_F^{2n-1}$ (Aronszajn/Kayalar–Weinert) ✓.
- (iii) β-block appended, IFT + SPD ✓.

## CS1 `cor:bipartite` — Bipartite spectral gap
- **Verified.** $\langle u,v\rangle_{\Omega^*}=\tilde u'\bB\tilde v$; substitution
  $\mathbf r=\bD_1^{1/2}\tilde u$ turns Friedrichs sup into $\sigma_2$ of
  $\bT=\bD_1^{-1/2}\bB\bD_2^{-1/2}$; top pair $(\bD_1^{1/2}\mathbf1,\bD_2^{1/2}\mathbf1)$
  has $\sigma_1=1$ (checked: $\bT\bD_2^{1/2}\mathbf1=\bD_1^{1/2}\mathbf1$ since
  $\bB\mathbf1=\bD_1\mathbf1$) ✓; Perron–Frobenius ⇒ $\sigma_2<1$ iff connected ✓.

## S3 `thm:multiG` — Exact rate & certificates, $G\ge2$
- Strategy: spectral analysis of level-space sweep operator. Provability: AS STATED. **Verified.**
- (ii) **Semisimplicity of eigenvalue 1** checked in full: with $\mathbf M-\bI=-(\mathcal D+\mathcal L)^{-1}\mathcal F$,
  a Jordan block gives $\mathcal Fv=-(\mathcal D+\mathcal L)f$, $f\in\mathrm{null}(\mathcal F)$;
  left-multiply by $f'$: LHS $=(\mathcal Ff)'v=0$; RHS uses $\mathcal Ff=0\Rightarrow
  f'\mathcal Df+2f'\mathcal Lf=0\Rightarrow f'(\mathcal D+\mathcal L)f=\tfrac12f'\mathcal Df>0$
  ($\mathcal D\succ0$) ⇒ RHS $<0\ne0$ — contradiction ✓✓. Rate = max non-unit eigenvalue
  (semiconvergence, Keller 1965/Berman–Plemmons) ✓. Identification equivalence ✓. Cost ✓.
- (iii) product-of-sines: multiplicative subspace correction $\mathbf E=(\bI-\bP_G)\cdots(\bI-\bP_1)$,
  Smith (1977) bound + Friedrichs complement duality $c(A,B)=c(A^\perp,B^\perp)$ ⇒
  $(1-\prod(1-c_g^2))^{1/2}$ ✓; correctly noted conservative for $G\ge3$.
- (iv) $G=2$ reduction to $c_F^2$ / bound $c_F$ ✓.

## S4 `thm:classical` — Classical asymptotics
- Strategy: reduction to Fahrmeir–Kaufmann (1985) + Gouriéroux et al (1984). **Conditionally
  verified.** The delegation is appropriate for a supplement; A-classical's conditions
  ($\lambda_{\min}(F_N)\to\infty$, $N^{-1}F_N\to F_\infty\succ0$, bounded covariates) are the
  natural match to F–K's divergence/normalization conditions.
- Note (S3, issue I-04): the exact correspondence to F–K's regularity hypotheses is asserted,
  not re-derived; a referee may want the specific conditions named.

## S5 `thm:diverging` — Diverging dimension
- Strategy: quadratic-minorization in intrinsic norm (self-contained). Provability: AS STATED.
  **Verified.** All five steps re-derived:
  Step1 $\E\|\nabla g(0)\|^2=\mathrm{tr}(F^{-1/2}\Var(U)F^{-1/2})=q_N$ (info equality) ✓.
  Step2 $\|F^{-1/2}\ba_i\|^2=h_i/\omega_i\le Cq_N/(N\underline\omega)$ ⇒
  $\max_i|\eta_i-\eta_{0i}|\le C_rq_N/\sqrt N=\varepsilon_N\to0$ (needs $q_N^2/N\to0$) ✓;
  Hessian sandwich $\delta_N=L\varepsilon_N/\underline\omega$ ✓ (uses $u'Fu\ge\underline\omega\sum(\ba_i'u)^2$).
  Step3 existence+rate $\|\hat\bs\|=O_p(\sqrt{q_N})$ via Taylor on sphere ✓.
  Step4 linearity remainder $=O_p(q_N^{3/2}/\sqrt N)=o_p(1)$ **iff $q_N^3/N\to0$** — matches
  the assumption exactly ✓. Step5 contrast Lindeberg CLT, $\max|c_{N,i}|=O(\sqrt{q_N/N})\to0$ ✓.
- Careful, tight, honest use of every assumption. No issues.

## LS2 `lem:conditional` — Profile = conditional (Poisson)
- **Verified.** $e^{\hat\alpha_j}=D_j/E_j$ ($D_j>0$); substitution ⇒ multinomial log-lik +
  $\gamma$-free terms (Poisson-splitting) ✓; concavity (log-sum-exp) ✓; Schur form at
  $\hat\alpha(\gamma)$: $M_j^{(0)}=D_j$, $\mathbf m_j=D_j\E_{\pi_j}[\bx]$,
  $\bQ-\sum_j\mathbf m_j\mathbf m_j'/M_j^{(0)}=\sum_j D_j\Var_{\pi_j}$ ✓.

## S6 `thm:proportional` — Proportional regime
- Strategy: convexity-lemma consistency + Lyapunov CLT + efficiency projection. Provability:
  AS STATED. **Verified.** Headline "no asymptotic bias" rests on $\E[\bs_j\mid D_j]=0$
  (multinomial score conditionally mean-zero) ⇒ $\E\bs_j=0$ exactly ∀J ✓ — this is the
  genuine mechanism, correctly identified.
- (i) $\Var(L_J)\le C/J$ (bounded $D_j$ moments), per-cell KL separation, convexity lemma
  (Andersen 1982/Pollard 1991) ⇒ uniform convergence + argmax consistency ✓.
- (ii) Lyapunov CLT (bounded third moments), Hessian uniform convergence ✓.
- (iv) **Efficiency projection re-derived in full**: projection residual of $\gamma$-score onto
  cell nuisance score $u_j=\sum_i(y_{ij}-\mu_{ij})$ equals $\sum_i y_{ij}\bx_{ij}-D_j\bar\bx_j^\pi
  =\bs_j$ (using $\bar\bx_j^\mu=\bar\bx_j^\pi$ since $\pi_{ij}\propto\mu_{ij}$) ✓✓ — the
  conditional score is the efficient score. Bound attained (van der Vaart Ch. 25) ✓.

## PS1 `prop:penalized` — Penalized proximal Newton
- **Verified.** Ridge: $\lambda_2\bI_p$ to $\gamma$-block ✓; lasso coordinate cycle $O((p{+}k)^2)$
  data-free after the pass ✓ (Friedman 2010 inner solver); proximal-Newton local quadratic
  convergence (Lee–Sun–Saunders 2014) ✓; KKT/active set functions of $(\bU,\bF)$ ✓.

## PS2 `prop:firth` — Tilted Firth
- **Verified.** $A_r=\tfrac12\mathrm{tr}(F^{-1}\partial F/\partial\theta_r)
  =\sum_i h_i\frac{b'''}{2b''}a_{ir}$ (checked the trace algebra) ✓; leverage partition
  $\ba_i'F^{-1}\ba_i=\bx_i'\bV_{\gamma\gamma}\bx_i+2\bx_i'\bc_{j(i)}+d_{j(i)}$ with
  $\bc_j=\bV_{\gamma\delta}\bw_j$, $d_j=\bw_j'\bV_{\delta\delta}\bw_j$ ✓; complexity ✓;
  $b'''/2b''=1/2$ (Poisson), $(1-2\mu_i)/2$ (logit) ✓; finiteness (Heinze–Schemper 2002) ✓.
