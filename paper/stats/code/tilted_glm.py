"""
tilted_glm.py -- Reference implementation of the tilted-cell-moment framework.

Accompanies the manuscript

    "Tilted Cell Moments: A One-Pass Exact Newton Method for
     Generalized Linear Models at Scale"

and its online supplement (absorption convergence theory, rates,
high-dimensional asymptotics, Firth/penalized/conditional extensions).

The module implements, for exponential dispersion families with canonical
or non-canonical log-concave links:

  * `fit_tilted`   -- exact Newton / Fisher-scoring iterations in which every
                      score and information block that touches the cell-level
                      design is accumulated from exponentially tilted
                      within-cell moments (Theorem 1 of the paper), in a
                      single O(N p_c) / O(N p_c^2) pass over the microdata;
                      optional absorption of additional high-dimensional
                      fixed-effect factors by nonlinear block Gauss--Seidel
                      (Algorithm S1; convergence theory in Supplement
                      Theorems S1-S3);
  * `fit_full`     -- textbook dense Newton on the full design, used to
                      verify that the tilted iterates coincide with the full
                      ones to machine precision (Section 5.1);
  * `profile_poisson` -- the concentrated (profile) Poisson likelihood for
                      one absorbed cell factor, whose score and Hessian are
                      Schur complements of tilted moments; this is the
                      estimator studied in the J/N -> c regime (Supplement
                      Theorem S6);
  * Firth's bias-reducing adjusted score integrated into the tilted
                      accumulation (Supplement Proposition S2);
  * observed-information, HC1-robust, and cluster-robust covariance
                      estimators computed cell-wise.

Only numpy (and scipy.linalg for symmetric solves) is required.  The code
favours clarity over micro-optimisation; all cell accumulations use
`np.bincount`, so every pass is O(N) regardless of the number of cells.
"""

from __future__ import annotations

import numpy as np
from dataclasses import dataclass, field
from scipy.linalg import cho_factor, cho_solve


# ----------------------------------------------------------------------
# Exponential dispersion families
# ----------------------------------------------------------------------
#
# Density  f(y; theta, phi) = exp{ (y*theta - b(theta))/phi + c(y, phi) }.
# With linear predictor eta = a'xi + offset and link g, mu = g^{-1}(eta).
# For the canonical link theta = eta, Newton weights are b''(eta) and
# observed and expected information coincide.  For a non-canonical link
# (gamma/log below) we use Fisher scoring with the standard IRLS weight
# omega = 1 / (V(mu) g'(mu)^2), which for the log link is mu^2/V(mu).


class Family:
    """Base class: canonical-link families override eta-parametrised maps."""

    name = "family"
    canonical = True

    def mean(self, eta):            # b'(eta) for canonical links
        raise NotImplementedError

    def weight(self, eta):          # Newton/Fisher weight per observation
        raise NotImplementedError

    def loglik(self, y, eta, scale=1.0):
        raise NotImplementedError

    def firth_ratio(self, eta):
        """b'''(eta)/(2 b''(eta)) -- the canonical-link Firth adjustment."""
        raise NotImplementedError


class Poisson(Family):
    """Poisson with log link (canonical).  `exposure` enters via offset."""

    name = "poisson"

    def mean(self, eta):
        return np.exp(eta)

    def weight(self, eta):
        return np.exp(eta)

    def loglik(self, y, eta, scale=1.0):
        # up to the additive constant -sum(log y!)
        return float(np.sum(y * eta - np.exp(eta)))

    def firth_ratio(self, eta):
        return 0.5  # b''' = b'' = exp(eta)


class Binomial(Family):
    """Bernoulli/binomial with logit link (canonical)."""

    name = "binomial"

    def mean(self, eta):
        return 1.0 / (1.0 + np.exp(-eta))

    def weight(self, eta):
        p = self.mean(eta)
        return p * (1.0 - p)

    def loglik(self, y, eta, scale=1.0):
        # y in [0,1] (proportions allowed with prior weights folded into y)
        return float(np.sum(y * eta - np.logaddexp(0.0, eta)))

    def firth_ratio(self, eta):
        p = self.mean(eta)
        return 0.5 * (1.0 - 2.0 * p)  # b''' = mu(1-mu)(1-2mu)


class GammaLog(Family):
    """Gamma with log link (non-canonical): Fisher scoring, weight == 1.

    V(mu) = mu^2 and g'(mu) = 1/mu give omega = 1 identically; the tilted
    moments therefore reduce to untilted cell moments for this family, an
    observation made in Section 3 of the paper.
    """

    name = "gamma"
    canonical = False

    def mean(self, eta):
        return np.exp(eta)

    def weight(self, eta):
        return np.ones_like(eta)

    def loglik(self, y, eta, scale=1.0):
        # log-likelihood kernel for fixed shape nu = 1/scale:
        #   -nu*(y/mu + log mu), up to terms in y and nu only
        nu = 1.0 / scale
        mu = np.exp(eta)
        return float(-nu * np.sum(y / mu + eta))

    def working_residual(self, y, eta):
        mu = np.exp(eta)
        return (y - mu) / mu  # (y-mu) * g'(mu)

    def score_residual(self, y, eta, scale=1.0):
        mu = np.exp(eta)
        return (y - mu) / mu / scale


FAMILIES = {"poisson": Poisson(), "binomial": Binomial(), "gamma": GammaLog()}


# ----------------------------------------------------------------------
# Design containers
# ----------------------------------------------------------------------


@dataclass
class CellDesign:
    """Model design split into cell-level and individual-level parts.

    Attributes
    ----------
    cells : (N,) int array of cell indices in 0..J-1
    W     : (J, k) cell-level design (rows indexed by cell id); typically
            cell dummies or a saturated categorical design plus intercept
    X     : (N, p) individual-level (continuous-involving) design; may have
            p = 0 columns
    offset: (N,) offset (log exposure for Poisson event-history data)
    """

    cells: np.ndarray
    W: np.ndarray
    X: np.ndarray
    offset: np.ndarray

    def __post_init__(self):
        self.N = self.cells.shape[0]
        self.J, self.k = self.W.shape
        self.p = self.X.shape[1] if self.X.ndim == 2 else 0

    def eta(self, gamma, delta, alpha_contrib=None):
        eta = self.W @ delta
        eta = eta[self.cells] + self.offset
        if self.p:
            eta = eta + self.X @ gamma
        if alpha_contrib is not None:
            eta = eta + alpha_contrib
        return eta


@dataclass
class FitResult:
    gamma: np.ndarray
    delta: np.ndarray
    loglik: float
    iterations: int
    converged: bool
    V: np.ndarray | None = None          # (p+k) x (p+k), gamma block first
    se: np.ndarray | None = None
    alpha: list = field(default_factory=list)   # absorbed effects per factor
    trace: list = field(default_factory=list)   # per-iteration diagnostics
    H: np.ndarray | None = None


# ----------------------------------------------------------------------
# Tilted-moment accumulation (the single O(N p^2) pass of Theorem 1)
# ----------------------------------------------------------------------


def _tilted_moments(y, eta, design: CellDesign, family: Family,
                    firth=False, Hinv_blocks=None):
    """One pass over the microdata: residual sums, tilted cell moments.

    Returns
    -------
    R  : (J,) cell sums of raw score residuals  sum_{i in j} (y_i - mu_i)
         (working-weighted for non-canonical links)
    M0 : (J,) tilted zeroth moments  sum_{i in j} omega_i
    M1 : (J, p) tilted first moments sum_{i in j} omega_i x_i
    Q  : (p, p) tilted Gram matrix   sum_i omega_i x_i x_i'
    Ux : (p,) individual-level score block
    """
    cells, X, J, p = design.cells, design.X, design.J, design.p
    omega = family.weight(eta)
    if family.canonical:
        resid = y - family.mean(eta)
    else:
        # Fisher scoring: score residual = (y - mu) g'(mu) * omega
        resid = family.score_residual(y, eta) * 1.0
        resid = omega * family.working_residual(y, eta)

    if firth:
        # Canonical-link Firth adjustment: residual -> residual + h * ratio,
        # with hat values h_i = omega_i a_i' H^{-1} a_i accumulated using the
        # partitioned inverse supplied in Hinv_blocks (see Proposition 4).
        h = _hat_values(design, omega, Hinv_blocks)
        resid = resid + h * family.firth_ratio(eta)

    R = np.bincount(cells, weights=resid, minlength=J)
    M0 = np.bincount(cells, weights=omega, minlength=J)
    if p:
        M1 = np.empty((J, p))
        for l in range(p):
            M1[:, l] = np.bincount(cells, weights=omega * X[:, l], minlength=J)
        Xw = X * omega[:, None]
        Q = X.T @ Xw
        Ux = X.T @ resid
    else:
        M1 = np.zeros((J, 0))
        Q = np.zeros((0, 0))
        Ux = np.zeros(0)
    return R, M0, M1, Q, Ux


def _hat_values(design: CellDesign, omega, Hinv_blocks):
    """Leverages h_i = omega_i a_i' H^{-1} a_i without forming the N x (p+k)
    design: a_i = (x_i', w_{j(i)}')', so h_i expands into three bilinear
    pieces, each computable with row-indexed cell quantities."""
    Vgg, Vgd, Vdd = Hinv_blocks       # blocks of H^{-1}: gamma, cross, delta
    X, W, cells, p = design.X, design.W, design.cells, design.p
    WV = W @ Vdd @ W.T if design.k else None
    # h = omega * ( x'Vgg x + 2 x'Vgd w + w'Vdd w )
    h = np.zeros(design.N)
    if p:
        h += np.einsum("ij,jk,ik->i", X, Vgg, X)
        Wrow = W[cells]                       # (N, k)
        h += 2.0 * np.einsum("ij,jk,ik->i", X, Vgd, Wrow)
        h += np.einsum("ij,jk,ik->i", Wrow, Vdd, Wrow)
    else:
        Wrow = W[cells]
        h += np.einsum("ij,jk,ik->i", Wrow, Vdd, Wrow)
    return omega * h


def _assemble(R, M0, M1, Q, Ux, W):
    """Assemble the (p+k) score and information from cell quantities."""
    k = W.shape[1]
    Ud = W.T @ R
    Hdd = (W * M0[:, None]).T @ W
    Hgd = M1.T @ W                        # (p, k)
    U = np.concatenate([Ux, Ud])
    p = Q.shape[0]
    H = np.zeros((p + k, p + k))
    H[:p, :p] = Q
    H[:p, p:] = Hgd
    H[p:, :p] = Hgd.T
    H[p:, p:] = Hdd
    return U, H


# ----------------------------------------------------------------------
# Absorbed factors: closed-form profile updates (Poisson) and
# weighted-demeaning alternating projections (all families)
# ----------------------------------------------------------------------


def _block_update(fam, y, eta_wo_g, groups, J_g, alpha0=None,
                  inner_newton=40, tol=1e-12):
    """Exact block maximiser for one absorbed factor given the linear
    predictor excluding that factor's effect.  For Poisson (canonical
    log) the maximiser is closed form; for any other canonical family it
    solves the per-group score equation sum_{i in g} b'(alpha_g + rest_i)
    = sum_{i in g} y_i by vectorised one-dimensional Newton (all groups
    at once).  This is the exact blockwise maximisation whose monotone
    ascent underlies Supplement Theorem S1 for log-concave links."""
    D = np.bincount(groups, weights=y, minlength=J_g)
    if fam.name == "poisson":
        den = np.bincount(groups, weights=np.exp(eta_wo_g), minlength=J_g)
        with np.errstate(divide="ignore"):
            return np.log(D) - np.log(den)
    a = np.zeros(J_g) if alpha0 is None else alpha0.copy()
    for _ in range(inner_newton):
        eta = a[groups] + eta_wo_g
        mu = fam.mean(eta)
        w = fam.weight(eta)
        f = np.bincount(groups, weights=mu, minlength=J_g) - D
        fp = np.bincount(groups, weights=w, minlength=J_g)
        step = f / np.where(fp > 1e-300, fp, 1e-300)
        a = a - step
        if np.max(np.abs(step)) < tol:
            break
    return a


def _solve_absorbed(fam, y, base_eta, absorb, alphas, tol=1e-10,
                    maxiter=1000):
    """Nonlinear block Gauss--Seidel over the absorbed factors at fixed
    (gamma, delta): Algorithm 1, inner loop.  Returns updated alphas and
    the total absorbed contribution.  Monotone in the likelihood by
    Supplement Theorem S1 (Poisson closed form; general canonical family
    by exact one-dimensional block maximisation)."""
    G = len(absorb)
    contrib = np.zeros_like(base_eta)
    for g in range(G):
        contrib += alphas[g][absorb[g]]
    for it in range(maxiter):
        max_move = 0.0
        for g in range(G):
            eta_wo_g = base_eta + contrib - alphas[g][absorb[g]]
            new = _block_update(fam, y, eta_wo_g, absorb[g],
                                len(alphas[g]), alpha0=alphas[g])
            move = np.max(np.abs(new - alphas[g])) if len(new) else 0.0
            max_move = max(max_move, move)
            contrib += (new - alphas[g])[absorb[g]]
            alphas[g] = new
        if max_move < tol:
            break
    return alphas, contrib, it + 1


# ----------------------------------------------------------------------
# Main fitting routines
# ----------------------------------------------------------------------


def fit_tilted(y, design: CellDesign, family="poisson", absorb=None,
               tol=1e-9, maxiter=100, firth=False, vce="oim",
               cluster=None, trace_iterates=False, scale=1.0,
               inner_tol=None):
    """Newton / Fisher scoring on the exact likelihood via tilted moments.

    Parameters
    ----------
    absorb : optional list of (N,) int arrays of group indices; each factor
             is absorbed by the nonlinear Gauss--Seidel of Algorithm 1
             (Poisson closed form; other canonical families by exact
             1-D block maximisation, matching Supplement Theorem S1).
    firth  : use Firth's bias-reducing adjusted score (canonical links).
    vce    : 'oim', 'robust', or 'cluster' (requires `cluster` indices).
    """
    fam = FAMILIES[family] if isinstance(family, str) else family
    absorb = absorb or []
    if absorb and firth:
        raise NotImplementedError("Firth with absorbed factors "
                                  "not implemented")
    if absorb and not fam.canonical:
        raise NotImplementedError("absorption assumes a canonical link")
    alphas = [np.zeros(int(a.max()) + 1) for a in absorb]

    p, k = design.p, design.k
    gamma = np.zeros(p)
    delta = np.zeros(k)
    # crude but effective starting value: match the overall mean rate
    if k and fam.name == "poisson":
        mu0 = max(np.mean(y), 1e-8) / max(np.mean(np.exp(design.offset)), 1e-12)
        if np.allclose(design.W[:, 0], 1.0):
            delta[0] = np.log(mu0)

    contrib = np.zeros(design.N)
    if absorb:
        alphas, contrib, _ = _solve_absorbed(
            fam, y, design.eta(gamma, delta), absorb, alphas,
            tol=inner_tol or tol)
    eta = design.eta(gamma, delta, contrib if absorb else None)
    ll = fam.loglik(y, eta, scale)

    Hinv_blocks = None
    trace = []
    converged = False
    H = None
    for it in range(1, maxiter + 1):
        if firth:
            # evaluate the leverage adjustment at the current iterate: one
            # extra unadjusted pass gives H(theta_t), whose partitioned
            # inverse feeds the hat values of the adjusted score U*(theta_t)
            R, M0, M1, Q, Ux = _tilted_moments(y, eta, design, fam)
            _, H0 = _assemble(R, M0, M1, Q, Ux, design.W)
            Hinv_blocks = _partitioned_inverse(H0, p)
        R, M0, M1, Q, Ux = _tilted_moments(y, eta, design, fam,
                                           firth=firth,
                                           Hinv_blocks=Hinv_blocks)
        U, H = _assemble(R, M0, M1, Q, Ux, design.W)
        try:
            step = np.linalg.solve(H, U)
        except np.linalg.LinAlgError:
            # information singular (e.g. vanishing weights under
            # separation): take the minimum-norm ascent direction and let
            # the iteration run to the cap, flagging non-convergence
            step = np.linalg.lstsq(H, U, rcond=None)[0]
        # Step-halving.  For the ordinary MLE we require (near-)ascent of
        # the log likelihood.  For Firth, the adjusted-score iteration is a
        # quasi-Fisher-scoring fixed point (Kosmidis & Firth 2010): the
        # step direction is not the exact gradient of the penalised
        # objective (the leverage adjustment is evaluated at the previous
        # iterate), so we follow standard practice (logistf, brglm2) and
        # halve only to keep the penalised likelihood finite and
        # non-decreasing up to a small slack, declaring convergence on the
        # step size.
        lam, accepted = 1.0, False
        while lam >= 2.0 ** -16:
            g_new = gamma + lam * step[:p]
            d_new = delta + lam * step[p:]
            if absorb:
                alphas_new = [a.copy() for a in alphas]
                alphas_new, contrib_new, _ = _solve_absorbed(
                    fam, y, design.eta(g_new, d_new), absorb, alphas_new,
                    tol=inner_tol or tol)
            else:
                contrib_new = None
            eta_new = design.eta(g_new, d_new, contrib_new)
            ll_new = fam.loglik(y, eta_new, scale)
            if firth:
                ll_new += 0.5 * _logdet_information(y, eta_new, design, fam)
                if it == 1 and lam == 1.0:
                    ll = fam.loglik(y, eta, scale) \
                        + 0.5 * _logdet_information(y, eta, design, fam)
                slack = 1e-6 * (1.0 + abs(ll))
            else:
                slack = 1e-12 * abs(ll)
            if np.isfinite(ll_new) and ll_new > ll - slack:
                accepted = True
                break
            lam *= 0.5
        if not accepted:
            break
        relmove = np.max(np.abs(np.concatenate([g_new - gamma,
                                                d_new - delta]))
                         / (1.0 + np.abs(np.concatenate([g_new, d_new]))))
        gamma, delta, ll, eta = g_new, d_new, ll_new, eta_new
        if absorb:
            alphas, contrib = alphas_new, contrib_new
        if trace_iterates:
            trace.append(np.concatenate([gamma, delta]).copy())
        if relmove < tol:
            converged = True
            break

    V, se = _vce(y, eta, design, fam, H, vce, cluster)
    return FitResult(gamma=gamma, delta=delta, loglik=ll, iterations=it,
                     converged=converged, V=V, se=se, alpha=alphas,
                     trace=trace, H=H)


def _partitioned_inverse(H, p):
    try:
        Hi = np.linalg.inv(H)
    except np.linalg.LinAlgError:
        Hi = np.linalg.pinv(H)
    return Hi[:p, :p], Hi[:p, p:], Hi[p:, p:]


def _logdet_information(y, eta, design, fam):
    R, M0, M1, Q, Ux = _tilted_moments(y, eta, design, fam)
    _, H = _assemble(R, M0, M1, Q, Ux, design.W)
    sign, logdet = np.linalg.slogdet(H)
    return logdet


def _vce(y, eta, design, fam, H, vce, cluster):
    if H is None:
        return None, None
    p, k, cells, W, X = design.p, design.k, design.cells, design.W, design.X
    try:
        Hinv = np.linalg.inv(H)
    except np.linalg.LinAlgError:
        Hinv = np.linalg.pinv(H)
    if vce == "oim":
        V = Hinv
    else:
        if fam.canonical:
            resid = y - fam.mean(eta)
        else:
            resid = fam.weight(eta) * fam.working_residual(y, eta)
        if vce == "robust":
            groups = np.arange(design.N)
        else:
            groups = np.asarray(cluster)
        Gn = int(groups.max()) + 1
        # per-cluster score sums, accumulated columnwise via bincount
        S = np.zeros((Gn, p + k))
        for l in range(p):
            S[:, l] = np.bincount(groups, weights=resid * X[:, l],
                                  minlength=Gn)
        Wrow = W[cells]
        for l in range(k):
            S[:, p + l] = np.bincount(groups, weights=resid * Wrow[:, l],
                                      minlength=Gn)
        meat = S.T @ S
        if vce == "robust":
            meat *= design.N / (design.N - 1.0)
        else:
            meat *= Gn / (Gn - 1.0)
        V = Hinv @ meat @ Hinv
    return V, np.sqrt(np.diag(V))


def fit_full(y, A, offset, family="poisson", tol=1e-9, maxiter=100,
             start=None, trace_iterates=False):
    """Dense textbook Newton on the full design A (N x (p+k)); used only to
    verify Theorem 1 numerically and to time the baseline."""
    fam = FAMILIES[family] if isinstance(family, str) else family
    n, q = A.shape
    theta = np.zeros(q) if start is None else start.copy()
    eta = A @ theta + offset
    ll = fam.loglik(y, eta)
    trace = []
    converged = False
    for it in range(1, maxiter + 1):
        omega = fam.weight(eta)
        if fam.canonical:
            resid = y - fam.mean(eta)
        else:
            resid = omega * fam.working_residual(y, eta)
        U = A.T @ resid
        H = (A * omega[:, None]).T @ A
        try:
            step = np.linalg.solve(H, U)
        except np.linalg.LinAlgError:
            step = np.linalg.lstsq(H, U, rcond=None)[0]
        lam, accepted = 1.0, False
        while lam >= 2.0 ** -16:
            theta_new = theta + lam * step
            eta_new = A @ theta_new + offset
            ll_new = fam.loglik(y, eta_new)
            if ll_new > ll - 1e-12 * abs(ll):
                accepted = True
                break
            lam *= 0.5
        if not accepted:
            break
        relmove = np.max(np.abs(theta_new - theta)
                         / (1.0 + np.abs(theta_new)))
        theta, eta, ll = theta_new, eta_new, ll_new
        if trace_iterates:
            trace.append(theta.copy())
        if relmove < tol:
            converged = True
            break
    return FitResult(gamma=theta, delta=np.zeros(0), loglik=ll,
                     iterations=it, converged=converged, trace=trace)


# ----------------------------------------------------------------------
# Profile (concentrated) Poisson likelihood for one absorbed cell factor:
# the estimator of the J/N -> c regime (Supplement Theorem S6)
# ----------------------------------------------------------------------


def profile_poisson(y, X, cells, offset, tol=1e-10, maxiter=50):
    """Maximise the Poisson profile likelihood in gamma after concentrating
    out one fixed effect per cell.  Score and Hessian are the Schur
    complements of the tilted moments; identical to conditional multinomial
    ML within cells (Lemma 1 of the paper).

    Returns gamma-hat, profile Hessian (negative), cluster-robust and
    profile-information covariance matrices.
    """
    N, p = X.shape
    J = int(cells.max()) + 1
    D = np.bincount(cells, weights=y, minlength=J)
    keep_cells = D > 0                    # cells with no events drop out
    gamma = np.zeros(p)
    for it in range(1, maxiter + 1):
        s = offset + X @ gamma
        E = np.bincount(cells, weights=np.exp(s), minlength=J)
        # mu_i with alpha-hat plugged in: mu_i = D_j * exp(s_i) / E_j
        mu = np.exp(s) * (D / np.where(E > 0, E, 1.0))[cells]
        resid = y - mu
        U = X.T @ resid
        M0 = np.bincount(cells, weights=mu, minlength=J)      # == D_j
        M1 = np.empty((J, p))
        for l in range(p):
            M1[:, l] = np.bincount(cells, weights=mu * X[:, l], minlength=J)
        Q = X.T @ (X * mu[:, None])
        with np.errstate(divide="ignore", invalid="ignore"):
            scal = np.where(M0 > 0, 1.0 / M0, 0.0)
        Hp = Q - M1.T @ (M1 * scal[:, None])   # Schur complement
        step = np.linalg.solve(Hp, U)
        gamma_new = gamma + step
        if np.max(np.abs(step) / (1.0 + np.abs(gamma_new))) < tol:
            gamma = gamma_new
            break
        gamma = gamma_new
    # covariance estimators at the optimum
    s = offset + X @ gamma
    E = np.bincount(cells, weights=np.exp(s), minlength=J)
    mu = np.exp(s) * (D / np.where(E > 0, E, 1.0))[cells]
    resid = y - mu
    M0 = np.bincount(cells, weights=mu, minlength=J)
    M1 = np.empty((J, p))
    for l in range(p):
        M1[:, l] = np.bincount(cells, weights=mu * X[:, l], minlength=J)
    Q = X.T @ (X * mu[:, None])
    with np.errstate(divide="ignore", invalid="ignore"):
        scal = np.where(M0 > 0, 1.0 / M0, 0.0)
    Hp = Q - M1.T @ (M1 * scal[:, None])
    Hpi = np.linalg.inv(Hp)
    # cluster (cell)-robust: profile scores sum within cells
    Sc = np.empty((J, p))
    for l in range(p):
        Sc[:, l] = np.bincount(cells, weights=resid * X[:, l], minlength=J)
    # centre by the cell-level projection: profile score per cell is already
    # the within-cell centred score (sum of resid within cell = 0), so Sc
    # rows are the exact per-cell profile scores.
    Gn = int(keep_cells.sum())
    meat = Sc.T @ Sc * (Gn / max(Gn - 1.0, 1.0))
    Vcl = Hpi @ meat @ Hpi
    return gamma, Hp, Hpi, Vcl


# ----------------------------------------------------------------------
# Weighted alternating projections (linear demeaning) + spectral theory
# ----------------------------------------------------------------------


def demean_ap(v, factors, weights, tol=1e-12, maxiter=5000,
              return_history=False):
    """Weighted alternating projections: residual of projecting v onto the
    orthogonal complement of span{factor dummies} in the diag(weights)
    inner product.  Gauss--Seidel sweeps over the factors; the linear rate
    for two factors is cos^2 of the Friedrichs angle (Supplement
    Theorem S2)."""
    v = v.astype(float).copy()
    wsums = [np.bincount(f, weights=weights) for f in factors]
    history = []
    for it in range(maxiter):
        move = 0.0
        sweep_change = 0.0
        for g, f in enumerate(factors):
            m = np.bincount(f, weights=weights * v) / wsums[g]
            v -= m[f]
            move = max(move, np.max(np.abs(m)))
            sweep_change += np.sum(wsums[g] * m * m)
        if return_history:
            # norm of the sweep's total adjustment: decays geometrically at
            # the spectral radius of the Gauss--Seidel operator
            history.append(np.sqrt(sweep_change))
        if move < tol:
            break
    return (v, history, it + 1) if return_history else v


def friedrichs_rate(f1, f2, weights):
    """Predicted asymptotic linear rate of the two-factor weighted
    alternating projections: sigma_2^2, with sigma_2 the second-largest
    singular value of the weighted bipartite transition kernel
    D_1^{-1/2} B D_2^{-1/2}  (Supplement Theorem S2 / Corollary S1)."""
    J1, J2 = int(f1.max()) + 1, int(f2.max()) + 1
    B = np.zeros((J1, J2))
    np.add.at(B, (f1, f2), weights)
    d1 = B.sum(axis=1)
    d2 = B.sum(axis=0)
    T = B / np.sqrt(np.outer(d1, d2))
    svals = np.linalg.svd(T, compute_uv=False)
    # largest singular value is 1 (the constants); the Friedrichs cosine is
    # the second one restricted to the complement of the intersection
    sigma2 = svals[1] if len(svals) > 1 else 0.0
    return sigma2 ** 2


# ----------------------------------------------------------------------
# Utilities
# ----------------------------------------------------------------------


def make_cell_dummies(cells, J=None, drop_first=True, intercept=True):
    """Cell-level design W (J x k): intercept + J-1 dummies by default."""
    J = J or int(cells.max()) + 1
    if intercept:
        k = J if drop_first else J + 1
        W = np.zeros((J, k))
        W[:, 0] = 1.0
        for j in range(1 if drop_first else 0, J):
            W[j, j if drop_first else j + 1] = 1.0
    else:
        W = np.eye(J)
    return W


def flops_per_iteration(N, p, k):
    """Leading-order FLOP counts of one Newton iteration (Corollary 1)."""
    full = N * (p + k) ** 2
    tilted = N * p ** 2 + N * p + (0 if k == 0 else 0) + (p + k) ** 3 / 3
    return full, tilted


# ----------------------------------------------------------------------
# G >= 3 absorbed factors: exact local rate and Friedrichs product bound
# (Supplement Theorem S3)
# ----------------------------------------------------------------------


def _level_space_blocks(factors, weights):
    """Assemble the level-space information blocks of the absorbed
    factors: diagonal weight totals D_g and pairwise weighted
    cross-tabulations B[g][h] (J_g x J_h).  These pairwise tables
    determine F_alpha,alpha completely -- the basis of Supplement Theorem S3(i)."""
    G = len(factors)
    J = [int(f.max()) + 1 for f in factors]
    D = [np.bincount(f, weights=weights, minlength=J[g])
         for g, f in enumerate(factors)]
    B = [[None] * G for _ in range(G)]
    for g in range(G):
        for h in range(g + 1, G):
            Bgh = np.zeros((J[g], J[h]))
            np.add.at(Bgh, (factors[g], factors[h]), weights)
            B[g][h] = Bgh
            B[h][g] = Bgh.T
    return J, D, B


def gs_spectral_rate(factors, weights, tol=1e-10):
    """Exact asymptotic per-sweep rate of the nonlinear Gauss--Seidel /
    weighted alternating projections over G >= 2 absorbed factors
    (Supplement Theorem S3(i)): the spectral radius of the block Gauss--Seidel
    operator of F_alpha,alpha on the quotient by its null space.

    The operator lives on the level space (dimension sum_g J_g) and is
    built solely from the pairwise weighted cross-tabulations; the flat
    directions (shifts between factors) are exact eigenvalue-1 fixed
    points and are removed before taking the maximum modulus.
    """
    G = len(factors)
    J, D, B = _level_space_blocks(factors, weights)
    n = sum(J)
    off = np.cumsum([0] + J)
    H = np.zeros((n, n))
    for g in range(G):
        H[off[g]:off[g + 1], off[g]:off[g + 1]] = np.diag(D[g])
        for h in range(G):
            if h != g:
                H[off[g]:off[g + 1], off[h]:off[h + 1]] = B[g][h]
    L = np.tril(H)            # D + strictly-lower blocks (D_g diagonal)
    U = np.triu(H, 1)
    M = -np.linalg.solve(L, U)
    lam = np.linalg.eigvals(M)
    # flats: eigenvalue exactly 1 with multiplicity = dim null(H)
    lam = lam[np.abs(lam - 1.0) > 1e-8]
    return float(np.max(np.abs(lam))) if len(lam) else 0.0


def friedrichs_product_bound(factors, weights):
    """A priori upper bound on the per-sweep contraction (Supplement Theorem S3(ii)):
    the one-sweep operator norm of the demeaning cycle obeys the
    Smith--Solmon--Wagner / Deutsch--Hundal product bound

        rho_GS <= ||sweep||_Omega <= [1 - prod_{g=1}^{G-1}(1 - c_g^2)]^{1/2},

    where c_g is the Friedrichs cosine between V_g and
    V_{g+1} + ... + V_G in the Omega inner product.  Each c_g is
    computed on the level space via a generalized Rayleigh quotient,
    with the exact-intersection directions (eigenvalue 1) removed.
    For G = 2 the bound returns c_F itself; the sharp asymptotic rate
    is then c_F^2 (Kayalar--Weinert), and in general the exact rate is
    `gs_spectral_rate`, for which this is a conservative certificate.
    """
    G = len(factors)
    J, D, B = _level_space_blocks(factors, weights)
    bound_complement = 1.0
    for g in range(G - 1):
        rest = list(range(g + 1, G))
        nR = sum(J[h] for h in rest)
        offR = np.cumsum([0] + [J[h] for h in rest])
        HR = np.zeros((nR, nR))
        BgR = np.zeros((J[g], nR))
        for a, h in enumerate(rest):
            HR[offR[a]:offR[a + 1], offR[a]:offR[a + 1]] = np.diag(D[h])
            BgR[:, offR[a]:offR[a + 1]] = B[g][h]
            for b, k in enumerate(rest):
                if k != h:
                    HR[offR[a]:offR[a + 1], offR[b]:offR[b + 1]] = B[h][k]
        dg = np.sqrt(D[g])
        K = (BgR / dg[:, None]).T
        K = np.linalg.pinv(HR, rcond=1e-12) @ K
        K = (BgR / dg[:, None]) @ K          # D^-1/2 B H_R^+ B' D^-1/2
        ev = np.linalg.eigvalsh((K + K.T) / 2)
        ev = ev[ev < 1.0 - 1e-8]             # remove exact intersections
        cg2 = float(np.max(ev)) if len(ev) else 0.0
        cg2 = min(max(cg2, 0.0), 1.0)
        bound_complement *= (1.0 - cg2)
    return float(np.sqrt(1.0 - bound_complement))
