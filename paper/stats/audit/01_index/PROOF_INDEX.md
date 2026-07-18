# Proof Index

- Script version: `1.0.0`
- Rules version: `1.0.0`
- Rules digest: `aded664a1907b349`
- Reference mode: `two-file`
- Supplement mode: `separate-self-contained`

## Summary
- Indexed units: 23
- Mechanical FAIL: 0
- Mechanical WARN: 0
- Topological layers: 1

## Proof Unit Inventory

| Type | Label | Line | Depends on | Summary |
|---|---|---|---|---|
| proposition | `prop:sufficiency` | 406 | — | [Failure of fixed collapse] Let the cumulant $b$ be real-analytic on $\R$ with $ |
| definition | `def:tilt` | 460 | — | [Tilted within-cell measures] Fix $\btheta = (\bgamma, \bdelta)$ and let $\omega |
| theorem | `thm:tilted` | 499 | — | [Exact tilted-moment representation] For the model \eqref{eq:edf}--\eqref{eq:lin |
| corollary | `cor:complexity` | 573 | — | [Scaling] One Newton iteration computed via Theorem~\ref{thm:tilted} costs \[ \u |
| proposition | `prop:noncanonical` | 625 | — | [Tilted Fisher scoring] With $\omega_i$ and $r_i$ as in \eqref{eq:fisher-weight} |
| theorem | `thm:global-main` | 864 | — | [Global convergence of the absorbed-factor iteration] Suppose the link is canoni |
| assumption | `ass:canonical` | 219 | — | The link is canonical and $b$ is strictly convex and twice continuously differen |
| assumption | `ass:mle-exists` | 224 | — | The MLE exists: $\sup_{\btheta} \ell$ is attained at some $\btheta^\ast$ (unique |
| assumption | `ass:blocks` | 229 | — | (i) Each absorbed group contains at least one observation with $y_i$ in the inte |
| theorem | `thm:global` | 247 | — | [Global convergence of nonlinear alternating projections; main-paper Theorem~2]  |
| theorem | `thm:rate` | 332 | — | [Local linear convergence] Let Assumptions~\ref{ass:canonical}--\ref{ass:blocks} |
| corollary | `cor:bipartite` | 360 | — | [Bipartite spectral gap] For $G = 2$, define the weighted biadjacency matrix $\b |
| theorem | `thm:multiG` | 406 | — | [Exact rate and certificates for $G \ge 2$ factors] Let the assumptions of Theor |
| assumption | `ass:classical` | 511 | — | (i) $q = p + k$ is fixed; covariates and offsets are uniformly bounded; (ii) $\l |
| theorem | `thm:classical` | 519 | — | [Classical asymptotics] Under Assumptions~\ref{ass:canonical} and \ref{ass:class |
| assumption | `ass:diverging` | 547 | — | (i) Covariates and offsets are uniformly bounded, and the true linear predictors |
| theorem | `thm:diverging` | 562 | — | [Diverging dimension] Under Assumptions~\ref{ass:canonical} and \ref{ass:divergi |
| lemma | `lem:conditional` | 610 | — | [Profile equals conditional] For every $\bgamma$ with all $D_j$-positive cells,  |
| assumption | `ass:proportional` | 647 | — | (i) Cell sizes satisfy $2 \le m_j \le \bar m$; $J/N \to c \in (0, 1)$; (ii) cova |
| theorem | `thm:proportional` | 656 | — | [Proportional regime] Under Assumption~\ref{ass:proportional}: (i) the profile M |
| proposition | `prop:penalized` | 766 | — | [Penalized tilted-moment proximal Newton] At each outer iteration, form the quad |
| proposition | `prop:firth` | 800 | — | [Tilted Firth] Given the partitioned inverse of $\bF$ (available from the Newton |
| lemma | `lem:levelsets` | 1090 | — | [Compact superlevel sets] Let $f : \R^d \to \R$ be strictly concave, continuous, |

## Suggested Check Order (topological layers)

- Layer 0: `ass:blocks`, `ass:canonical`, `ass:classical`, `ass:diverging`, `ass:mle-exists`, `ass:proportional`, `cor:bipartite`, `cor:complexity`, `def:tilt`, `lem:conditional`, `lem:levelsets`, `prop:firth`, `prop:noncanonical`, `prop:penalized`, `prop:sufficiency`, `thm:classical`, `thm:diverging`, `thm:global`, `thm:global-main`, `thm:multiG`, `thm:proportional`, `thm:rate`, `thm:tilted`
