# Deterministic design sweep over panel geometry (replaces Task 3)

**Not a simulation.** No outcome variable, no DGP, no sampling of data, no
replications. Aggregation weights are exact deterministic functions of the design
(panel length `T`, adoption-period vector `g`, cohort sizes `n_g`, control pool
`N_c`). The core is `design_weights()` in `R/design_sweep_core.R`; every figure is
that function enumerated over configurations. (The "random-placement arm" samples
*designs*, not data.)

**Core validated** (`design_sweep_core.R` self-test): CS `simple`/`group`/`dynamic`/
`calendar` cohort weights reproduce the committed P&P `aggte` weights to 1e-7; the
exact stacked `V_s` (within-stack two-way demean) matches an `felm` residualization
to 1e-15.

---

## Figure C — weight on uncheckable cohorts (`design_sweep_figureC.R`)
`figures/design_uncheckable_share.*`, `output/design_sweep_grid.csv` (1,872 configs),
`output/design_uncheckable_named.csv`. Uncheckable = `k/n_g > 1` (B1/B2);
`k = pre_g − 1`. Share of total aggregation weight on uncheckable cohorts, per scheme.

**P&P share: simple 0.65, group 0.80, dynamic 0.53, calendar 0.53, stacked 0.91.**
So 53–91% of P&P's aggregation weight lands on cohorts whose parallel-trends
assumption cannot be checked — under stacked, 91%. P&P's uncheckable cohorts are
2003, 2004, 2008, 2009, 2012–2020 (all but 2001, 2002).

**P&P's position: the 74th–80th percentile** of the swept distribution (per scheme).
Reported, not resolved: **the problem is common, not exotic** — most designs put a
large share of weight on uncheckable cohorts (the distributions are centered well
above zero, with a spike at 1.0 for late-placement designs), and P&P sits high in
that distribution but short of the tail. What pushes P&P up is its many late-adopting
small cohorts (2008–2020, mostly 1–3 treated units), which have large `k` and small
`n_g`. Driver, one parameter at a time (in the grid): **late placement** drives the
share toward 1.0; **skewed sizes** (single-unit cohorts) raise it; larger `G` raises
it; larger `T` raises it (more pre-periods per late cohort).

## Figure A — weight vs pre-period count (`design_sweep_figureAB.R`)
`figures/design_weight_vs_pre.*`. Fixed T=30, evenly spaced cohorts. Expected shapes
**confirmed**: CS `simple` strictly decreasing and linear in `pre_g` (weight ∝
`n·(T−pre)`); CS `group` flat (equal sizes); stacked peaks mid-panel (`q(1−q)` maxed
at `q=0.5`). **Departure to report: `dynamic` and `calendar` coincide exactly** (one
line hidden under the other) — true here and in the P&P weights. The unequal-sizes
panel separates size from timing: `group` then tracks `n_g`, and the stacked peak
shifts toward the larger cohorts.

## Figure B — agreement between weight vectors (`design_sweep_figureAB.R`)
`figures/design_weight_agreement.*`, `output/design_sweep_corL1.csv`. Ex-ante
diagnostic: `cor(w_CS_simple, w_stacked)` and `L1 = Σ|w_CS_simple − w_stacked|` vs
cohort placement (centroid). cor → **+1 for late placement** (both schemes pile onto
the late small cohorts); **negative for 171 of 688 configs**, concentrated in random
and early-clustered placements at mid-centroid. **Bound:** divergence attributable to
weighting alone is `L1·(b−a)/2` for cohort effects in `[a,b]`; e.g. `[−0.5, 0.5]`
gives up to **0.67** across the sweep. **P&P:** cor 0.65, L1 0.65, centroid 0.56 →
bound 0.32 for `[−0.5, 0.5]` (positively correlated weights, but a large absolute
gap).

## Validation — the closed form (`design_sweep_validation.R`)
`output/design_validation_pp.csv`. Derivation: for a balanced stack with an all-zero
common control pool, `resid_{it} = (a_i − p̄)(b_t − q̄)`, so
`V_s = N·p(1−p)·T·q(1−q)` **exactly**.
- **(0) Confirmed:** max relative `|exact − closed|/exact` = **5.4e-16** over 200
  balanced common-pool designs → the derivation is correct; the formula may be
  presented **with the balance + common-pool condition attached**.
- **(1) Imbalance** (missing cells): 5% → 4% median deviation, 20% → 28%, 40% → 74%.
- **(2) Balancing weights** (controls downweighted): weight 0.5 → 10%, 0.1 → 91%.
- **P&P (11 stacks):** `w_exact` vs `w_closed-form` correlation 0.986, L1 0.123 —
  good except (a) **2002** has no controls (`p=1`, the formula gives 0, but
  `V_s=5.79` from the within-treated adopter contrast — outside the formula's scope),
  and (b) **2009** over-predicted ~14%. **2009-vs-2013 `V_s` ratio: closed-form 2.75
  vs exact 2.41**, reproducing the earlier 2.8-vs-2.4 check. The looseness is the
  real stack's departures from the clean design — the mixed control pool (700
  always-`no.req=1` units perturb the effective shares) and, for 2002, the absence of
  controls — not an error in the derivation.

## Identification bookkeeping
`design_weights()` carries each cohort's status per estimator (`ok`,
`single_treated`, `collinear:no_pre`, `collinear:FE`, `no_contrast`) and sets the
stacked weight to NA with the reason where `V_s = 0`; the grid records the count of
unidentified cohorts per config. No outcome variable appears anywhere in this task.

---

## Weight vs. checkability (MDE), swept over geometry (`design_sweep_checkability.R`)

The same weight-vs-checkability comparison as `weight_checkability.R` / Figure C,
now run across the design space. Deterministic: no outcomes, no DGP. Per cohort,
per scheme, per config we compute the aggregation weight and the analytic
**minimum detectable pre-trend (MDE)** at 80% power / α=0.05, expressed in units
of σ (the outcome noise SD): `MDE/σ = sqrt(λ*(k)/[a'(I+11')⁻¹a]) · sqrt(1/n_g + 1/n_c)`,
a pure function of the geometry (`k`, `n_g`, `n_c`) — σ cancels, so no outcome is
needed. Gates: `no_pre` (k<1), `size` (k/n_g>1, the B1/B2 rank gate). "Resolvable"
= checkable **and** `MDE/σ ≤ 0.303` (= a P&P-scale |β|=0.10 pre-trend at σ=0.3297).
**Control pools are scheme-specific and realistic**, as in the data: CS family + SA
use the never-treated pool (41); stacked uses the not-yet-treated pool (741). The
sweep varies geometry (T∈{15,21,30}, placement, sizes) — 1,872 configs, identical
universe to Figure C. Output `MDE/σ` validated against the committed
`weight_checkability_cohorts.csv` to 4e-5 (uniroot tol).

**The control pool, not the geometry, gates CS/SA.** Under the never-treated pool
(n_c=41) the term `sqrt(1/n_c)` alone is 0.156, so `MDE/σ ≥ C(k)·0.156 ≥ 0.62` for
every k — above the 0.303 resolvable line. **No panel geometry in the sweep makes a
CS/SA pre-trend resolvable at a P&P-scale effect**: grid median resolvable share = 0
for simple/group/dynamic/calendar, and P&P's 0.00 is the generic case, not a quirk.
The binding constraint is the small never-treated pool, and it cannot be fixed by
re-timing or re-sizing cohorts.

**Only the stacked estimator's large pool can resolve anything, and P&P sits at the
favorable extreme.** With n_c=741 the pool floor drops to `sqrt(1/741)=0.037`, so a
large mid-panel cohort can clear 0.303. Across the sweep the median stacked
resolvable share is still 0 (most geometries load high-k cohorts that fail the size
gate), but **P&P's stacked resolvable share is 0.369 — the ~98th percentile**, and
its weight-weighted median `MDE/σ` is 0.178 — the **0th percentile (the single most
checkable design in the entire sweep)**. Both are driven by one cohort, 2009
(n_g=8, k=8, pool 741). Even so, **63% of stacked weight still lands on uncheckable
cohorts.**

| scheme | grid median resolvable | grid median uncheckable | grid median wtd-med MDE/σ | **P&P** resolvable | **P&P** uncheckable | **P&P** wtd-med MDE/σ |
|---|---|---|---|---|---|---|
| simple   | 0.000 | 0.605 | 1.313 | 0.000 | 0.710 | 1.877 |
| group    | 0.000 | 0.722 | 1.040 | 0.000 | 0.833 | 1.877 |
| dynamic  | 0.000 | 0.552 | 1.313 | 0.000 | 0.647 | 1.877 |
| calendar | 0.000 | 0.552 | 1.313 | 0.000 | 0.647 | 1.877 |
| stacked  | 0.000 | 0.733 | 0.854 | **0.369** | 0.631 | **0.178** |

P&P percentiles (per scheme): resolvable share — CS/SA ~98th (tied at the 0 mass
point, i.e. generic); stacked 98th (exceptional). wtd-med MDE/σ — CS/SA 70–82nd
(P&P's checkable cohorts are *harder* to check than the median design, since its one
checkable CS cohort is the small 2002); stacked 0th (most checkable in the sweep).

**Figures.** `design_checkability_scatter.*` — the comparison itself, pooled across
all grid geometries (weight vs MDE/σ, by scheme), P&P cohorts overlaid: every CS/SA
geometry sits right of the 0.303 line, P&P's stacked 2009 is the lone point left of
it. `design_checkability_share.*` — distribution of the resolvable-weight share
(spike at 0 everywhere; P&P stacked alone in the right tail at 0.37).
`design_checkability_wmed.*` — distribution of weight-weighted median MDE/σ (P&P CS
in the right/hard tail, P&P stacked at the extreme left/easy edge).

**Abstractions.** The balanced grid gives every cell a count and one `n_g` per
cohort shared across schemes, so it cannot represent two real-data wrinkles carried
exactly in the P&P anchor (read from `weight_checkability_cohorts.csv`): Memphis
making stacked's 2009 n_g=8 vs CS's 7, and the 2002 stack having no clean controls.
Outputs: `output/design_checkability_grid.csv`, `output/design_checkability_cohorts_grid.csv`,
`output/tab_design_checkability.{csv,tex}`.
