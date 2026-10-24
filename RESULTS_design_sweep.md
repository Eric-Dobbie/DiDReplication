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

> **SUPERSEDED (commit 55ddc68).** This section has three defects, corrected in the
> follow-up revision: (1) it compared a per-period *slope* MDE against a *level*
> threshold (0.303); (2) its control pools were hard-coded at 41/741 rather than
> derived from each configuration's geometry; (3) the claim "√(1/41) forces MDE/σ ≥
> 0.62 for every k" and the "resolvable share" headline follow from those two
> errors. See the revised section below for the corrected analysis.

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

---

## Weight vs. checkability (MDE), swept over geometry — REVISED (supersedes the section above)

Corrects the three defects flagged against 55ddc68. Script `design_sweep_checkability.R`,
2,016 configs (216 geometry configs × never-treated count N_c ∈ {10,41,120}, plus an
1,800-config random-placement arm at N_c=41). Deterministic; no outcomes.

**Fix 1 — units (level-on-level).** The treatment effect is a level, so the MDE is now
a level. We report (a) the **jump MDE** — the smallest detectable *uniform level*
pre-trend, `MDE_jump/σ = √(λ*(k)·(k+1)/k)·√(1/n_g+1/n_c)` — as the headline level
measure, and (b) the slope MDE converted to an implied mean post bias
(`slope·h`, `h=(T−gp)/2`). The earlier "trend MDE" was a per-period *slope*; comparing it
to the level threshold is what produced the spurious "2009 is resolvable". **Under the
correct level comparison nothing changes sign except that conclusion:** stacked 2009's
detectable pre-trend is **jump 1.56σ / trend-bias 1.04σ** (geometry pool) — ≈ 4.9× / 3.3×
the |β|=0.10 effect — not 0.18σ. `weight_checkability.R` was updated the same way and its
CSV/table/figure regenerated; its stacked-2009 row now reads jump MDE 0.482 (1.46σ,
4.87× |eff|), trend-bias 0.322 (0.98σ, 3.26× |eff|).

**Is anything resolvable after fix 1? No.** Resolvable (MDE_jump/σ ≤ 0.303) requires
`√(1/n_g+1/n_c) ≤ 0.303/C(k)`; even at k=1, `C(1)=3.96`, so it needs n_g ≈ 170. No swept
cohort (max n_g=8) comes close. Resolvable share is **0 in every scheme, every geometry,
and for P&P** — so it is dropped as a headline (fix 7).

**Fix 7 — headline is the weight-weighted MDE distribution.** Across the sweep the
weighted-median jump MDE sits at **1.5–3σ** (≈ 5–9× the effect) for every scheme. P&P is
**typical, not exotic: 42nd–54th percentile** (wtd-median 1.84–1.88σ vs grid-median
1.86–1.88σ). See `design_checkability_wmed.*`; `design_checkability_scatter.*` shows every
cohort — grid and P&P — right of the resolvable line.

**Fix 2 — no pool-floor claim; binding reasons tabulated instead.** For each cohort the
*binding* reason is one of `no_pre` (k<1), `size` (k/n_g>1), `underpowered` (checkable but
jump MDE>0.303), `resolvable`. Weight share by reason (`design_checkability_binding.csv`):

| scheme | P&P no_pre | P&P size | P&P underpowered | P&P resolvable | grid-median size |
|---|---|---|---|---|---|
| simple   | 0.06 | 0.65 | 0.29 | 0 | 0.40 |
| group    | 0.03 | 0.80 | 0.17 | 0 | 0.56 |
| dynamic  | 0.12 | 0.53 | 0.35 | 0 | 0.31 |
| calendar | 0.12 | 0.53 | 0.35 | 0 | 0.31 |
| stacked  | 0.01 | 0.90 | 0.09 | 0 | 0.68 |

The binding constraint is overwhelmingly the **size** gate (k/n_g>1: too many pre-leads
relative to treated units), not a control-pool floor; the remaining checkable weight is
**underpowered**, never resolvable.

**Fix 4 — pools from geometry, not 41/741.** CS family + SA use the never-treated count
N_c; stacked/CS-NYT use N_c + not-yet-treated at g. P&P's geometry-derived stacked pool
runs **70 (2001) down to 41 (2020)** — e.g. 2009 → 54 — far below the committed 741,
because the 741 folds in ~700 *always*-treated (always-no-requirement) units that are not
valid not-yet-treated controls. Using the honest smaller pool makes 2009 slightly *less*
checkable (jump 1.56σ vs 1.46σ at 741), reinforcing the conclusion.

**Fix 3 — gate sensitivity (k/n_g ≤ 1 vs < 1).** The boundary case is k/n_g = 1 exactly.
In the grid, 80 of 2,430 cohort-cells (3.3%) sit on it. For P&P the single flip is
**stacked 2009** (k=8, n_g=8 → 1.00): *checkable* under ≤1, *uncheckable* under <1; all
other P&P cohorts are off the boundary. **Memphis excluded** (its double 2004/2009 switch
dropped, n_g→7): k/n_g = 8/7 > 1, so 2009 fails the size gate outright under either
convention. So 2009's checkability hinges entirely on one agency counted in two stacks and
on a weak-inequality gate.

**Fix 5 — P&P anchor via the same grid functions** (its timing vector + CS cohort sizes +
N_c=41), so it is comparable and placed in the distribution with a percentile. The
committed real-data point (741 pool) appears only as a separate reference marker, with no
percentile.

**Fix 6 — dynamic vs calendar.** Not a bug: `did::aggte` places **identical weight on each
cohort** under dynamic and calendar (confirmed independently by numerically differentiating
`overall.att` w.r.t. each att(g,t) atom — cohort marginals match to 3e-13, reproducing
`task2_cohort_weights.csv`). They diverge only in *within-cohort, across-period* atom
weights (e.g. 2009's 2009 vs 2020 post-years get 0.012/0.021 under dynamic, reversed under
calendar), so the overall ATTs differ (dynamic 0.0549, calendar 0.0486) while the cohort
weights do not. In a balanced panel the two coincide even at the atom level.

Outputs: `output/design_checkability_grid.csv`, `output/design_checkability_cohorts_grid.csv`,
`output/design_checkability_binding.csv`, `output/tab_design_checkability.{csv,tex}`,
`figures/design_checkability_{wmed,scatter}.*`.
