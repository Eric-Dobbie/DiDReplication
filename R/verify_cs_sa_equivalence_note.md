# Note: Why Callaway–Sant'Anna and Sun–Abraham coincide in our figures

**Summary.** In the fatal-encounters analysis, the Callaway & Sant'Anna (CS) and
Sun & Abraham (SA) estimators return numerically identical results — the same
aggregate ATT (+0.045) and coincident cohort- and component-level points on the
weight-vs-estimate figures (`plot_cohort_level.R`, `plot_weight_vs_beta.R`). This
is **not a coincidence of our data**. It is a known algebraic property of the two
estimators that holds whenever (a) both use **never-treated** units as the control
group and (b) the panel is **balanced**. Under those two conditions the estimators
are built from the *same* cohort-by-period 2×2 DiD comparisons and their
aggregation weights collapse to the *same* values, so both the disaggregated
"atoms" and the aggregate ATT match to machine precision. The claim is verified
below with a reproducible check (`R/verify_cs_sa_equivalence.R`), and confirmed
against the package documentation and source.

---

## 1. The two estimators build the same 2×2 comparisons

**Same building block.** For a treated cohort `g` and calendar period `t ≥ g`,
CS's group-time effect and SA's cohort-event effect are the identical
never-treated 2×2 difference-in-differences,

```
ATT(g, t)  ==  CATT(g, e = t − g)
           =  [ Ȳ_g(t) − Ȳ_g(g−1) ]  −  [ Ȳ_N(t) − Ȳ_N(g−1) ],
```

where `Ȳ_g` averages the treated cohort, `Ȳ_N` averages the never-treated
controls, and `g−1` is the common base period.

- **CS.** `did::att_gt(..., control_group = "nevertreated")`. The docs:
  *"The default is `nevertreated` which sets the control group to be the group of
  units that never participate in the treatment. This group does not change across
  groups or time periods"* (`?att_gt`). For post-treatment cells the base period is
  `g−1`; the docs note the base-period choice *"results in the same post-treatment
  estimates of ATT(g,t)'s"* — so the post-period cells are unconditional on the
  `varying`/`universal` setting (confirmed empirically below).
- **SA.** `fixest::feols(y ~ sunab(cohort, period) | id + period)`. The
  `sunab` docs: *"By default the never treated cohorts are taken as reference and
  the always treated are excluded from the estimation"* and the reference relative
  period is `ref.p = -1`, i.e. `g−1` (`?sunab`). So SA's interacted regression
  identifies each `CATT(g,e)` off the *same* never-treated reference and the *same*
  base period `g−1` as CS.

Same treated cohort, same control pool (never-treated), same base period ⇒ the CS
and SA per-cell estimates are algebraically the same number.

## 2. On a balanced panel the aggregation weights also coincide

The aggregate ATT is a weighted average of the post-treatment atoms; the two
packages define those weights differently, but a balanced panel makes them equal.

- **CS "simple" aggregation** (`did:::compute.aggte`, `type == "simple"`):
  ```
  simple.att <- sum(att[keepers] * pg[keepers]) / sum(pg[keepers])
  ```
  where `keepers` are the post cells (`group ≤ t`) and `pg = P(G = g)` is the
  *unit share* of cohort `g`. So each post cell `(g,t)` is weighted ∝ `N_g`
  (the number of units in cohort `g`).
- **SA `agg = "att"` aggregation** (`fixest:::aggregate.fixest`):
  ```
  shares <- colSums(sign(mm[, e>=0 cols]));  shares <- shares/sum(shares)
  ATT    <- sum(shares * coef)
  ```
  i.e. each post cell `(g,e)` is weighted by the *number of observations* loading
  on that cohort×event column. In a balanced panel every event `e` of cohort `g`
  is populated by exactly `N_g` units, so `shares(g,e) = N_g`.

Both schemes therefore place weight ∝ `N_g` on the *same* set of post cells, so
after normalization the weights are identical — hence the aggregate ATTs match,
not merely approximately. (This is the mechanism the atom-extraction comments in
`extract_atoms.R` describe for each estimator; here we show the two reduce to the
same thing.)

## 3. Reproducible check (`R/verify_cs_sa_equivalence.R`)

The script simulates a **balanced** staggered panel (cohorts 2004/2007/2010 plus
never-treated controls, 200 units × 15 years, heterogeneous dynamic effects so the
estimators are not degenerately equal), runs CS and SA, and joins the atoms on
`(cohort g, event time e)`.

**Baseline — balanced + never-treated:**

| quantity | result |
|---|---|
| cells compared (all `e`) | 42 |
| max \|ATT(g,t) − CATT(g,e)\| over all cells | **5.4 × 10⁻¹³** |
| max over post cells (`e ≥ 0`) | 5.4 × 10⁻¹³ |
| aggregate ATT, CS vs SA | 0.2511660546 vs 0.2511660546 (\|diff\| = **1.7 × 10⁻¹³**) |

The match is at machine precision, **cell by cell**, not just in the aggregate.

**Base-period nuance.** Post-treatment cells (`e ≥ 0`, the ones that drive the ATT)
match SA under *both* `base_period = "varying"` (the paper's default) and
`"universal"`. The *pre-period placebos* (`e < 0`) only line up with SA's
event-study coefficients under `base_period = "universal"`, which fixes the
reference at `g−1` exactly as `sunab` does; under `"varying"` they use a rolling
`t−1` base and differ by construction. This matters only for the pre-trend panel of
an event-study plot, not for the ATT.

## 4. Breaking each assumption (divergence appears exactly where theory says)

| scenario | CS ATT | SA ATT | \|agg diff\| | cell max diff | expected |
|---|---|---|---|---|---|
| balanced + never-treated | 0.25117 | 0.25117 | 1.7 × 10⁻¹³ | 5.4 × 10⁻¹³ | **equal** |
| **unbalanced** + never-treated | 0.11681 | 0.25951 | 1.4 × 10⁻¹ | — | diverge |
| balanced + **not-yet-treated** (CS) | 0.24985 | 0.25117 | 1.3 × 10⁻³ | 1.3 × 10⁻¹ | diverge |

- **Unbalance the panel** → the aggregate moves by ~0.14. Two channels break at
  once: `did::att_gt` defaults to `allow_unbalanced_panel = FALSE`, so it
  *rebalances by dropping* any unit not seen in every period (it prints
  *"…converting to balanced panel by dropping them"*), while `sunab` keeps all
  rows — so they no longer use the same observations; and the cohort×event obs
  counts become uneven, so SA's obs-count weights no longer equal CS's fixed
  group-share weights. Either channel alone breaks the equivalence.
- **Switch CS to not-yet-treated** → individual post cells diverge by up to 0.13,
  because CS now compares each cohort against a *different, larger* control pool
  (later-treated units + never-treated) while `sunab` still identifies off the
  never-treated reference. Note the *aggregate* here happens to stay close
  (\|diff\| = 0.0013) even though the *cells* differ substantially — which is
  exactly why the check is done at the cell level, not just on the pooled ATT.

The script ends with `stopifnot(...)` assertions on all four outcomes, so it fails
loudly if a package update ever breaks the equivalence.

## 5. External confirmation

The equivalence is also stated in the methodological literature: the
de Chaisemartin & D'Haultfœuille survey notes that the Sun & Abraham estimator is
identical to Callaway & Sant'Anna when both use the same (never-treated) control
group, and Sun & Abraham (2021) discuss the relationship of their interaction-
weighted estimator to Callaway & Sant'Anna directly.

---

### Suggested footnote for the paper

> Under a balanced panel with never-treated controls, the Callaway–Sant'Anna and
> Sun–Abraham estimators are numerically identical: both are built from the same
> cohort-by-period never-treated 2×2 comparisons (same treated cohort, same control
> pool, same `g−1` base period), and on a balanced panel their aggregation weights
> both reduce to the cohort's unit share, so the disaggregated components and the
> aggregate ATT coincide to machine precision. The apparent single "CS/SA" point in
> Figures X–Y is therefore two overlapping estimators, not one. The equivalence
> breaks — as expected — if the panel is unbalanced or if a different control group
> (not-yet-treated) is used. See `R/verify_cs_sa_equivalence.R`.

### Sources consulted

- `?att_gt` (package `did` 2.5.1): `control_group` and `base_period` argument
  descriptions.
- `did:::compute.aggte`, `type == "simple"` branch: `simple.att <- sum(att[keepers]
  * pg[keepers]) / sum(pg[keepers])`.
- `?sunab` (package `fixest` 0.14.2): `ref.c` (never-treated reference) and
  `ref.p = -1` defaults.
- `fixest:::aggregate.fixest`, `agg = "att"` branch: `shares <-
  colSums(sign(mm)); shares <- shares/sum(shares); sum(shares * coef)`.
- Reproducible check and machine-precision results: `R/verify_cs_sa_equivalence.R`,
  outputs `output/cs_sa_cell_equivalence.csv` and
  `output/cs_sa_equivalence_summary.csv`.
