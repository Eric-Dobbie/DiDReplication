# D6 — weight vs. checkability (one panel per estimator)

`R/d6_weight_checkability.R` → `weight_checkability_table.csv`,
`figures/fig_weight_checkability_d6_{v1,v2,preg}.{pdf,png}` (+ `output/*.png`).

**Claim under test.** Under every modern estimator, the cohorts carrying the most
aggregation weight are the cohorts whose parallel-trends assumption is least
checkable. x = checkability (`mde_g`, log; companion with `pre_g`); y = cohort
aggregation weight; one panel per estimator (CS simple / stacked R2 / imputation BJS).

## Headline result

**Supported for CS and imputation; *not* supported for stacked.**

| estimator | sample | ATT | untestable weight share | ρ(weight, pre_g) |
|---|---|---|---|---|
| CS (simple)      | deck     | +0.045 | 0.35 | −0.78 |
| CS (simple)      | droppers | +0.065 | 0.37 | −0.80 |
| stacked (R2 FE)  | deck     | −0.099 | 0.06 | +0.22 |
| stacked (R2 FE)  | droppers | −0.119 | 0.00 | +0.31 |
| imputation (BJS) | droppers | −0.047 | 0.32 | −0.85 |

- **CS & imputation** put ρ ≈ −0.8 on (weight, pre_g) and **32–37% of all weight on
  untestable cohorts**. The single largest-weight cohort, **2002** (pre_g = 2, < 3
  pre-periods), is untestable in both; **2009** is the other big weight and sits at a
  high `mde_g`. So weight piles onto the least-checkable cohorts — the claim holds.
- **Stacked contradicts the pattern.** Its weight peaks on **2009 and 2013** — mid-panel,
  *testable* cohorts — and its untestable-weight share is ≈ 0 with ρ = +0.3. The FWL
  weight `V_s = p(1−p)·q(1−q)·N` is maximized at mid-sample adoption (q ≈ ½), which is
  exactly where pre-periods are plentiful. The stacked estimator, by its own weighting
  geometry, *avoids* the early untestable cohorts rather than loading them. This is the
  one place the expected pattern breaks, and it is a property of the stacked weights, not
  of the data. (It is also the headline from the design sweep: stacked weight concentrates
  mid-panel.)

Two figure versions: **(i)** CS/stacked on the deck sample + imputation on droppers;
**(ii)** all three on droppers. The pattern is identical across versions.

## Validation (all pass)

1. **Weights sum to 1** within every estimator × sample (to 1e-8).
2. **∑ w_g τ̂_g reproduces each aggregate ATT:** CS deck +0.044765 (target +0.045, diff
   −2.3e-4 = rounding), stacked deck −0.098915 (target −0.099, diff +8.5e-5 = rounding),
   imputation droppers −0.046504 (exact, this is the BJS ATT itself). CS droppers +0.065072,
   stacked droppers −0.119049 reproduce their own re-estimated ATTs exactly.
3. **CS vs. imputation simple weights (droppers) — close but not identical** (max abs diff
   0.027, on the two big cohorts 2002 and 2009). Attribution, as expected: (a) **cohort
   2001 is dropped by imputation but kept by CS** — its sole member Denver CO has only one
   untreated pre-period (2000), a fixed-effect singleton that Stage 1 removes, leaving no
   agency intercept to impute Y(0); CS estimates its ATT directly against never-treated
   controls under a universal base period and needs no within-unit pre-history; (b)
   **panel imbalance** (missing `any.fatal`, observed 2000–2020 only) makes imputation's
   treated-observation shares deviate slightly from CS's group-size simple weights. In a
   balanced panel with no dropped cohorts the two menus coincide.
4. **Cohort 2000 is dropped by CS and imputation** (confirmed absent). Reason: its members
   adopt in the first outcome year (2000), so they have zero untreated pre-period
   observations — `att_gt` drops already-treated-in-first-period units, and imputation has
   no untreated obs to fit their FE.

## Checkability construction

- `pre_g = g − 2000` (observed untreated pre-periods); `post_g = 2020 − g + 1`.
- `mde_g = (z_0.975 + z_0.80)·SE(β̂_trend) = 2.802·SE`, where β̂_trend is the differential
  (treated − comparison) linear pre-trend slope from `feols(any.fatal ~ year*D)` on
  **pre-periods only**, clustered by agency. Comparison group is the estimator's own:
  CS → never-treated; stacked → that stack's clean controls (the 741 No-Change agencies);
  imputation → never-treated ∪ not-yet-treated at g (Liu block). No post-period outcome
  ever enters.
- `testable_g = FALSE` if `pre_g < 3`, or the estimator has no comparison group in its own
  design (stacked 2000–2003 have zero clean controls), or the cohort is dropped (2000).

## Weight movement, deck → droppers

Dropping the 6 reversers/adopters moves cohort weights materially: **CS** max |Δw| = 0.068
(∑|Δ| = 0.27), largest at 2009 (+0.068), 2003 (−0.055), 2004 (−0.038); **stacked** max |Δw|
= 0.086 (∑|Δ| = 0.35), largest at 2009 (+0.086), 2002 (−0.060), 2012 (−0.053). Removing the
adopters concentrates still more weight on 2009 in both estimators.

## Things that contradicted / qualified the expected pattern

1. **Stacked does not show the pattern** (above) — the central caveat.
2. **CS and imputation ATTs disagree in sign** on the same droppers sample (+0.065 vs
   −0.047) despite near-identical cohort *weights*. The disagreement is entirely in the
   per-cohort estimates τ̂_g (Layer 2 / counterfactual), not the weights (Layer 1) — which
   is the motivation for the optional counterfactual-mass figure (not built here).
3. **Stacked 2002 is untestable yet carries weight 0.060** (deck): it has zero clean
   controls, so its V_s comes from the within-treated adopter contrast, not a
   never-treated comparison — hence untestable but non-zero weight. Footnoted in Panel B;
   do not conflate it with CS cohort 2002 (41 never-treated controls).

## Flags / deviations from the brief

- **Package substitution:** `did2s` / `didimputation` could not be installed (CRAN
  unreachable through the proxy; R 4.3.3 predates the binaries). Panel C is a from-scratch
  BJS/Gardner imputation (agency + year FE on untreated obs → impute Y(0) → average τ over
  treated obs); under R5 (no covariates) this targets the identical estimand, and ∑ w_g τ_g
  reproduces its ATT exactly.
- **File names:** the brief asked for `fig_weight_checkability.{pdf,png}` /
  `fig_weight_checkability_preg.pdf`, but those collide with the committed outputs of the
  separate `weight_checkability.R` deliverable. To avoid clobbering it, D6 figures use the
  `fig_weight_checkability_d6_` prefix (`_v1`, `_v2`, `_preg`).
- **Two companions:** only one `pre_g` companion was requested; it is built on the droppers
  sample (version ii) so all three panels are on a common footing.
- **Optional Layer-2 counterfactual figure:** not built (left as the stated optional step).
