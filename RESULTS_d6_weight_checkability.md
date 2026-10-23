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
| CS (simple)        | deck     | +0.045 | 0.35 | −0.78 |
| CS (simple)        | droppers | +0.065 | 0.37 | −0.80 |
| stacked (R2 FE)    | deck     | −0.099 | 0.06 | +0.22 |
| stacked (R2 FE)    | droppers | −0.119 | 0.00 | +0.31 |
| stacked Cengiz ±4 (NT) | deck | −0.089 | 0.00 | −0.48 |
| stacked Cengiz ±4 (NT) | droppers | −0.093 | 0.00 | +0.22 |
| imputation (BJS)   | droppers | −0.043 | 0.37 | −0.80 |

(imputation ATT and shares updated after the item-1 singleton fix, which restores cohort 2001.)

- **CS & imputation** put ρ ≈ −0.8 on (weight, pre_g) and **35–37% of all weight on
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
3. **CS vs. imputation simple weights (droppers) — now EXACTLY identical** (max abs diff
   0.0000 across all 11 cohorts) after the item-1 singleton fix. Before the fix they
   differed by up to 0.027 because imputation dropped cohort 2001 (Denver CO's lone
   untreated pre-period, 2000, was removed as a fixed-effect singleton under the default
   `fixef.rm="perfect_fit"`). With `fixef.rm="none"` that one observation identifies
   Denver's agency intercept, cohort 2001 is retained, and the two simple-weight menus
   coincide to machine precision — confirming Layer 1 is the same menu. (The per-cohort
   *estimates* still differ — see the bridge below.)
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
   −0.043) despite *exactly identical* cohort weights. The disagreement is entirely in the
   per-cohort estimates τ̂_g (Layer 2 / counterfactual), not the weights (Layer 1), and the
   bridge (follow-up 3) localizes it to the CS→imputation step rather than the control
   group — the motivation for the optional counterfactual-mass figure (not built here).
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

---

## Follow-ups (items 1–4)

**1. Imputation singleton fix (done).** Stage 1 now runs with `fixef.rm="none"`, so a
lone untreated observation still identifies its agency intercept (α_i = Y − γ_t). Cohort
2001 (Denver) is retained; imputation ATT moves −0.047 → **−0.043**; CS and imputation
simple weights become exactly identical (validation 3 above). Panel C and the table are
updated.

**2. Panel B2 — stacked Cengiz ±4 with never-treated controls (robustness, done).** A
balanced ±4 window (cohorts 2004–2016 with full windows), controls = the 41 never-treated
only (not P&P's 741 never+always), Cengiz interactive FE (agency×stack + year×stack =
the R2 spec), w_s = V_s/V. **ATT −0.089 (deck) / −0.093 (droppers)** — close to the
full-panel stacked −0.099 / −0.119 — and the weighting pattern is the same: all weight on
mid-panel *testable* cohorts (2009 dominant, then 2013), untestable share 0.00. So Panel
B's finding (stacked weight avoids the untestable early cohorts) is **robust** to using
clean never-treated controls and a balanced window; it is a property of the V_s geometry,
not of P&P's control pool. Added as the fourth panel in all three figures.

**3. Bridge CS-NT → CS-NYT → imputation (droppers).** `d6_bridge_cs_imputation.csv`.
Overall ATT: **CS-NT +0.065 → CS-NYT +0.077 → imputation −0.043.** The control group
(never-treated → not-yet-treated) barely moves the estimate (+0.065 → +0.077), and the
per-cohort τ̂_g are almost unchanged between CS-NT and CS-NYT. **The entire CS–imputation
gap is the CS→imputation step** — i.e. the estimand/counterfactual (imputed agency+year FE
Y(0) vs. aggregated att(g,t)), not the comparison group. It localizes to a few small
cohorts whose τ̂ differ sharply: 2013 (CS +0.30 vs imp +0.03), 2017 (+0.74 vs +0.02),
2018 (+0.01 vs −0.61) — these flip the weighted sum's sign. This is the Layer-2 channel
the optional counterfactual figure would visualize.

**4. Pre-tests for the highest-weight testable cohorts, 2009 and 2013.**
`d6_pretests_2009_2013.csv`. Each estimator's own comparison group; pre-periods only.

| cohort | test | CS (NT) | stacked (No-Change) | imputation | Cengiz (NT) |
|---|---|---|---|---|---|
| 2009 | linear slope (p) | −0.006 (0.71) | −0.018 (0.18) | −0.005 (0.72) | −0.006 (0.71) |
| 2009 | joint Wald (p) | 3.70 (0.00) | 4.33 (0.00) | 4.61 (0.00) | 3.70 (0.00) |
| 2009 | MDE (slope) | 0.041 | 0.037 | 0.040 | 0.041 |
| 2013 | linear slope (p) | −0.023 (0.019) | −0.021 (0.014) | −0.021 (0.026) | −0.023 (0.019) |
| 2013 | joint Wald (p) | 3.81 (0.00) | 6.25 (0.00) | 4.17 (0.00) | 3.81 (0.00) |
| 2013 | MDE (slope) | 0.026 | 0.024 | 0.026 | 0.026 |

- **2009**: the linear pre-trend slope is small and **not significant** (p 0.18–0.72) under
  every comparison group — passes the slope test — **but the joint pre-trend Wald rejects**
  (p≈0).
- **2013**: a significant negative pre-trend slope (≈ −0.02/yr, p 0.01–0.03) — **fails the
  slope test** — and the joint Wald also rejects.

**Caveat on the joint Wald.** It rejects (p≈0) for *every* cohort/estimator, with 8–12
restrictions against only 3–7 treated clusters. This is exactly the over-rejection regime
documented in the B-chain (`RESULTS_verification_B.md`): joint pre-trend Wald tests
massively over-reject when restrictions ≳ treated clusters, so these joint p-values are not
trustworthy evidence of a real pre-trend. The **one-restriction linear-slope test is the
reliable read**: 2009 is flat (passes), 2013 trends (fails). Either way, the two
highest-weight *testable* stacked cohorts do **not** both pass a clean pre-test — the
weight that isn't on untestable cohorts sits on cohorts that are themselves borderline
(2013 fails) or fail the joint test (2009), so "testable" is not "passes".

New/updated outputs: `d6_bridge_cs_imputation.csv`, `d6_pretests_2009_2013.csv`; the table
and all three figures now carry the fourth (Cengiz) panel.
