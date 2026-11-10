# Remaining-runs results (Tasks 1, 2, 4) — status

**Environment note.** After a container reset wiped the raw data and R, the data
was re-uploaded and the R toolchain rebuilt (conda-forge: R 4.3.3, `did` 2.1.2,
`fixest` 0.13.2, `lfe` 3.1.1). All tasks below are now run under R with real
re-fits; Task 1's pure-recombination numbers reproduced exactly on re-run.

Common cohort set (the 11 the stacked spec identifies, `w_s > 0`): **2002, 2004,
2008, 2009, 2012, 2013, 2014, 2017, 2018, 2019, 2020**.

---

## Task 1 — Counterfactual reweighting  ✅
`R/task1_counterfactual_reweight.R` → `output/task1_reweight_2x2.csv`,
`output/task1_cohort_inputs.csv`.

| Atoms | Weights | Estimate |
|---|---|---|
| CS | CS simple | **+0.056** |
| CS | stacked FWL | +0.075 |
| stacked | CS simple | −0.061 |
| stacked | stacked FWL | **−0.099** |

Diagonal reproduces the known aggregates on the common set (stacked −0.099 exactly;
CS +0.056 restricted vs +0.045 unrestricted — reported side by side). Gap (−0.155) =
**atom channel −0.117** + weight channel +0.019 + interaction −0.057.

**Which channel dominates.** The CS-vs-stacked disagreement is overwhelmingly an
**atom-level** disagreement, not a weighting artifact: swapping only the per-cohort
effects (weights held at CS) moves the aggregate −0.117 (¾ of the gap); swapping only
the weights moves it +0.019. CS and stacked estimate different objects per cohort
(2002 +0.222 vs +0.007; 2017 +0.722 vs −0.059; 2018 +0.230 vs −0.292). The gap cannot
be reweighted away — it lives in the atoms.

---

## Task 2 — Aggregation-scheme sensitivity  ✅
`R/cs_aggregation_sensitivity.R` → `output/task2_aggregation_sensitivity.csv`,
`output/task2_cohort_weights.csv`. Cell weights reconstructed and verified against
`aggte`'s own `overall.att` (|diff| ≤ 1e-17 for all four schemes).

| scheme | overall | SE | w(2002) | w(2009) |
|---|---|---|---|---|
| simple | +0.045 | 0.065 | 0.290 | 0.256 |
| group | +0.106 | 0.054 | 0.167 | 0.233 |
| dynamic | +0.055 | 0.067 | 0.353 | 0.196 |
| calendar | +0.049 | 0.055 | 0.353 | 0.196 |

All four CS aggregates are **positive and insignificant** (SE 0.054–0.067); the CS
sign is robust to the aggregation scheme. `simple`/`dynamic`/`calendar` put the most
weight on 2002 (0.29–0.35); `group` spreads weight onto the noisy large-positive late
cohorts and rises to +0.106. The sign fragility is cohort-specific (2002), not
scheme-specific (Task 4.1).

---

## Task 4.1 — Leave-one-cohort-out  ✅
CS: `R/cs_loo.R` → `output/cs_loo_se.csv`; stacked: `output/stacked_loo_both.csv`.

- **CS: dropping 2002 flips the sign** (+0.045 → **−0.028**, SE 0.081). No other
  cohort flips it; dropping 2009 raises it most (+0.119). CS's positive aggregate
  rests entirely on 2002.
- **Stacked: no cohort flips the sign** (plain corrected-FE stays in [−0.115, −0.082];
  most influential drop 2018 → −0.082). Robust where CS is not.
- Several cohorts have a single treated unit (flagged in `n_treated`): cluster-robust
  SEs are unreliable there.

## Task 4.2 — Pre-trend × weight  ⚠️ SUPERSEDED (see `RESULTS_verification_A.md`)
The original p-value figure (`pretrend_vs_weight`, `task4_pretrend_weight.*`) is
**retired**. Verification (Tasks A1–A4) found its CS p-values were a `V_analytical`
scaling artifact (used without `/n`, inflating the vcov 71×), and that the CS run
used `base_period="varying"` — a different pre-treatment estimand from SA. The
corrected analysis (`taskA4_pretrend_magnitude.R` → `figures/pretrend_magnitude.*`)
plots pre-trend coefficient **magnitude with CI**, not p-values, and shows that
pre-trend support is weakest exactly where aggregation weight is highest (2002,
2009) — a power/magnitude statement, not the "CS supported, SA not" cross-estimator
claim, which does not hold.

---

## Task 3 — held for review
Per the brief, the boundary-asymmetry simulation waits until Task 1 is reviewed. Task
1's headline (atom channel dominates the CS-vs-stacked gap) is the thing to check
before the simulation design is fixed. The P&P adoption-timing gradient the brief asks
for (slope of `beta_s` on `g`) will be computed as part of Task 3.
