# Verification runs before Task 3 (Tasks A1–A4)

Context: Tasks 1, 2, 4.1 accepted. Task 4.2 **rejected** pending verification of the
pre-trend × weight results (a joint Wald p = 1.000 looked like an artifact). Run in
order; each can change the next. Environment: R 4.3.3, `did` 2.1.2, `fixest` 0.13.2.

---

## A1 — Wald test diagnostics  → the p-values were an artifact (confirmed), but not the suspected mechanism
`R/taskA1_wald_diagnostics.R` → `output/taskA1_wald_diagnostics.csv`.

- **Not rank deficiency / not ginv.** Every CS (varying) and SA pre-block is full
  rank (0/12 rank-deficient); `ginv` tolerance 1.49e-8 dropped no singular values
  (CS-2009 smallest sv 0.026 ≫ tol·max = 2e-7). So the original explanation (ginv on
  a rank-deficient vcov) is wrong.
- **The real bug: `V_analytical` scaling.** did returns `V_analytical` as a sparse
  matrix whose SE is `sqrt(diag(V)/n)` (n = 71 units; ratio to did's own analytical
  SE = 1.0000). The Task-4.2 code used `diag(V)` **without `/n`**, inflating the vcov
  71×, shrinking every CS Wald statistic 71×, and sending p→1. With the `/n` fix,
  CS(varying) pre-trends **reject for most cohorts** (2009 Wald 31.4 p=0.0001 vs SA
  25.1 p=0.0015; 2008/2012–2020 p≈0). The 57× CS-vs-SA Wald gap **disappears**.
- **Rank deficiency does appear — in CS *universal*.** The universal base period
  includes the mechanically-zero reference cell (one exactly-zero singular value), so
  those blocks are rank-deficient and the p is reported as **NA** by the decision
  rule (not a pseudo-inverse p). Excluding the base cell recovers SA's block.
- Since agency = unit here, the correctly-scaled analytical vcov is already
  agency-clustered.

## A2 — Are CS and SA pre-estimates the same object?  → yes under universal, no as-run
`R/taskA2_base_period.R` → `output/taskA2_base_period.csv`.

- **CS as run in Task 4.2 used `base_period = "varying"`.** CS-varying ≠ SA:
  max |difference| = **1.146**. Different pre-treatment estimand.
- **CS-universal = SA exactly**: max |difference| = **3.9e-13**. Same estimand.
- SE ratio (CS analytical `/n` ÷ SA cluster-robust) = **0.895** (was a spurious 7.54
  before the `/n` fix). So under the matched convention the SEs are comparable too.
- **Verdict:** the Task-4.2 "identical estimates, opposite verdicts, inference-only"
  framing was wrong on two counts — the figure compared different estimands
  (varying vs SA) **and** the opposite verdicts were the vcov scaling bug. Under the
  matched (universal) convention CS and SA agree in estimate, SE, and pre-trend
  verdict. **The inference-only story is dropped.**

## A3 — Per-cohort atom comparison, CS vs stacked (gates Task 3)  → ambiguous, leans "noise concentrated by weighting"
`R/taskA3_atom_bootstrap.R` → `output/taskA3_atom_bootstrap.csv`,
`output/taskA3_atom_channel_draws.csv`. Paired agency-cluster bootstrap, B=500; fast
estimators verified to reproduce the committed atoms exactly; weights held at
full-sample values.

- **6 of 11 cohorts** show a CS−stacked atom difference distinguishable from 0
  (p<0.05): 2004, 2012, 2013, 2017, 2018, 2019 — all **low/moderate weight**.
- **The two highest-weight cohorts are NOT distinguishable:** 2002 (weight 0.328,
  diff +0.215, **p=0.061**) and 2009 (weight 0.290, diff −0.054, **p=0.664**). They
  carry 62% of the CS aggregation weight.
- **Atom channel:** point −0.117; bootstrap 95% CI **[−0.184, +0.018] — crosses
  zero** (complete-draws-only CI [−0.170, +0.001] also reaches zero). The channel is
  dominated by 2002 (≈−0.071 of −0.117), whose difference is only marginal.
- **Reading (reported, not resolved):** genuine identification differences exist in
  several low-weight cohorts, so it is **not pure noise**; but the headline −0.117
  atom channel is **not distinguishable from zero** and rests on the highest-weight
  cohort (2002) whose difference is not significant. This leans toward the brief's
  second framing — *cross-estimator disagreement here is noise concentrated by the
  weighting scheme* (on 2002/2009) — rather than "the estimators cleanly identify
  different objects." **Task 3 should be designed to demonstrate the second thing.**

## A4 — Pre-trend figure rebuilt on magnitude  → replaces the retired p-value figure
`R/taskA4_pretrend_magnitude.R` → `output/taskA4_pretrend_magnitude.csv`,
`figures/pretrend_magnitude.{pdf,png}`. y = largest |pre-treatment coefficient| with
95% CI; x = aggregation weight; size = pre-periods available; label = cohort (n =
treated agencies); facet per estimator; CS uses universal base with the corrected
`/n` SE.

- Confirms A2 visually: CS(universal) and SA panels coincide in point estimates.
- The claim it supports: **heavy weight ≠ checkability.** Cohort 2002 (CS weight
  0.328) has a large pre-coefficient (0.40) with a CI crossing zero on a **single**
  pre-period; 2009 (weight 0.290) has a wide CI crossing zero. The cohorts with tight
  CIs clearly off zero (2008, 2012, 2017, 2019 — real pre-trend violations) carry
  almost no weight. A tight interval around zero and a wide interval around zero are
  now visually distinct, which the p-value version hid.
- The p-value figure (`pretrend_vs_weight`, `task4_pretrend_weight.*`) is **retired**
  (A1 showed its p-values were unsound).

---

### Net effect on the earlier write-up
`RESULTS_remaining_runs.md` Task 4.2 is superseded: there is **no** "CS supported /
SA not, inference-only" finding. The corrected statement is that pre-trend support is
**weakest exactly where aggregation weight is highest** (2002, 2009), under all three
estimators, and that is a magnitude/power statement, not a cross-estimator one.
