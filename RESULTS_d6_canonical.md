# D6 — canonical numbers of record

One file for every number the paper quotes from D6. Produced by `R/d6_verify.R`;
machine-readable copy in `output/d6_canonical_numbers.csv` (160 rows), with
`output/d6_att_untestable_rho.csv` and `output/d6_pretests_all.csv` as the backing tables.
No new estimators or specifications — committed per-cohort estimates are reused and the
2009/2013 pre-tests reproduce the committed values **exactly** (see validation).

Conventions: R5 (never-treated controls unless stated, no covariates, universal base, no
binning); `deck` = all cohorts, `droppers` = absorbing only (excludes the 6 adopters + 2
reversals); imputation stage 1 uses `fixef.rm="none"`. **Final testability rule:**
`testable_g = FALSE` if `pre_g < 3`, or no comparison group in the estimator's own design, or
the cohort is dropped by the estimator, **or `n_treated_clusters = 1`.**

## Validation (all pass)

1. **Reproduction.** Every round-1–5 number in this file reproduces to ≤ 0.001, except one,
   listed below with its cause. The 2009/2013 pre-test p-values (CR1 / WCR / RI-unmatched,
   slope and joint) reproduce the committed `d6_pretests_2009_2013_robust.csv` to **0.000**.
2. **Weights.** `max|Σw − 1| = 9.0e-12`; `max|Σ w_g τ̂_g − ATT| = 7.1e-13` across all
   estimator × sample panels.
3. **s=1 gap identity.** 2001 gap = −2.7e-11 (= 0), 2002 gap = −0.2500; holds to 1e-8.

**The one number that moved:** imputation **ρ(weight, pre_g) = −0.80** (not the −0.85 of
round 3). Cause: the cohort-2001 fix (`fixef.rm="none"`) restored cohort 2001 to the
imputation menu, changing the weight-vs-pre_g rank correlation. (Round 5 already reported
−0.80; the −0.85 predates the fix.) New this pass: standard errors for every ATT (CS via
`aggte`, stacked/Cengiz via clustered `felm`, imputation via an agency cluster bootstrap,
B=499) — these were not computed in earlier rounds.

## ATTs (± SE)

| estimator | sample | ATT | SE |
|---|---|---|---|
| CS-NT | deck | +0.0448 | 0.0652 |
| CS-NT | droppers | +0.0651 | 0.0823 |
| CS-NYT | deck | +0.0562 | 0.0669 |
| CS-NYT | droppers | +0.0765 | 0.0829 |
| imputation | droppers | −0.0429 | 0.0582 |
| stacked (agency×stack FE) | deck | −0.0989 | 0.0398 |
| stacked (agency×stack FE) | droppers | −0.1190 | 0.0473 |
| Cengiz ±4 (NT) | deck | −0.0888 | 0.0678 |
| Cengiz ±4 (NT) | droppers | −0.0933 | 0.0783 |

## Item 1 — untestable weight shares, and the sample question

> **Updated in round 6:** the checkability rule changed pre_g ≥ 3 → pre_g ≥ 2, so the
> untestable shares and n_testable below are superseded for CS-NT/CS-NYT/imputation (2002
> becomes testable): CS/CS-NYT/imputation untestable drop to 0.213 (deck) / 0.243 (droppers);
> stacked and Cengiz are unchanged. The round-6 values are the current canonical numbers (see
> `RESULTS_d6_round6.md` and `d6_round6_changelog.csv`). ATTs, ρ and gap numbers are unchanged.

| estimator | sample | untestable share | n_testable | testable cohorts |
|---|---|---|---|---|
| CS-NT | deck | 0.503 | 5 | 2004, 2009, 2013, 2017, 2018 |
| CS-NT | droppers | **0.537** | 3 | 2009, 2013, 2017 |
| CS-NYT | deck | 0.503 | 5 | 2004, 2009, 2013, 2017, 2018 |
| CS-NYT | droppers | **0.537** | 3 | 2009, 2013, 2017 |
| imputation | droppers | **0.537** | 3 | 2009, 2013, 2017 |
| stacked | deck | 0.241 | 5 | 2004, 2009, 2013, 2017, 2018 |
| stacked | droppers | **0.232** | 3 | 2009, 2013, 2017 |
| Cengiz ±4 | deck | 0.214 | 3 | 2004, 2009, 2013 |
| Cengiz ±4 | droppers | **0.250** | 2 | 2009, 2013 |

**Answer to the sample question:** the round-5 values **0.54 / 0.23 / 0.25 are the DROPPERS
sample** (CS/imputation 0.537, stacked 0.232, Cengiz 0.250). The n-cluster figure
`fig_weight_checkability_d6_v2_ncluster.*` is the **v2 = droppers** figure, so yes, it uses
that sample. (CS-NT and CS-NYT shares coincide because the only cohort where NYT controls
would change availability, 2020, is already single-treated-unit and hence untestable.)

## Item 2 — Spearman ρ(weight, pre_g)  (n = cohorts with positive weight)

| estimator | sample | ρ | n |
|---|---|---|---|
| CS-NT | deck | −0.78 | 13 |
| CS-NT | droppers | −0.80 | 11 |
| CS-NYT | deck | −0.77 | 13 |
| CS-NYT | droppers | −0.80 | 11 |
| imputation | droppers | **−0.80** | 11 |
| stacked | deck | +0.22 | 11 |
| stacked | droppers | +0.31 | 9 |
| Cengiz ±4 | deck | −0.48 | 6 |
| Cengiz ±4 | droppers | +0.22 | 5 |

CS/CS-NYT/imputation put weight on low-pre_g cohorts (ρ ≈ −0.8); stacked/Cengiz-droppers put
it on mid-panel cohorts (ρ > 0).

## CS–imputation gap (droppers)

- **s=1 gap per cohort** (τ̂_imp(g,g) − τ̂_CS-NYT(g,g)): 2001 **0.000**, 2002 **−0.250**,
  2004 +0.114, 2008 −0.106, 2009 +0.058, 2013 −0.250, 2014 −0.087, 2017 −0.727, 2018 −0.615,
  2019 −0.305, 2020 −0.059.
- **Aggregate gap** (common imputation weights): **CS-NYT − imputation = +0.1194**
  ( = +0.0765 − (−0.0429); imputation − CS-NYT = −0.1194).
- **Per-cohort contributions** w_g·gap_avg_g: 2002 −0.071, 2017 −0.033, 2013 −0.025, 2009
  +0.016, 2004 +0.008, 2018 −0.007, others ≤|0.005|. **2002 + 2017 + 2013 = 75.2%** of the
  total |contribution|.

## Item 3 — pre-tests, 2009 / 2013 / 2017 (each estimator's own comparison group)

Slope test (differential linear pre-trend; p-values) and joint pre-trend test, under four
inference methods: CR1 cluster-robust, restricted wild cluster bootstrap (Webb, B=999 slope
/ 499 joint), randomization inference with unmatched never-treated placebos (RIu), and RI
with placebos matched on pre-period outcome mean (RIm). `MDE` = Zc · WCR bootstrap SE of the
slope.

| cohort | comparison | slope | sl CR1 | sl WCR | sl RIu | sl RIm | MDE | jt CR1 | jt WCR | jt RIu | jt RIm |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 2009 | never | −0.006 | 0.710 | 0.711 | 0.744 | 0.800 | 0.042 | 0.000 | 0.314 | 0.006 | 0.074 |
| 2009 | No-Change | −0.018 | 0.182 | 0.237 | 0.329 | 0.445 | 0.040 | 0.000 | 0.496 | 0.003 | 0.034 |
| 2009 | never+not-yet | −0.005 | 0.719 | 0.696 | 0.744 | 0.800 | 0.039 | 0.000 | 0.186 | 0.002 | 0.024 |
| 2013 | never | −0.023 | 0.019 | 0.092 | 0.131 | 0.570 | 0.043 | 0.000 | 0.538 | 0.001 | 0.022 |
| 2013 | No-Change | −0.021 | 0.014 | 0.125 | 0.160 | 0.605 | 0.041 | 0.000 | 0.568 | 0.001 | 0.001 |
| 2013 | never+not-yet | −0.021 | 0.026 | 0.105 | 0.150 | 0.605 | 0.042 | 0.000 | 0.546 | 0.001 | 0.001 |
| **2017** | never | −0.027 | 0.000 | **0.037** | **0.008** | 0.060 | 0.041 | 0.000 | 0.156 | 0.001 | 0.001 |
| **2017** | No-Change | −0.025 | 0.000 | **0.026** | **0.013** | 0.089 | 0.041 | 0.000 | 0.114 | 0.001 | 0.001 |
| **2017** | never+not-yet | −0.027 | 0.000 | **0.031** | **0.008** | 0.060 | 0.042 | 0.000 | 0.178 | 0.001 | 0.001 |

**2017 result.** Of the three testable (multi-cluster) cohorts, **2017 is the one with a
pre-trend that survives few-cluster-robust inference.** Its differential slope (≈ −0.027/yr)
is significant under CR1 (p≈0), **WCR (0.026–0.037)**, and RI-unmatched (0.008–0.013); only
RI-mean-matched softens it to 0.060–0.089 (borderline). This contrasts with **2009** (flat
slope under every method) and **2013** (CR1-significant but **not** under WCR, 0.09–0.13). As
before, the joint test rejects under CR1 and RI but not under WCR (the documented
over-rejection when restrictions ≫ treated clusters), so the slope test is the reliable read.
Bottom line: the weight that remains "testable" after the single-cluster gate still includes
a cohort (2017) that fails its own pre-trend check, and 2017 is also the 2nd-largest
contributor to the CS–imputation gap.

## Out of scope (not done, per brief)

Leave-one-out on the stacks; the recombination gap; new figures. No figure was regenerated —
no figure's numbers changed (the moved number, imputation ρ, already matched the committed
n-cluster figure, which is droppers).
