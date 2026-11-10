# D6 round 6 — checkability threshold pre_g ≥ 2 and the unchecked-weight diagnostic

`R/d6_round6.R`. Checkability rule changes pre_g ≥ 3 → **pre_g ≥ 2**; the diagnostic becomes
**unchecked weight** `U^k(δ̄) = Σ_g W^k_g · 1{MDE^k_g > δ̄}` (paper §3.2). No new estimators —
committed weights/estimates and the round-4 WCR machinery are reused. MDE is the WCR-based
minimum detectable pre-trend slope: `MDE^k_g = (z_.975+z_.80)·SE_WCR(δ̂^k_g)` if pre_g ≥ 2, a
comparison group exists, and n_treated_clusters ≥ 2; else ∞.

## Headline

**Unchecked weight U(∞) — old (pre_g ≥ 3) vs new (pre_g ≥ 2):**

| estimator | sample | U(∞) old | U(∞) new | Δ = weight on newly testable |
|---|---|---|---|---|
| CS-NT | deck | 0.503 | **0.213** | 0.290 (2002) |
| CS-NT | droppers | 0.537 | **0.243** | 0.293 (2002) |
| CS-NYT | deck | 0.503 | **0.213** | 0.290 (2002) |
| CS-NYT | droppers | 0.537 | **0.243** | 0.293 (2002) |
| imputation | droppers | 0.537 | **0.243** | 0.293 (2002) |
| stacked | deck | 0.241 | 0.241 | 0 |
| stacked | droppers | 0.232 | 0.232 | 0 |
| Cengiz ±4 | deck | 0.214 | 0.214 | 0 |
| Cengiz ±4 | droppers | 0.250 | 0.250 | 0 |

Only CS-NT / CS-NYT / imputation move, and only because **2002 (pre_g = 2) becomes testable**;
its weight (≈ 0.29) leaves the untestable set. Stacked and Cengiz are unchanged: stacked 2002
has no clean controls and Cengiz's ±4 window excludes 2002, so both keep MDE = ∞ there.

**Unchecked weight at the break-even tolerance** `δ̄_g = θ*/a_g` (the slope whose extrapolated
bias equals θ*):

| estimator | sample | U(θ*=0.05) | U(θ*=0.10) |
|---|---|---|---|
| CS-NT | deck | **1.000** | 0.973 |
| CS-NT | droppers | **1.000** | **1.000** |
| CS-NYT | deck | **1.000** | 0.973 |
| CS-NYT | droppers | **1.000** | **1.000** |
| imputation | droppers | **1.000** | **1.000** |
| stacked | deck/droppers | **1.000** | **1.000** |
| Cengiz ±4 | deck/droppers | **1.000** | **1.000** |

**At θ* = 0.05, no cohort anywhere is checked — U = 1 for every estimator × sample.** At
θ* = 0.10 only one cohort checks (2018, deck, under CS), leaving U = 0.973. The reason is the
gap between what the pre-test can resolve (best MDE ≈ 0.039) and the break-even slope
(δ̄_g ≤ 0.04, and ≤ 0.0095 for imputation/stacked where a_g = 10.5): a slope small enough to
rule out matters, but too small to detect. Testability (U(∞) = 0.21–0.54) overstates how much
weight is actually checkable; the power-aware diagnostic puts it at essentially zero.

## 2002 MDE and its reliability flag

| estimator | sample | 2002 MDE (WCR) |
|---|---|---|
| CS-NT | deck / droppers | 0.786 / 0.962 |
| CS-NYT | deck / droppers | 0.768 / 0.953 |
| imputation | droppers | 0.953 |
| stacked | deck / droppers | ∞ (no clean controls — confirmed) |

**Reliability flag.** 2002 has **4 treated clusters** in droppers (5 in deck). Total Webb draws
are 6^G over all G clusters (effectively unlimited), but the bootstrap's variation is driven by
the treated side: only **6^4 = 1296** distinct treated-cluster weight patterns. With < ~10
treated clusters the restricted WCR is unreliable (MacKinnon–Webb; severe below 5), so 2002's
MDE (≈ 0.77–0.96) is reported but **should not be trusted as a precise number**. Either way it
is enormous — two pre-periods give almost no power — so 2002 is "testable" by rule yet unchecked
at any usable tolerance.

## Break-even bias a_g (generic = closed form, flag: ok)

`a_g = mean(post event times) − mean(baseline event times)`, computed generically from each
estimator's baseline; it matches the closed forms exactly for every cohort (flag `ok`):

- CS-NT, CS-NYT (g−1 baseline): `a_g = (post_g+1)/2` — e.g. 2009 = 6.5, 2013 = 4.5, 2017 = 2.5.
- imputation, full-panel stacked (pre-average baseline): `a_g = (pre_g+post_g)/2 = 21/2 = 10.5`
  for **every** cohort (the 2000–2020 panel is balanced, so pre_g+post_g ≡ 21). Worth noting:
  the pre-average baseline makes a unit slope's bias constant and large (×10.5).
- Cengiz ±4 (window [g−4, g+4]): `a_g = k + 1/2 = 4.5` for every cohort.

## Validation (all pass)

1. **MDE reproduction.** Cohorts already testable under pre_g ≥ 3 reproduce `d6_mde_wcr.csv`
   **exactly** (n = 11 droppers cells, max |diff| = 4.9e-17).
2. **U curve** is non-increasing in δ̄ and ≥ U(∞) at every grid point (δ̄ ∈ [0, 0.06], step 0.001).
3. New U(∞) ≤ old for every cell, and the difference **equals the weight on newly testable
   cohorts** (2002) to 1e-9.
4. Weights still sum to 1 (max |Σw − 1| = 9.0e-12); ATTs unchanged.

**One edge case flagged (not a failure):** stacked **deck 2004** (n_g = 2, w = 0.057) is
rule-testable (pre_g = 4, has controls) but its WCR SE is uncomputable (NA). It is counted
**testable** for U(∞) (matching the canonical untestable_share 0.241), but since its MDE is
undefined it is **never "checked"** in U(δ̄) — i.e. it stays in the unchecked set at all finite
δ̄. No other cell is affected.

## Canonical numbers that changed (full log in `d6_round6_changelog.csv`)

Only the checkability-rule-dependent quantities moved; **ATTs, SEs, ρ, and all gap numbers are
unchanged.**

| key | old | new | cause |
|---|---|---|---|
| untestable_share CS-NT\|deck | 0.503 | 0.213 | rule pre_g≥3 → ≥2 (2002 testable) |
| untestable_share CS-NT\|droppers | 0.537 | 0.243 | " |
| untestable_share CS-NYT\|deck | 0.503 | 0.213 | " |
| untestable_share CS-NYT\|droppers | 0.537 | 0.243 | " |
| untestable_share imputation\|droppers | 0.537 | 0.243 | " |
| n_testable CS-NT\|deck | 5 | 6 | +2002 |
| n_testable CS-NT\|droppers | 3 | 4 | +2002 |
| n_testable CS-NYT\|deck | 5 | 6 | +2002 |
| n_testable CS-NYT\|droppers | 3 | 4 | +2002 |
| n_testable imputation\|droppers | 3 | 4 | +2002 |

New canonical sections: `U_inf_old`, `U_inf_new`, `U_breakeven` (θ*∈{0.05,0.10}), `MDE_2002`,
`a_g`.

## Figure 1 refresh

`fig_weight_checkability_d6_v2_ncluster.*` regenerated with pre_g ≥ 2: **2002 moves from the
"not testable" strip onto the MDE axis** (blue) in the CS and imputation panels, at its CR1
MDE (≈ 0.74), far to the right; untestable weight drops from 0.54 to 0.24 there. Stacked and
Cengiz 2002 remain in the strip. Everything else is unchanged.

## Deliverables / outputs

`R/d6_round6.R`; updated `output/d6_canonical_numbers.csv`; `output/d6_unchecked_weight.csv`
(items 2–3), `output/d6_unchecked_percohort.csv` (per-cohort W_g/pre_g/post_g/a_g/δ̄_g/MDE/checked
behind item 3), `output/d6_unchecked_curve.csv` (item 4, CSV only), `output/d6_round6_mde.csv`
(all MDEs), `output/d6_round6_changelog.csv`; refreshed Figure 1.
