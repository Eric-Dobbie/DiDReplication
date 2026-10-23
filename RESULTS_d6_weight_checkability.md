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

---

## Round 3 (items 1–3)

**1. Robust inference for the 2009/2013 pre-tests — CR1 vs restricted wild cluster
bootstrap (Webb, B=999 slope / 499 joint) vs randomization inference over never-treated
placebo cohorts (B=999).** `d6_pretests_2009_2013_robust.csv`. (fwildclusterboot is not
installable here; WCR is implemented by hand — restricted bootstrap-t for the slope,
restricted Wald bootstrap for the joint test. Cengiz ≡ CS, both never-treated.)

| cohort | comparison | slope | slope CR1 | slope WCR | slope RI | joint CR1 | joint WCR | joint RI |
|---|---|---|---|---|---|---|---|---|
| 2009 | CS/Cengiz (never) | −0.006 | 0.71 | 0.71 | 0.74 | 0.00 | 0.31 | 0.006 |
| 2009 | stacked (No-Change) | −0.018 | 0.18 | 0.24 | 0.33 | 0.00 | 0.50 | 0.003 |
| 2009 | imputation (never+not-yet) | −0.005 | 0.72 | 0.70 | 0.74 | 0.00 | 0.19 | 0.002 |
| 2013 | CS/Cengiz (never) | −0.023 | 0.019 | 0.092 | 0.131 | 0.00 | 0.54 | 0.001 |
| 2013 | stacked (No-Change) | −0.021 | 0.014 | 0.125 | 0.160 | 0.00 | 0.57 | 0.001 |
| 2013 | imputation (never+not-yet) | −0.021 | 0.026 | 0.105 | 0.150 | 0.00 | 0.55 | 0.001 |

- **The CR1 joint rejection is confirmed to be an over-rejection artifact.** CR1 gives
  joint p = 0 everywhere; the restricted wild cluster bootstrap gives **0.19–0.57 (no
  rejection)** for all six. This is precisely the B-chain finding (joint Wald over-rejects
  when restrictions ≫ treated clusters), now shown directly on these cohorts.
- **2013's "significant" slope does not survive either.** CR1 p = 0.014–0.026 rises to
  **WCR 0.09–0.13 and RI 0.13–0.16** — not significant at 5% once the 3-treated-cluster
  count is handled. 2009's slope is non-significant under every method.
- **RI-joint rejects (p ≈ 0.001–0.006), but that is a mismatched-reference artifact, not
  evidence of a pre-trend:** the placebo null draws from the 41 never-treated agencies, so
  its Wald distribution is far tighter than the real 3–7-treated-cluster sampling
  distribution, making the observed Wald look extreme. The slope RI (a single, scale-free
  coefficient) is better behaved and agrees with WCR.
- **Bottom line:** under few-cluster-correct inference (restricted WCR, the recommended
  procedure), **neither 2009 nor 2013 shows a significant pre-trend — slope or joint.** The
  round-2 CR1 rejections (2013 slope; all joint) were small-cluster artifacts. So the two
  highest-weight *testable* stacked cohorts do pass their pre-tests once inference accounts
  for the handful of treated clusters — which strengthens, not weakens, the reliance on
  those few cohorts, since the untestable/early cohorts still carry no check at all.

**2. Predicted vs. actual CS–imputation per-cohort gap (droppers).**
`d6_gap_predicted_vs_actual.csv`, `figures/d6_gap_predicted_vs_actual.*`. Predicted gap =
`slope_g × (pre_g − 1)/2`, with `slope_g` the cohort's differential linear pre-trend vs
never-treated; actual = `τ̂_imp − τ̂_CS` (cohort averages). **Correlation 0.50 (n = 9;**
cohorts 2001–2002 excluded, `pre_g < 3` so no estimable slope). The direction is right for
most cohorts and the biggest actual gaps (2017, 2018, 2013) are correctly signed but
*under-predicted* in magnitude; 2019 is the one sign flip. So an extrapolated linear
pre-trend accounts for roughly half of the per-cohort CS–imputation disagreement — real
and in the predicted direction, but not the whole story (the remainder is nonlinear
pre-trends, effect dynamics, and small-cohort noise). This is consistent with the bridge
result that the gap is a Layer-2 (counterfactual) phenomenon.

**3. Cengiz panel caption (done).** All three figures now carry: "Panel B2 (Cengiz ±4):
only cohorts 2004–2016 have a full ±4 window in 2000–2020; cohorts with pre_g < 4
(2001–2003) are EXCLUDED by the window, not down-weighted to zero."

New outputs: `d6_pretests_2009_2013_robust.csv`, `d6_gap_predicted_vs_actual.csv`,
`figures/d6_gap_predicted_vs_actual.*`.

---

## Round 4 (items 1–3)

**1. Exact CS–imputation gap at s=1 (event time 0) + own/inherited decomposition.**
`d6_gap_identity_s1.csv`, `d6_liu_weights.csv`, main figure `d6_gap_decomposition.*`.
Δ̂^Imp_{g,0} := τ̂_imp(g,g) − τ̂_CS-NYT(g,g) (imputation FE counterfactual vs Callaway–
Sant'Anna with not-yet-treated controls, varying base = g−1).

- **Verified:** cohort **2001 gap = −2.7e-11** (≈ 0 to 1e-8) — Denver has a single
  untreated pre-period, so its imputation residual at g−1 is identically zero. 2002 gap =
  **−0.250**. The own/inherited split sums to the gap **exactly** (max |own+inh−gap| = 0).
- The gap is **overwhelmingly the own last-pre-period block bias**; the inherited
  adjustment-cohort term (imputation's γ̂ being pooled over all untreated obs rather than
  cohort g's initial control group) is small everywhere — e.g. 2017: own −0.681 vs
  inherited −0.046; 2002: own −0.234 vs inherited −0.016; |inherited| ≤ 6% of the gap.
- The Liu (2025) Props 2.1/2.2 weights `N_k/(Σ_{j≥k}N_j + N_∞)` are constructed from cohort
  sizes (`d6_liu_weights.csv`, e.g. 2009 → 0.121, 2002 → 0.063) as the theoretical
  attribution of the inherited term across later (adjustment) cohorts. **Flag:** Liu (2025)
  is not available in this environment, so the total decomposition is exact by construction
  while the per-cohort inherited *attribution* follows the stated formula as a
  reconstruction — not independently verified against the paper.
- The round-3 linear-slope predicted-vs-actual plot is moved to the appendix
  (`d6_appendix_gap_linear_slope.*`); this block-bias decomposition is the main gap figure.

**2. mde_g via slope-only restricted wild cluster bootstrap (Webb, B=999), all cohorts.**
`d6_mde_wcr.csv`, figure `fig_weight_checkability_d6_v2_wcrmde.*`; `n_treated_clusters` added.
mde = Zc·SE_WCR, SE_WCR = SD of restricted bootstrap slopes.

- The CR1 mde **badly understates** the detectable slope for **single-treated-cluster**
  cohorts: e.g. CS 2018 0.007 → **0.094** (14×), 2008 0.022 → 0.082, 2019 0.006 → 0.081;
  but for the 7-cluster cohort **2009 it barely moves** (0.041 → 0.042). With one treated
  cluster the cluster-robust analytic SE is essentially meaningless; WCR restores an honest
  (much wider) interval.
- **Ranking reshuffle (checkability, low mde = more checkable):** under WCR the high-weight
  **2009 jumps from least-checkable (CR1 rank #8) to #3–4**, while the single-unit late
  cohorts (2018, 2019) become the *least* checkable. Every estimator's ordering changes for
  most cohorts (full list in the run log / `d6_mde_wcr.csv`). The qualitative picture holds
  — the top-weight untestable cohort 2002 is still off the axis — but the mde-based ordering
  among *testable* cohorts is an artifact of CR1 until WCR corrects it, and 2009's weight
  then sits on a moderately-checkable (not extreme-mde) cohort.

**3. RI diagnostic + size-matched joint RI (presented, not adjudicated).**
`d6_ri_diagnostics.csv`, `d6_ri_sizematched.csv`.

| cohort | treated mean_y / var_y / officers | never-pool mean_y / var_y / officers |
|---|---|---|
| 2009 | 0.49 / 0.25 / 521 | 0.20 / 0.16 / 950 |
| 2013 | 0.51 / 0.26 / 798 | 0.24 / 0.18 / 950 |

Treated cohorts have **higher pre-period outcome mean and variance** (and fewer officers)
than the never-treated placebo pool. Redrawing the placebo cohorts from **size-matched**
never-treated agencies does **not** close the RI rejection:

| cohort | joint p — WCR (mean) | joint RI unmatched | joint RI size-matched |
|---|---|---|---|
| 2009 | 0.332 | 0.006 | 0.002 |
| 2013 | 0.551 | 0.006 | 0.002 |

So the RI–WCR disagreement is not a size-composition artifact — it survives matching on
agency size, pointing instead to the outcome level/variance gap between treated cohorts and
the never-treated pool. WCR (resampling the actual treated clusters) does not reject; RI
(against a calmer, higher-N never-treated null) does. Both are reported; neither is declared
correct here.

New outputs: `d6_gap_identity_s1.csv`, `d6_liu_weights.csv`, `d6_gap_decomposition.*`,
`d6_mde_wcr.csv`, `fig_weight_checkability_d6_v2_wcrmde.*`, `d6_ri_diagnostics.csv`,
`d6_ri_sizematched.csv`, `d6_appendix_gap_linear_slope.*`.
