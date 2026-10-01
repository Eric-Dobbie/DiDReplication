# Remaining-runs results (Tasks 1, 2, 4) — status

**Environment note.** The session container was reset: the raw data (`dta.csv`,
`stacked_fatal.csv`, both gitignored) and the R toolchain were wiped. Tasks that
are **pure recombination of the committed atoms/weights** were completed from the
committed CSVs (`output/atoms_long.csv`, `output/fwl_decomp_unweighted.csv`) with a
stdlib computation and are reproduced by committed R scripts. Tasks that need a
**re-fit** (standard errors from `aggte`, CS leave-one-out SEs, dynamic/calendar
aggregates, the pre-trend figure, and Task 3) are **blocked pending re-upload of the
replication data + an R rebuild**. Point estimates throughout are exact (CS's
never-treated-control `ATT(g,t)` are independent across cohorts, and the stacked
decomposition is committed).

Common cohort set (the 11 the stacked spec identifies, `w_s > 0`): **2002, 2004,
2008, 2009, 2012, 2013, 2014, 2017, 2018, 2019, 2020**.

---

## Task 1 — Counterfactual reweighting  ✅ complete
`R/task1_counterfactual_reweight.R` → `output/task1_reweight_2x2.csv`,
`output/task1_cohort_inputs.csv`.

**2×2 (restricted to the common 11 cohorts):**

| Atoms | Weights | Estimate |
|---|---|---|
| CS | CS simple | **+0.056** |
| CS | stacked FWL | +0.075 |
| stacked | CS simple | −0.061 |
| stacked | stacked FWL | **−0.099** |

Diagonal check: stacked (restricted) = **−0.099** = unrestricted stacked (the
stacked spec only identifies these 11, so restriction is a no-op for it). CS
(restricted) = **+0.056** vs **+0.045** unrestricted (all 13 CS cohorts) — the
restriction to the common set raises CS by +0.011, reported side by side rather
than absorbed.

**Gap decomposition** (stacked_restr − CS_restr = **−0.155**):

| Channel | Value |
|---|---|
| weight channel  `Σ(W_st−W_cs)·A_cs` | **+0.019** |
| atom channel    `Σ W_cs·(A_st−A_cs)` | **−0.117** |
| interaction     `Σ(W_st−W_cs)(A_st−A_cs)` | **−0.057** |
| sum | −0.155 ✓ |

**Which channel dominates (one paragraph).** The disagreement between CS (+0.056)
and stacked (−0.099) on the common cohort set is overwhelmingly an **atom-level**
disagreement, not a weighting artifact. Holding weights fixed at CS's and swapping
only the per-cohort effects (CS → stacked) moves the aggregate by **−0.117**, three-
quarters of the total −0.155 gap; holding atoms fixed and swapping only the weights
moves it by just **+0.019**. The interaction (−0.057) is itself large because the
two estimators disagree *most* on exactly the cohorts whose weights also differ most
(2002, 2017, 2018). The per-cohort inputs make the atom gap concrete: CS and stacked
estimate different objects per cohort — e.g. 2002 is +0.222 under CS (41 never-
treated clean controls) but +0.007 under stacked (within-treated, no clean
controls); 2017 is +0.722 vs −0.059; 2018 is +0.230 vs −0.292. So the CS-vs-stacked
disagreement cannot be reweighted away: it lives in the atoms.

---

## Task 2 — Aggregation-scheme sensitivity  ◐ partial (point estimates + weights; SEs blocked)
`output/task2_cs_aggregation_partial.csv`. `simple` and `group` are exact from the
committed CS atoms; **`dynamic`, `calendar`, and all SEs require re-fitting `aggte`**
(blocked).

| Scheme | CS overall | 2002 weight | 2009 weight |
|---|---|---|---|
| simple (size × #post) | **+0.045** | 0.290 | 0.256 |
| group (size only) | **+0.106** | 0.167 | 0.233 |
| dynamic | *pending refit* | | |
| calendar | *pending refit* | | |

`simple` overweights early cohorts (2002 has 19 post periods → weight 0.29, the
single largest); `group` shifts weight onto the noisy, large-positive late cohorts
(2017 +0.72, 2018 +0.23, 2013 +0.30 each go from ~0.03 to 0.10), pushing CS from
+0.045 to +0.106. **The CS sign does not flip between simple and group — both are
positive.** The sign fragility is cohort-specific (2002), not scheme-specific (see
Task 4.1), which is the sharper statement for the main text.

---

## Task 4.1 — Leave-one-cohort-out  ✅ point estimates (CS SEs blocked)
`output/task4_cs_loo.csv` (CS); `output/stacked_loo_both.csv` (stacked, committed).

- **CS: dropping cohort 2002 flips the sign** (+0.045 → **−0.028**, Δ −0.072). No
  other cohort flips it; dropping 2009 raises it most (+0.119). CS's positive
  aggregate rests entirely on 2002.
- **Stacked: no single cohort flips the sign** — the plain corrected-FE aggregate
  stays in [−0.115, −0.082] across all leave-one-out runs (most influential: drop
  2018 → −0.082, Δ +0.017; drop 2013 → −0.115). The stacked negative is robust where
  the CS positive is not.

## Task 4.2 — Pre-trend × weight joint plot  ⧗ blocked (figure)
The stacked per-cohort pre-trend statistics are committed (`stacked_pretrends.csv`)
and cohort weights are available, but the figure (weight on x, pre-trend stat on y,
one series per estimator) and the joint CS/SA pre-trend tests need R/plotting and a
re-fit for the CS/SA joint Wald statistics. Deferred with Task 2's SEs.

---

## Blocked, pending data re-upload + R rebuild
1. Task 2 SEs and the `dynamic`/`calendar` aggregates (`aggte` re-fit).
2. CS leave-one-out SEs (att_gt re-fit per dropped cohort).
3. Task 4.2 pre-trend × weight figure (plot + CS/SA joint pre-trend tests).
4. **Task 3** boundary-asymmetry simulation — also gated on Task 1 review per the
   brief, so held regardless.

To unblock, re-upload the replication data (`dta.csv`, `stacked_fatal.csv`); the R
environment will be rebuilt as before.
