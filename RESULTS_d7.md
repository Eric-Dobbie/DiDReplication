# D7 — analytic geometry sweep (paper §3.3 / appendix)

`R/d7_geometry_sweep.R`. Fully analytic: weights are closed-form functions of geometry, the
MDE is the iid closed form. No outcomes, no estimation, no random draws. Fixed 12-geometry
grid (T∈{10,25} × timing{early,spread,late} × size{equal,concentrated}); N_c=40 never-treated,
N_T=40 treated over K=4 cohorts, σ=1, each unit a cluster. Outputs: `d7_sweep_cohort.csv`
(geometry×estimator×cohort), `d7_sweep_summary.csv` (geometry×estimator), `fig_d7_geometry_sweep.*`.

Adoption times (after dedup/shift — none needed): early T10 = 2,3,4,5 / T25 = 2,6,9,13;
spread T10 = 2,5,7,10 / T25 = 2,10,17,25; late T10 = 5,7,8,10 / T25 = 13,17,21,25.

## Which geometries discriminate CS vs stacked (gap in U(∞) = CS − stacked)

| geometry | U(∞) CS | U(∞) stacked | **gap** |
|---|---|---|---|
| T25 spread equal | 0.48 | 0.07 | **0.41** |
| T10 spread equal | 0.45 | 0.14 | **0.31** |
| T25 early equal | 0.32 | 0.06 | **0.27** |
| T25 spread concentrated | 0.21 | 0.03 | 0.18 |
| T10 early equal | 0.30 | 0.13 | 0.17 |
| T10 spread concentrated | 0.20 | 0.07 | 0.13 |
| T25 early concentrated | 0.16 | 0.04 | 0.12 |
| T10 early concentrated | 0.15 | 0.08 | 0.06 |
| all **late** geometries | 0.00 | 0.00 | **0.00** |

**The diagnostic discriminates on testability (U(∞)) in early- and spread-timed geometries,
most sharply with equal cohort sizes** — the largest gap is 0.41 (T25, spread, equal). It does
**not** discriminate in late geometries (U(∞) = 0 for every estimator). Concentrated sizes
shrink the gap: putting 25 of the 40 treated units in the 2nd-earliest cohort (which adopts
late enough to be testable) moves the dominant weight onto a checkable cohort, roughly halving
CS's U(∞) versus equal.

## Does early-heavy adoption produce the expected pattern? Yes.

With early (and spread) adoption the first cohort lands at t = 2 (pre_g = 1 < 2 → MDE = ∞), and
**CS loads weight there**: its simple weight W_g ∝ n_g·post_g rewards early adopters (large
post_g), so ρ(weight, pre_g) = −0.8 to −1.0 and U(∞) reaches 0.30–0.48 (equal sizes). The
**stacked V_s ∝ p(1−p)·q(1−q) peaks at mid-panel (q ≈ ½)**, so it places little weight on the
earliest cohort (ρ = +0.2 to +1.0) and its U(∞) is 0.04–0.14 — exactly the contrast the
diagnostic is meant to expose. **Late-heavy adoption narrows it to nothing**: with no cohort at
t = 2 (earliest adoption t = 5 or 13 → pre_g ≥ 4), every cohort is testable and U(∞) = 0 for
all three estimators.

**A reversal at the break-even tolerance.** At δ̄_g = θ*/a_g with θ* = 0.25, U ≈ 1 for *every*
estimator in early/spread geometries — the weight that is "testable" is still not *checkable*,
because the detectable slope (best MDE ≈ 0.07 here) far exceeds the break-even slope. So the
estimators separate on testability but **not** at the power-aware tolerance in early/spread —
both are essentially P&P-like (fully unchecked). The one place U at break-even separates them is
**late geometries, where CS does better than stacked**: CS's g−1 baseline gives a small
accumulation factor a_g = (post_g+1)/2 for late cohorts (few post periods), so its break-even
tolerance δ̄_g is larger and more cohorts clear it (e.g. T25 late equal: U = 0.79 CS vs 1.00
stacked; T25 late concentrated: 0.20 vs 1.00). Imputation and full-panel stacked share
a_g = T/2, so their break-even tolerance is tiny and they stay near U = 1 everywhere.

## P&P validation

- **CS weights reproduce the canonical CS-NT simple cohort weights to 2.7e-12** (≤ 1e-8). ✓
- **Stacked** closed form vs canonical FWL w_s: correlation **0.980**, L1 **0.162**, max |diff|
  0.060 — close but not exact, as expected. Attribution of the gaps:
  1. **2002** (diff −0.060): **controls definition** — the closed form uses n_c = 0 (no clean
     controls) → V_s = 0, while the real FWL gets V_s > 0 from the within-treated adopter
     contrast.
  2. **2016** (diff +0.038): a **single-treated reversal** (Portsmouth's 2nd switch) — the real
     stacked fit gives V_s = 0 (collinear, no within-stack contrast), the closed form gives
     nonzero.
  3. **2009, 2004** (diff +0.035, +0.008): **imbalance** — the closed form assumes a balanced
     common-pool stack, but P&P's real stack mixes 41 never-treated with ~700 always-treated
     control units, so the effective shares differ.

**P&P is the extreme case this sweep contrasts with:** every one of its weighted cohorts is
unchecked (U ≈ 1 at any usable tolerance), which in the sweep corresponds to the early/spread,
equal-size corner — the discriminating region — pushed to its limit.

## Caveat (state in the paper)

The MDE used here is the **iid closed form** — it assumes independent errors and ignores the
few-cluster problem (it uses the normal critical value, not a wild-bootstrap one). It is
therefore **optimistic**: real detectable slopes are larger, so every reported unchecked weight
U is a **lower bound**. The true unchecked share is at least this large.

## Validation summary (all pass)

1. P&P: CS exact (2.7e-12); stacked gap reported and attributed (above).
2. Weights sum to 1 within each geometry × estimator (max |Σw − 1| = 1.1e-16).
3. U at break-even is non-increasing in θ* and ≥ U(∞) everywhere.

## Noted but not pursued (scope capped at the grid)

- At the break-even tolerance the diagnostic barely separates CS from stacked except in late
  geometries; a finer θ* grid or an SA panel would sharpen the break-even comparison, but the
  brief caps the scope, so this is left as a note.
