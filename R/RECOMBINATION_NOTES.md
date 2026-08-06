# Recombination notes — CS, Sun-Abraham, and Stacked atom decomposition

Outcome: `any.fatalities` (probability of a fatal civilian encounter), Payson &
Parinandi (2024), fatal-encounters analysis (Table 3). Not the racial-diversity
analysis.

## What was reconstructed

Their stacked specification was recovered from `stacked_fatal.csv` and
`stacked_weights.R` (their code), not from the instruction file's defaults:

| Role | Their column |
|---|---|
| outcome | `any.fatalities` |
| treatment indicator | `no.req` (=1 once a residency requirement is absent; "requirement dropped") |
| unit FE | `agency.id` |
| stack-specific time FE | `year.cohort` (year × cohort) |
| stack / cohort id | `cohort` (= treatment year) |
| balancing weights | `weights` (entropy-balanced; treated = 1) |

- **Control group = never-treated "clean controls."** Each cohort-`g` stack =
  the agencies treated in year `g` (droppers *and* adopters) + the same pool of
  741 never-changing ("No Change") control agencies. No not-yet-treated
  borrowing.
- **Window.** `scaled.year` is defined only on ±4 and is used to build the
  entropy-balancing weights and event-study plots; it does **not** trim the
  regression. The pooled and cohort regressions run on all outcome-available
  rows (2000–2020). We replicate that (no ±4 trim).
- **Fixed effects.** The paper's m3/m4 use `agency.id + year.cohort` — a
  *shared* agency FE (not stack-specific) plus a stack-specific time FE
  (`year.cohort` = year×cohort). Our decomposition instead uses the
  **leakage-free** `agency_stack (= agency.id × cohort) + year.cohort` FE — see
  the FE-correction section below.

## Recombination results

`sum(weight_i * estimate_i)` vs each estimator's package-reported pooled value:

| Estimator | Reconstructed | Package pooled | Abs diff | Rel diff |
|---|---|---|---|---|
| CS (`aggte` simple) | 0.04476502 | 0.04476502 | 6.9e-18 | 1.6e-16 |
| SA (`sunab` agg="att") | 0.04476502 | 0.04476502 | 0 | 0 |
| Stacked (FWL, est-sample resid) | −0.10439876 | −0.10439722 | 1.5e-6 | 1.5e-5 |
| Stacked (naive, resid on ALL rows) | −0.09075827 | −0.10439722 | 1.4e-2 | **13.1%** |

CS and SA match to machine precision: the exact aggregation weights are
extracted from the fitted objects — `pg` (group shares) from `did`'s internal
one-row-per-unit analysis data, and `colSums(sign(model.matrix))` (cohort ×
period observation counts) from `fixest`.

## Reversible treatment and the CS ≡ SA equivalence

Two agencies — **memphis TN** and **portsmouth NH** — both *adopted* and later
*dropped* a residency requirement, so `year.changed` is not unit-constant. The
paper itself flags exactly these two (`callaway_replication.R` drops them). We
collapse each unit to its **first treatment year** (standard staggered-adoption
timing; `did` ≥ 2.2 errors on reversible `gname`). This reproduces the CS
estimate that `did` 2.1.2 produced implicitly, and applies the *same* cohort
definition to SA.

With a consistent cohort definition, **CS and SA are identical to machine
precision** (overall and every cohort, max diff ~4e-13). This is expected: for a
**balanced panel with a never-treated control group**, the CS "simple" and SA
"att" aggregations put the same weight (∝ cohort size) on the same 2×2 DiD
cells, so the two estimators coincide. An earlier version double-counted the two
reversible units in SA (cohort assigned row-wise), which spuriously pushed SA to
0.0517; the corrected value is 0.0448 = CS. The only figure where CS and SA
differ is the calendar-time vs. event-time view, because that reflects the
aggregation *dimension*, not the estimator.

## Where the stacked check "failed", and the diagnosed reason

The naive reconstruction (which mirrors their `stacked_weights.R` verbatim)
misses the pooled coefficient by **13.1%**. This is a **sample mismatch in the
FWL residualization**, not clustering, covariates, or a GMM weighting matrix:

- `any.fatalities` is only observed 2000–2020 (Fatal Encounters DB); the panel
  starts in 1987, so rows 1987–1999 have a missing outcome.
- Their code residualizes `no.req ~ 1 | agency.id + year.cohort` over **all
  278,324 rows**, including the 1987–1999 rows. Those rows contribute
  residualized-treatment variance to the **weights** but are dropped from the
  **coefficient** regression (which needs a non-missing outcome). Weights and
  coefficients are then computed on different samples, so the FWL identity
  `beta_pooled = sum_g w_g * beta_g` breaks.
- Residualizing on the **estimation sample** (non-missing outcome — the rows
  `felm` actually uses) restores the identity to 1.5e-6, with identical FEs,
  clustering, and no covariates. Only the residualization sample changes.

`extract_stacked_atoms(..., resid_on = "estimation")` (the default) is the
correct decomposition; `resid_on = "all"` reproduces the naive, biased version
for comparison. The driver reports both.

Note: the naive number (−0.0908) coincidentally lands on the paper's Table 3
headline (−0.0908, which is the *covariate-adjusted, entropy-weighted* spec).
That is a coincidence; the plain-TWFE pooled coefficient the FWL weights target
is −0.1044.

## Fixed-effects correction (leakage-free `agency × stack` unit FE)

The paper's stacked FE `agency.id + year.cohort` pools each agency's fixed effect
across every stack it appears in. Agency reuse is pervasive — of 776 agencies,
**743 (95.7%) appear in >1 stack**, driven entirely by the 741 never-treated
controls (each in 11 stacks); treated units appear in 1 (the 2 reversibles in 2).
The correct stacked spec interacts the unit FE with the stack:
`agency_stack (= agency.id × cohort) + year.cohort`. `fe_diagnostic.R` quantifies
the effect on the P&P full panel:

| spec | shared `agency.id` | interacted `agency×stack` | Δ |
|---|---|---|---|
| m3 (covariates) | −0.1026 | −0.0971 | +0.0055 |
| plain (no cov) | −0.1044 | −0.0989 | +0.0055 |
| m4 (cov+weights, Table 3) | −0.0908 | −0.0891 | +0.0017 |

**The pooled ATT barely moves** (~0.005 / ~0.002, well within the SE ≈ 0.038) —
a non-issue for the paper's headline, smaller than the recombination gap (0.014)
and far smaller than the window sensitivity (~0.06). But the **per-cohort effects
are artificially homogenized** under the shared FE: every stack comes out ≈ −0.11.
Under the correct interacted FE they reveal their true heterogeneity (+0.02 to
−0.37); it averages out in the pooled estimate because cohort 2009 (weight 0.37)
is stable. `run_extraction.R` therefore uses the interacted FE for the
decomposition, and the atom/cohort figures reflect that leakage-free spec (pooled
stacked = −0.0989).

## Identified stacks (under the leakage-free `agency × stack` FE)

Cohort **2009** dominates the weight (0.37, 8 treated agencies) at β ≈ −0.12; the
other cohorts are genuinely heterogeneous once the FE leakage is removed (e.g.
2018 ≈ −0.29, 2014 ≈ −0.22, 2012 ≈ +0.02). Two structural atoms:

- **2002** (weight 0.059, β = +0.007) is identified *without clean controls* —
  purely off a dropper-vs-adopter contrast within the stack. A genuinely
  different kind of atom; it pulls the pooled estimate toward zero. (Its estimate
  is unchanged by the FE correction, since it has no controls to leak across.)
- **2016** now drops out under the interacted FE: its only treated unit is the
  reversible portsmouth agency, which is collinear within its own stack once the
  unit FE is stack-specific. (Under the shared FE it carried ≈0 weight anyway.)

Cohorts 2000, 2001, 2003 drop out entirely (NA/NaN coefficient: no clean
controls and no within-stack treatment contrast).

## Reconstructing the stacked dataset from first principles

We do not have the authors' stack-construction code, so `reconstruct_stacked.R`
rebuilds `stacked_fatal.csv` from `dta.csv` to test how well the rule is
recovered. It matches **exactly**: same 278,324 rows, same cohorts, identical
row-set, and zero mismatches on `treat`, `no.req`, `any.fatalities`,
`scaled.year`, and `year.cohort`.

The recovered construction rule:
- Stacks are indexed by treatment cohort `g` (a policy CHANGE year); cohorts
  with `g >= 2000` are kept (the outcome is only observed from 2000).
- Treated rows of stack `g` = the `dta` rows with `year.changed == g`. Because
  `year.changed` is stored row-wise, the two reversible agencies (memphis,
  portsmouth) are split at their change point automatically — the pre-change
  segment goes to the adopt-cohort, the post-change segment to the drop-cohort
  (reproducing the 22/12 and 29/5 partial panels).
- Clean controls = the 741 never-changing agencies, full panels, attached only
  when `g - 4 >= 2000` (=> cohorts 2004–2020).

**The "implicit window" is ±4 (4 pre, 4 post), but it does NOT trim rows.** Every
unit keeps its full 1987–2020 panel; the window only (a) gates control
eligibility (a cohort needs 4 observable pre-periods) and (b) defines
`scaled.year` and the entropy-balancing sample. A literal ±4-windowed rebuild
(9 rows/unit) would therefore NOT match — the provided file is full-panel.

The single component `reconstruct_stacked.R` cannot reproduce is the `weights`
column: treated weights are 1, but the control weights are entropy-balancing
(`ebal`) weights that require the authors' balancing specification (covariates,
moments, normalization — they sum to ~70 per cohort).

## CS/SA note

For CS and SA the outcome window (2000–2020) forces the effective sample:
`att_gt` automatically drops 716 units "already treated in the first period"
(all always-no-requirement cities coded to 1987, plus cohorts ≤ 2000), leaving
never-treated controls (gname = 0, 41 units) + cohorts 2001–2020. The SA sample
is matched to this. This differs from the stacked design, which treats **all**
741 never-changing cities as clean controls — a genuine cross-estimator
difference in the control pool, not an error.
