# assumptions.md

Complete record of every assumption, mapping, sample restriction, transformation,
and analysis decision made so far in this project, with rationale and the code
that implements them. Nothing about the analysis should depend on knowledge that
is not written here or in the referenced repo files.

Project: replication + estimator decomposition for **Payson & Parinandi (2024),
"Residency Blues: The Unintended Consequences of Police Residency Requirements."**
We analyze the **fatal-encounters** outcome only (their Table 3), **not** the
racial-diversity outcome.

Repo: `Eric-Dobbie/DiDReplication`, working branch
`claude/fatal-stacked-regression-decomp-1t9frt`.

---

## 1. Goal

Decompose three staggered-adoption DiD estimators into their component "atoms"
(disaggregated estimates) and aggregation weights, verify each recombines to its
own pooled estimate, and stress-test the stacked design (window length, control
group, fixed-effects specification). The three estimators:

1. **CS** — Callaway & Sant'Anna, `did::att_gt` → ATT(g,t), aggregated `simple`.
2. **SA** — Sun & Abraham, `fixest::sunab` → CATT(g,e), aggregated `att`.
3. **Stacked** — cohort sub-experiments (the paper's own approach), decomposed via
   Frisch–Waugh–Lovell (FWL) into cohort effects + variance-share weights.

---

## 2. Environment & reproducibility

- **Data is NOT in the repo.** The three CSVs (`dta.csv`, `stacked_fatal.csv`,
  `coord.csv`) come from the authors' replication zip. In this workspace they live
  at a scratch path; scripts locate them via the **`DIDREP_DATA`** environment
  variable (default `./data`). `.gitignore` excludes `/data/`. Outputs go to
  **`DIDREP_OUT`** (default `./output`), figures to **`DIDREP_FIG`** (default
  `./figures`).
- **R via conda-forge** (CRAN is blocked by the egress proxy in this environment;
  installed R + packages through micromamba/conda-forge). Interpreter used:
  a conda env's `Rscript`.
- **Package versions** (record — behavior depends on these):
  - R **4.3.3**
  - **`did` 2.5.1**  ← was 2.1.2 earlier in the session; auto-upgraded when other
    packages were installed. **This matters:** `did` ≥ 2.2 *errors* on reversible
    `gname` (see §5). `did` 2.1.2 silently used first-treatment; our recode now
    makes that explicit so results are version-stable.
  - `fixest` 0.13.2, `lfe` 3.1.1, `ggplot2` 3.5.2, `ggrepel`, `patchwork`, `broom`,
    `dplyr`, `tidyr`.
- **CS bootstrap is seeded** (`set.seed(0)` inside `extract_cs_atoms`): `att_gt`
  uses a multiplier bootstrap for SEs by default. Point estimates and weights are
  deterministic; only SEs jittered run-to-run until we seeded it. `atoms_long.csv`
  is now fully reproducible.

---

## 3. Data files & schema

### `dta.csv` — main long panel (used for CS and SA)
- **26,758 rows = 787 agencies × 34 years (1987–2020)**, balanced in rows.
- Key columns and our mapping:

| Role | Column | Notes |
|---|---|---|
| Outcome | `any.fatalities` | **binary**: any fatal civilian encounter that year. Observed **only 2000–2020** (Fatal Encounters DB); **NA for 1987–1999**. |
| Treatment | `no.req` | = 1 when the agency has **no** residency requirement. The event of interest ("requirement dropped") is `no.req` 0→1. Coefficient labeled "Requirement Dropped." |
| Unit id (string) | `agency.id` | real-world agency. |
| Unit id (numeric) | `agency.num` | `as.numeric(factor(agency.id))` — `did` needs a numeric id. |
| Time | `year` | 1987–2020. |
| Cohort / first-treatment | `year.changed` | year the policy changed; **NA for never-changers**. |
| Covariates (paper Table 3 / `m2`) | `log.pop`, `log.med.inc`, `pct.white`, `pct.white.officers.imputed` | `pct.white.officers.imputed` is missing on some rows (see §14). |
| Change type | `change.type` | `"No Change"` / `"Dropped"` / `"Adopted"`. |

### `stacked_fatal.csv` — the authors' pre-built stacked dataset (used for the stacked estimator)
- **278,324 rows, 59 columns.** This is the "P&P full panel" stack.
- Additional/derived columns and our mapping:

| Role | Column | Notes |
|---|---|---|
| Stack / sub-experiment id | `cohort` | = treatment year of that sub-experiment. |
| Stack-specific time FE | `year.cohort` | **= `as.numeric(paste0(year, cohort))`** — verified; it IS year×cohort (a clean stack×time indicator). 505 distinct values. |
| Treated-unit flag | `treat` | 1 = treated cohort member, 0 = clean control. |
| Balancing weights | `weights` | entropy-balancing (`ebal`): treated = 1, controls = balanced (sum ≈ 70 per stack). **We did not reconstruct these** — they need the authors' `ebal` call. |
| Event time | `scaled.year` | = `year - cohort`, defined **only on [-4, +4]**, else NA. |
- Carries the same `agency.id`, `year`, `no.req`, `any.fatalities`, covariates, `change.type` as `dta.csv` (verified identical values on shared keys).

---

## 4. The 2000 boundary (critical, load-bearing)

- **The outcome `any.fatalities` exists only 2000–2020** (21 years); it is NA for
  1987–1999. Every sample decision flows from this:
  - Any regression drops the 1987–1999 rows automatically (missing outcome).
  - A treatment cohort is only observable if its treatment year ≥ 2000.
  - The stacked design attaches clean controls only where a **4-year pre-window**
    is observable: `cohort − 4 ≥ 2000` ⇒ cohorts **2004–2020** get controls
    (cohorts 2000–2003 are treated-only, no controls).
  - The full panel carries the 1987–1999 rows anyway (no trimming): of the
    277,134 control rows, **105,963 (38%) are pre-2000 with a missing outcome** and
    never enter any regression.

---

## 5. Reversible-treatment units (memphis, portsmouth)

- **`"memphis tennessee"`** (adopted 2004, dropped 2009) and
  **`"portsmouth new hampshire"`** (adopted 2003, dropped 2016) both changed policy
  **twice**, so their `year.changed` is not unit-constant. The paper itself flags
  and drops exactly these two (`callaway_replication.R`).
- **Handling differs by estimator, deliberately:**
  - **CS / SA (`run_extraction.R`):** collapse each unit to its **first** treatment
    year via `ave(year.changed, agency.id, FUN = min-of-nonNA)`. Rationale: this is
    the standard staggered-adoption timing, it **reproduces `did` 2.1.2's implicit
    behavior exactly** (CS overall = +0.044765 either way), and it makes CS and SA
    use a consistent cohort assignment. Required because `did` ≥ 2.2 errors on
    reversible `gname`.
  - **Stacked (P&P file):** the provided file handles them by **splitting each
    unit's panel at its change point** (memphis: 1987–2008 in the 2004 stack,
    2009–2020 in the 2009 stack). Our reconstruction reproduces this exactly using
    the row-wise `year.changed == g` rule (§13).
  - **Cengiz windowed designs:** the two are **dropped entirely** (their two events
    are < 8 years apart, so their ±k windows overlap and would reuse rows).

---

## 6. Sample restrictions per estimator

- **CS (`att_gt`):** input is the full recoded `dta`. `att_gt` internally drops
  1987–1999 (missing outcome) and **716 units "already treated in the first period"**
  (all always-no-requirement agencies coded to 1987, plus cohorts ≤ 2000). Result:
  control group = **41 never-treated agencies** (`gname == 0`, always-had-requirement),
  treated cohorts **2001–2020**. `control_group = "nevertreated"`.
- **SA (`sunab`):** restricted to `!is.na(any.fatalities) & (year.changed == 0 |
  year.changed >= 2001)` — never-treated controls + cohorts estimable in the outcome
  window; drops the pre-2000 already-treated set to match CS's effective sample.
  Never-treated recoded to `cohort = 10000` so `sunab` treats them as controls.
- **Stacked:** uses the provided `stacked_fatal.csv` as-is (all 278,324 rows;
  regressions drop the missing-outcome rows). 11 control-bearing cohorts (2004–2020)
  + 4 treated-only early cohorts (2000–2003).

**Note on control-group asymmetry (documented, not a bug):** the stacked design uses
**all 741 never-changing ("No Change") agencies** as clean controls, whereas CS/SA
use only the **41** always-had-requirement never-treated agencies. This is because
the paper's `callaway_replication.R` codes always-**no**-requirement agencies as
`gname = 1987` (already-treated), removing them from the CS control pool, while the
stack keeps every non-changer as a control. This is a genuine cross-estimator
difference in the control pool (41 vs 741), and it — together with the weighting —
drives the sign disagreement (CS/SA positive, stacked negative).

---

## 7. Transformations applied (exhaustive)

1. `agency.num <- as.numeric(factor(agency.id))`.
2. **First-treatment collapse** of `year.changed` (reversible units) — CS/SA only.
3. **Never-treated recode:** `year.changed`: `NA & no.req==0 → 0`;
   `NA & no.req==1 → 1987`; else keep. (From `callaway_replication.R`.)
4. **SA never-treated sentinel:** `cohort_sa = ifelse(year.changed==0, 10000, year.changed)`.
5. **Stacked corrected FE id:** `agency_stack = paste(agency.id, cohort)` (§9).
6. **Cengiz stacks:** window trim to `[g-k, g+k]`, `stackunit = paste(agency.id, cohort)`,
   `stacktime = paste(year, cohort)`.
No other transformations. Outcome, treatment, and covariates are used as-is from the
CSVs.

---

## 8. Estimator specifications & headline numbers

All effects are on **P(fatal encounter)**. Positive = dropping a requirement raises
fatal encounters.

### CS / SA (never-treated, no covariates, comparable sample)
- **CS overall (`aggte` simple) = SA overall (`sunab` att) = +0.044765**, identical
  to machine precision. **Rationale for the equality:** for a *balanced panel with a
  never-treated control group*, CS-simple and SA-att place the same (cohort-size)
  weight on the same 2×2 DiD cells, so they coincide. (An earlier SA value of +0.052
  was an artifact of the reversible units being double-counted across cohorts in
  `sunab`; fixed by the first-treatment collapse.)

### CS control-group × covariate variants (`cs_control_variants.R`)
| control | covariates | overall ATT |
|---|---|---|
| never-treated | none | +0.0448 |
| never-treated | paper's 4 (`m2`) | +0.0818 |
| not-yet-treated | none | +0.0562 |
| not-yet-treated | paper's 4 | +0.0770 |
All four are positive and statistically insignificant (SE ≈ 0.06–0.09, few treated
units, ≤41-unit control pool). Adding covariates ≈ doubles the estimate; the
control-group toggle moves it modestly. With covariates, CS ≠ SA (CS uses
doubly-robust adjustment inside each ATT(g,t); SA uses global linear controls) —
they diverge up to ~0.15 at the cohort level.

### Stacked — the paper's spec vs. the correct spec
Paper's `m3`/`m4` use **`agency.id + year.cohort`** FE (shared agency FE). The
**correct** spec interacts the unit FE with the stack:
**`agency_stack (= agency.id × cohort) + year.cohort`**. `fe_diagnostic.R`:

| spec | shared `agency.id` | interacted `agency×stack` |
|---|---|---|
| plain (no cov) | −0.10440 | **−0.09892** |
| m3 (covariates) | −0.10257 | −0.09707 |
| m4 (cov + weights, Table 3) | −0.09081 | −0.08914 |

**Decision:** the decomposition (`run_extraction.R`) uses the **interacted FE**
(pooled stacked = **−0.0989**). Rationale in §9.

---

## 9. The fixed-effects decision (leakage-free `agency × stack`)

- **`year.cohort` already is year×cohort** (stack-specific time FE) — verified.
- **`agency.id` is pooled across stacks:** 776 agencies, **743 (95.7%) appear in
  >1 stack** — entirely the 741 never-treated controls (each in 11 stacks); treated
  units appear in 1 (the 2 reversibles in 2). A shared agency FE therefore leaks a
  control's level across the 11 stacks it sits in.
- **Empirical effect:** the **pooled ATT barely moves** (~0.005 plain/m3, ~0.002 m4 —
  well within SE ≈ 0.038; smaller than the recombination gap 0.014 and far smaller
  than the window sensitivity ~0.06) → a **non-issue for the paper's headline**. But
  the **per-cohort atoms are artificially homogenized to ≈ −0.11 under the shared
  FE** and reveal true heterogeneity (+0.02 to −0.37) under the interacted FE. It
  averages out because cohort **2009** (weight 0.37) is stable.
- **We therefore use the interacted FE for the decomposition** so the cohort-atom
  figures are leakage-free; the Cengiz windowed figures still reference the paper's
  actual full-panel value (−0.104) as "what P&P got," since the 0.005 difference is
  immaterial to the window story.
- Under the interacted FE, cohort **2016 drops out** (its only treated unit,
  reversible portsmouth, is collinear within its own stack) → **11 stacked atoms**
  (2002, 2004, 2008, 2009, 2012, 2013, 2014, 2017, 2018, 2019, 2020).

---

## 10. FWL weights & recombination

- **Stacked weights = Frisch–Waugh–Lovell variance shares:** residualize `no.req`
  on the FEs (`idname + stacktime`), then `weight_g = Σ_g resid² / Σ resid²`.
  Cohort effects come from the fully-interacted regression
  `y ~ no.req:factor(cohort) | idname + stacktime`.
- **Residualize on the ESTIMATION sample** (non-missing outcome), *not* all rows.
  Rationale: the paper's `stacked_weights.R` residualizes over all 278k rows
  (including pre-2000 missing-outcome rows) but the coefficient regression drops
  them — a sample mismatch that makes their `sum(w·β)` miss the pooled coefficient by
  **13%**. Restricting the residualization to the estimation sample restores the FWL
  identity to ~1e-6. Both versions are reported (`resid_on = "estimation"` default;
  `"all"` for the diagnostic row).
- **CS/SA weights are extracted from the fitted objects** (exact, not reconstructed):
  CS `pg` (group shares) from `did`'s internal one-row-per-unit data; SA
  `colSums(sign(model.matrix))` (cohort×period obs counts) from `fixest`.
- **Recombination check** (`output/recombination_check.csv`), current:
  - CS: recon = pkg = +0.044765, |diff| 6.9e-18
  - SA: recon = pkg = +0.044765, |diff| 0
  - stacked (interacted FE): recon = pkg = −0.098915, |diff| 2.4e-16
  - stacked_naive (resid on all rows): −0.08228 vs −0.09892, |diff| 1.66e-2 (the
    diagnosed sample-mismatch failure).

---

## 11. Cengiz windowed designs (`cengiz_*.R`)

Literal Cengiz–Dube–Lindner–Zipperer (2019) stacks, distinct from the P&P full panel:
- **Balanced ±k window: exactly 2k+1 rows/unit.** Eligible cohorts need the full
  window observed: `g−k ≥ 2000 AND g+k ≤ 2020`.
- **Stack-specific unit AND time FE** (`agency×cohort + year×cohort`).
- **Reversible units dropped.**
- **±4 result: −0.0615 (never-treated) ≈ −0.0616 (not-yet-treated)**, wide CI,
  insignificant. Recombination exact.
- **Control-group choice is immaterial here** (never vs not-yet): the 741
  never-treated dominate; not-yet adds only ~9 future-treated units/stack. What moves
  the estimate is the **window** (−0.062 windowed vs −0.104 full-panel).
- **Stack-specific vs shared agency FE give identical results in the windowed case**
  (controls never switch `no.req`, so are absorbed either way) — but NOT in the
  full-panel case (§9), because there covariate/time-FE demeaning differs.
- **Window sensitivity** (`cengiz_window_sensitivity.R`, never-treated): coefficient
  drifts from −0.056 (k=2) toward −0.124 (k=9), approaching the full-panel −0.104 as
  k grows — but eligible cohorts collapse from 9 to 1 (k=10 empty), so the movement is
  largely **cohort composition**, not a clean long-run trace. All windowed estimates
  are imprecise.
- **Dimensionality** (`stacked_dimensionality.R` → `output/stacked_dimensionality.tex`):
  rows by control group × window; non-monotonic in k (2k+1 rows/unit grows while
  cohorts shrink); not-yet-treated adds only a few hundred rows over never-treated.

---

## 12. Control-group definitions (precise)

- **never-treated (stacked / CS-as-coded):** stacked uses all 741 `change.type ==
  "No Change"` agencies (any direction). CS-as-coded uses only the 41 with
  `gname == 0` (always-had-requirement), because always-no-requirement non-changers
  are coded `gname = 1987`.
- **not-yet-treated (Cengiz):** for cohort `g`, the 741 never-treated **plus** units
  with `year.changed > g + k` (treated after the window), excluding anything that
  changes inside `[g−k, g+k]`, excluding reversibles.

---

## 13. Reconstruction of the stacked dataset (`reconstruct_stacked.R`)

Rebuilds `stacked_fatal.csv` from `dta.csv` with **no access to the authors'
build code**, and matches it **exactly** (278,324 rows, identical cohorts and
row-set, 0 mismatches on `treat`, `no.req`, `any.fatalities`, `scaled.year`,
`year.cohort`). Rule:
- Cohorts = distinct `year.changed` values with **`g ≥ 2000`**.
- **Treated rows of stack g = `dta` rows where `year.changed == g`** (row-wise,
  which auto-splits the reversibles at their change point).
- Controls = the 741 never-changers, **full panels**, attached only when
  **`g − 4 ≥ 2000`**.
- Derived: `treat`, `scaled.year = year − g` on [-4,4] else NA, `year.cohort =
  paste0(year, g)`.
- **The ±4 "window" gates control eligibility and defines `scaled.year`/the
  balancing sample — it does NOT trim rows.** Full 1987–2020 panels are kept, so a
  literal ±4-windowed rebuild would NOT match (far fewer rows).
- **Only unreproduced column: `weights`** (entropy balancing — needs their `ebal`
  spec).

---

## 14. Effective sample size (verified)

- File rows: 278,324. Rows with observed outcome: **171,906** = 171,171 control
  (741 × 21 × 11) + **735 treated**.
- **`felm` N:** plain = **171,906**; the paper's `m3`/`m4` (with covariates) =
  **171,444**, because **462** outcome-observed rows are dropped for missing
  `pct.white.officers.imputed`.
- So the model's effective N is **171,444** (their spec) or 171,906 (plain).
  **171,171 is the control-only count** — it omits the 735 treated rows. For `m4`
  (weighted) the Kish effective N is smaller still (not yet computed).

---

## 15. Repo inventory

### `R/` scripts (all committed)
| File | Purpose |
|---|---|
| `extract_atoms.R` | data-agnostic extractors `extract_cs_atoms`, `extract_sa_atoms`, `extract_stacked_atoms`, `fwl_stack_weights`, `build_stacks`. **Engine — embedded in §17.** |
| `run_extraction.R` | driver: hardcodes P&P column names, recodes, runs the three extractors + recombination, writes `atoms_long.csv`, `recombination_check.csv`. **Embedded in §17.** |
| `reconstruct_stacked.R` | rebuild `stacked_fatal.csv` from `dta.csv` and verify the exact match (§13). |
| `fe_diagnostic.R` | Step 1–4 FE-spec diagnostic (§9). |
| `cs_control_variants.R` | CS × {never/not-yet} × {no cov / paper cov}; coefficient plot. |
| `cengiz_stacked.R` | literal Cengiz ±4 stack (never-treated) + comparison to P&P. |
| `cengiz_notyet.R` | Cengiz ±4 with not-yet-treated controls. |
| `cengiz_window_sensitivity.R` | Cengiz coefficient vs window k = 2..10. |
| `stacked_dimensionality.R` | row counts by control group × window → `.tex`. |
| `cohort_level_variants.R` | cohort-level CS/SA under control/covariate variants. |
| `plot_weight_vs_beta.R` | combined weight-vs-β scatter (3 estimators + hull). |
| `plot_weight_vs_beta_facets.R` | small-multiples, per-estimator free weight axis. |
| `plot_cohort_level.R` | one point per cohort, 3 estimators. |
| `plot_calendar_event.R` | CS→calendar-time, SA→event-time profiles. |
| `RECOMBINATION_NOTES.md` | narrative notes (recombination, FE correction, reconstruction). |

### `output/` (all committed)
`atoms_long.csv` (531 atoms: CS 260, SA 260, stacked 11), `recombination_check.csv`,
`cs_variants_{summary,atoms}.csv`, `cengiz_stacked_{atoms,summary}.csv`,
`cengiz_notyet_atoms.csv`, `cengiz_control_groups_summary.csv`,
`cengiz_window_sensitivity.csv`, `stacked_dimensionality.{csv,tex}`.

### `figures/`
`weight_vs_beta_decomposition`, `weight_vs_beta_smallmultiples`,
`cohort_level_estimates`, `cohort_level_variants`, `calendar_and_event_time`,
`cs_variants_comparison`, `cengiz_vs_pp_stacked`, `cengiz_control_groups`,
`cengiz_window_sensitivity` (each `.png` + `.pdf`).

---

## 16. Open items / known caveats

- **`weights` (ebal) column not reconstructed** — would need the authors' `ebal`
  covariates/moments/normalization.
- **Cengiz figures still reference the paper's shared-FE full-panel −0.104** (not
  the corrected −0.099); intentional, but not yet made consistent.
- **`m4` Kish effective N** not computed.
- **CS/SA covariate divergence** characterized at the cohort level but not decomposed
  atom-by-atom.
- Data path is environment-specific (`DIDREP_DATA`); scripts assume `dta.csv` and
  `stacked_fatal.csv` sit in that one folder.

---

## 17. Core code (embedded verbatim)

The full analytical engine. Everything else is in `R/` (§15).

### `R/extract_atoms.R`

```r
# =============================================================================
# extract_atoms.R
#
# Data-agnostic extraction of the disaggregated "atoms" (component estimates)
# and their aggregation weights for three staggered-adoption DiD estimators:
#
#   1. Callaway & Sant'Anna (CS)   -- ATT(g,t)
#   2. Sun & Abraham        (SA)   -- CATT(g,e)
#   3. Stacked regression   (stacked) -- cohort-specific stack effects
#
# Every function takes column-name arguments explicitly; NO variable names from
# any particular dataset are hardcoded here. Dataset-specific preparation
# (recoding never-treated, choosing the estimation sample, etc.) lives in the
# driver script (run_extraction.R), not in these functions.
#
# Each extractor returns a tidy data frame with the common schema
#   estimator, atom_id, group, time, estimate, se, weight
# and carries the estimator's own package-reported pooled/overall estimate as
# attributes:  attr(df, "pooled_estimate"), attr(df, "pooled_se"),
#              attr(df, "pooled_label").
# The recombination check (Step 5) compares sum(weight * estimate) to that
# attribute.
# =============================================================================

suppressMessages({
  library(did)
  library(fixest)
  library(lfe)
  library(dplyr)
})

# -----------------------------------------------------------------------------
# Step 2 -- Callaway & Sant'Anna atoms
# -----------------------------------------------------------------------------
# yname, tname, idname, gname : column names (gname = first-treatment period,
#   with never-treated coded 0). xformla : covariate formula (default ~1).
# control_group, clustervars, ... : passed through to did::att_gt.
#
# Weight approach:  the "simple" aggregation in did::aggte weights each ATT(g,t)
# with t >= g by  pg_g / sum_{keepers} pg  , where pg_g = P(G = g) is the share
# of units in group g (see did:::compute.aggte, type == "simple":
#   simple.att <- sum(att[keepers] * pg[keepers]) / sum(pg[keepers]) ).
# We EXTRACT pg directly from the fitted object's internal one-row-per-unit
# analysis data (dp$data) rather than reconstructing group sizes from the raw
# input, because att_gt drops already-treated units and rows with missing
# outcomes -- so raw group counts do not match the counts did actually uses.
# This makes the recombination exact (to ~1e-16) instead of merely approximate.
# -----------------------------------------------------------------------------
extract_cs_atoms <- function(data, yname, tname, idname, gname,
                             xformla = ~1, control_group = "nevertreated",
                             clustervars = NULL, seed = 0, ...) {

  # att_gt uses a multiplier bootstrap for inference by default; seed it so the
  # reported SEs are reproducible (point estimates/weights are deterministic).
  set.seed(seed)
  mp <- did::att_gt(yname = yname, tname = tname, idname = idname, gname = gname,
                    xformla = xformla, control_group = control_group,
                    clustervars = clustervars, data = data, ...)

  agg <- did::aggte(mp, type = "simple", na.rm = TRUE)

  # --- exact simple-aggregation weights: pull pg from the fitted object -------
  dp     <- mp$DIDparams
  idata  <- dp$data
  gcol   <- dp$gname
  icol   <- dp$idname
  wcol   <- if (".w" %in% names(idata)) idata$.w else rep(1, nrow(idata))
  keepu  <- !duplicated(idata[[icol]])          # one row per unit
  ug     <- idata[[gcol]][keepu]                 # unit-level group
  uw     <- wcol[keepu]                          # unit-level weight (1 if none)

  atoms <- data.frame(group = mp$group, time = mp$t,
                      estimate = mp$att, se = mp$se,
                      stringsAsFactors = FALSE)
  atoms <- atoms[!is.na(atoms$estimate), , drop = FALSE]   # na.rm, as in did

  # pg_g = weighted share of units in group g (denominator cancels on normalize)
  pg <- sapply(atoms$group, function(g) sum(uw * (ug == g)) / sum(uw))
  keepers <- atoms$time >= atoms$group                     # post-treatment cells
  w <- rep(0, nrow(atoms))
  w[keepers] <- pg[keepers] / sum(pg[keepers])
  atoms$weight <- w

  out <- data.frame(
    estimator = "CS",
    atom_id   = sprintf("g%d_t%d", atoms$group, atoms$time),
    group     = atoms$group,
    time      = atoms$time,
    estimate  = atoms$estimate,
    se        = atoms$se,
    weight    = atoms$weight,
    stringsAsFactors = FALSE
  )
  attr(out, "pooled_estimate") <- as.numeric(agg$overall.att)
  attr(out, "pooled_se")       <- as.numeric(agg$overall.se)
  attr(out, "pooled_label")    <- "did::aggte(type='simple') overall ATT"
  out
}

# -----------------------------------------------------------------------------
# Step 3 -- Sun & Abraham atoms
# -----------------------------------------------------------------------------
# gname = cohort/first-treatment period; never-treated units must be flagged so
# they act as pure controls. Pass `never_value` = the code used for never-treated
# in gname (default 0); it is internally recoded to `ref_c` (a period far beyond
# the sample) so fixest::sunab treats those units as never treated.
#
# The atoms are the disaggregated cohort x relative-time coefficients CATT(g,e),
# read from summary(res, agg = FALSE). The overall ATT is summary(res,
# agg = "att").
#
# Weight approach:  fixest's agg = "att" (see fixest:::aggregate.fixest) sets
#   shares = colSums(sign(model_matrix[, e>=0 cols]))  (weighted by obs weights
#   if present), then shares/sum(shares). I.e. each CATT(g,e) with e >= 0 is
# weighted by the number of estimation-sample observations loading on that
# cohort x relative-time cell. We reconstruct exactly that from model.matrix().
# -----------------------------------------------------------------------------
extract_sa_atoms <- function(data, yname, idname, tname, gname,
                             never_value = 0, ref_c = 10000,
                             cluster = NULL, covars = NULL, ...) {

  d <- data
  d[["._coh"]] <- ifelse(d[[gname]] == never_value, ref_c, d[[gname]])

  # optional covariates enter as linear controls alongside the sunab term
  rhs_cov <- if (is.null(covars)) "" else paste0(paste(covars, collapse = " + "), " + ")
  fml <- stats::as.formula(sprintf(
    "%s ~ %ssunab(._coh, %s) | %s + %s", yname, rhs_cov, tname, idname, tname))
  cl  <- if (is.null(cluster)) stats::as.formula(paste0("~", idname)) else cluster
  res <- fixest::feols(fml, data = d, cluster = cl, ...)

  att <- fixest::coeftable(summary(res, agg = "att"))
  att_val <- att["ATT", 1]; att_se <- att["ATT", 2]

  ct  <- summary(res, agg = FALSE)$coeftable          # disaggregated CATT(g,e)
  cn  <- rownames(ct)
  est <- ct[, 1]; se <- ct[, 2]

  # robust parse of names like "year::-3:cohort::2009" (sub() backrefs are
  # unreliable with '::' separators here, so use regexec/regmatches)
  mm_parse <- regmatches(cn, regexec("::(-?[0-9]+):.*::([0-9]+)$", cn))
  e <- as.integer(vapply(mm_parse, function(z) z[2], character(1)))
  g <- as.integer(vapply(mm_parse, function(z) z[3], character(1)))

  # keep only the sunab cohort x period coefficients (covariate rows, if any,
  # do not parse to an (e, g) and are dropped)
  is_atom <- !is.na(e) & !is.na(g)
  e <- e[is_atom]; g <- g[is_atom]; est <- est[is_atom]; se <- se[is_atom]
  cn <- cn[is_atom]

  # exact fixest att weights: observation counts per cohort x period cell
  mmx    <- stats::model.matrix(res)
  shares <- colSums(abs(sign(mmx)))[cn]
  sel    <- (e >= 0) & !is.na(est)
  w      <- rep(0, length(est))
  w[sel] <- shares[sel] / sum(shares[sel])

  out <- data.frame(
    estimator = "SA",
    atom_id   = sprintf("g%d_e%d", g, e),
    group     = g,
    time      = e,                     # relative event time
    estimate  = as.numeric(est),
    se        = as.numeric(se),
    weight    = w,
    stringsAsFactors = FALSE
  )
  attr(out, "pooled_estimate") <- as.numeric(att_val)
  attr(out, "pooled_se")       <- as.numeric(att_se)
  attr(out, "pooled_label")    <- "fixest::sunab(agg='att') overall ATT"
  out
}

# -----------------------------------------------------------------------------
# FWL weight helper (used by the stacked extractor)
# -----------------------------------------------------------------------------
# Standard OLS/TWFE weighting result: the pooled coefficient on a treatment
# indicator equals a weighted average of stack-specific coefficients, with
# weights proportional to each stack's share of the total variance of the
# treatment indicator AFTER residualizing on the model's fixed effects
# (Frisch-Waugh-Lovell). We implement it explicitly:
#
#   1. residualize `treatname` on the fixed effects (idname + stacktime), via
#      OLS with no other covariates -- resid_treat.
#   2. weight_g = sum_g(resid_treat^2) / sum(resid_treat^2).
#
# IMPORTANT: the residualization MUST use the same rows the outcome regression
# uses. If the pooled/cohort regression drops rows with a missing outcome, those
# rows must be dropped here too, otherwise the weights are computed on a
# different sample than the coefficients and the FWL identity fails. `data`
# should therefore be the ESTIMATION sample (rows with non-missing outcome).
# -----------------------------------------------------------------------------
fwl_stack_weights <- function(data, treatname, idname, stacktime, stackid) {
  fml <- stats::as.formula(sprintf("%s ~ 1 | %s + %s", treatname, idname, stacktime))
  rt  <- lfe::felm(fml, data = data)$residuals[, 1]
  df  <- data.frame(stackid = data[[stackid]], rt2 = rt^2)
  w   <- df %>%
    dplyr::group_by(stackid) %>%
    dplyr::summarize(w_num = sum(rt2), .groups = "drop") %>%
    dplyr::mutate(weight = w_num / sum(w_num))
  stats::setNames(w$weight, as.character(w$stackid))
}

# -----------------------------------------------------------------------------
# Step 4 -- Stacked regression atoms
# -----------------------------------------------------------------------------
# Operates on an ALREADY-STACKED dataset -- i.e. one "sub-experiment" per
# treatment cohort, each combining that cohort's treated units with a clean
# control pool, stacked into one long frame with a stack id. (In the Payson &
# Parinandi data this is `stacked_fatal.csv`, built upstream with never-treated
# "clean controls" and a +/-4 balancing window; see run_extraction.R and the
# accompanying note. `build_stacks()` below reconstructs such a frame from a raw
# panel when one is not provided.)
#
#   stacktime : stack-specific time FE (e.g. year x cohort)   -> tname role
#   stackid   : sub-experiment id (e.g. treatment cohort)     -> gname role
#   treatname : 0/1 treatment indicator used in the pooled spec
#
# Atoms: stack-specific coefficients from the fully-interacted regression
#   y ~ treat:factor(stackid) | idname + stacktime         (cluster by idname).
# Weights: FWL variance shares from fwl_stack_weights(), computed on the
# estimation sample (resid_on = "estimation", the default) or, to reproduce a
# naive implementation, on all rows (resid_on = "all").
#
# Stacks whose interacted coefficient is NA/NaN (no identifying variation, e.g.
# a stack with no clean controls and no within-stack treatment contrast) are
# dropped, and the surviving weights are renormalized to sum to 1.
# -----------------------------------------------------------------------------
extract_stacked_atoms <- function(data, yname, idname, stacktime, stackid,
                                  treatname, cluster = idname,
                                  resid_on = c("estimation", "all"),
                                  window = NULL, control_group = "never_treated",
                                  ...) {
  resid_on <- match.arg(resid_on)

  # estimation sample = rows with a non-missing outcome (what felm actually uses)
  est <- data[!is.na(data[[yname]]), , drop = FALSE]

  # --- stack-specific (fully interacted) coefficients -----------------------
  fml <- stats::as.formula(sprintf(
    "%s ~ %s:factor(%s) | %s + %s | 0 | %s",
    yname, treatname, stackid, idname, stacktime, cluster))
  reg <- lfe::felm(fml, data = data)
  co  <- coef(reg); se_all <- reg$se
  nm  <- gsub(sprintf("%s:factor(%s)", treatname, stackid), "",
              names(co), fixed = TRUE)
  atoms <- data.frame(group = as.integer(nm),
                      estimate = as.numeric(co),
                      se = as.numeric(se_all[names(co)]),
                      stringsAsFactors = FALSE)

  # --- FWL variance-share weights -------------------------------------------
  wdat <- if (resid_on == "estimation") est else data
  wvec <- fwl_stack_weights(wdat, treatname, idname, stacktime, stackid)
  atoms$weight <- wvec[as.character(atoms$group)]

  # drop unidentified stacks, renormalize weights over survivors
  atoms <- atoms[!is.na(atoms$estimate) & !is.nan(atoms$estimate), , drop = FALSE]
  atoms$weight <- atoms$weight / sum(atoms$weight)

  # --- pooled target: the plain (unweighted, no-covariate) stacked coef ------
  # This is the spec the FWL decomposition targets (matches the residualization).
  pf <- stats::as.formula(sprintf("%s ~ %s | %s + %s | 0 | %s",
                                  yname, treatname, idname, stacktime, cluster))
  preg <- lfe::felm(pf, data = data)
  pooled <- coef(preg)[treatname]
  pooled_se <- preg$se[treatname]

  out <- data.frame(
    estimator = "stacked",
    atom_id   = sprintf("stack%d", atoms$group),
    group     = atoms$group,
    time      = NA_integer_,           # stack effect is not time-specific
    estimate  = atoms$estimate,
    se        = atoms$se,
    weight    = atoms$weight,
    stringsAsFactors = FALSE
  )
  attr(out, "pooled_estimate") <- as.numeric(pooled)
  attr(out, "pooled_se")       <- as.numeric(pooled_se)
  attr(out, "pooled_label")    <- sprintf(
    "pooled felm(%s ~ %s | %s + %s) [unweighted, no covariates]",
    yname, treatname, idname, stacktime)
  out
}

# -----------------------------------------------------------------------------
# build_stacks() -- reconstruct a stacked frame from a raw long panel
# -----------------------------------------------------------------------------
# Provided for the data-agnostic "build the stack" requirement. NOT used by the
# P&P driver (which reads their pre-built stacked_fatal.csv), but documents /
# implements the stack-construction logic:
#   * for each treatment cohort g, take units first-treated in g PLUS a clean
#     control pool
#       - control_group = "never_treated": units never treated (gname == never_value)
#       - control_group = "not_yet_treated": units with gname > g (+ never-treated)
#   * optionally restrict each sub-experiment to a symmetric +/-`window` around g
#   * define stack-specific time id `stacktime` = paste(time, g), a within-stack
#     treatment indicator, and stack id `stackid` = g.
# Returns the stacked data frame; feed it to extract_stacked_atoms().
# -----------------------------------------------------------------------------
build_stacks <- function(data, idname, tname, gname, yname,
                         never_value = 0, control_group = "never_treated",
                         window = NULL) {
  control_group <- match.arg(control_group, c("never_treated", "not_yet_treated"))
  cohorts <- sort(unique(data[[gname]][data[[gname]] != never_value]))
  never   <- data[data[[gname]] == never_value, , drop = FALSE]

  pieces <- lapply(cohorts, function(g) {
    treated <- data[data[[gname]] == g, , drop = FALSE]
    if (control_group == "never_treated") {
      ctrl <- never
    } else {
      nyt  <- data[data[[gname]] > g & data[[gname]] != never_value, , drop = FALSE]
      ctrl <- rbind(nyt, never)
    }
    sub <- rbind(treated, ctrl)
    if (!is.null(window)) {
      sub <- sub[abs(sub[[tname]] - g) <= window, , drop = FALSE]
    }
    sub[["._stackid"]]   <- g
    sub[["._stacktime"]] <- paste0(sub[[tname]], "_", g)
    sub[["._treat"]]     <- as.integer(sub[[gname]] == g & sub[[tname]] >= g)
    sub
  })
  do.call(rbind, pieces)
}
```

### `R/run_extraction.R`

```r
# =============================================================================
# run_extraction.R  --  driver for the atom/weight extraction
#
# Loads the Payson & Parinandi (2024) replication data, supplies the ACTUAL
# column names (hardcoded here, not in the functions), calls the three
# data-agnostic extractors in extract_atoms.R, runs the recombination check,
# and writes two CSVs:
#     output/atoms_long.csv          -- all atoms + weights, all three estimators
#     output/recombination_check.csv -- reconstructed vs package pooled estimate
#
# Outcome of interest: any.fatalities (probability of a fatal civilian
# encounter). This is the fatal-encounters analysis (Table 3), NOT the racial
# diversity analysis.
#
# Data location: set DATA_DIR to the folder holding dta.csv and
# stacked_fatal.csv (defaults to ./data, override with the DIDREP_DATA env var).
# =============================================================================

# locate extract_atoms.R next to this script (works under Rscript and source())
.this_file <- tryCatch({
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) normalizePath(f) else sys.frame(1)$ofile
}, error = function(e) NULL)
.here <- if (!is.null(.this_file)) dirname(.this_file) else "R"
source(file.path(.here, "extract_atoms.R"))

DATA_DIR <- Sys.getenv("DIDREP_DATA", unset = "data")
OUT_DIR  <- Sys.getenv("DIDREP_OUT",  unset = "output")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- P&P column names (the ONLY place these are hardcoded) ------------------
Y      <- "any.fatalities"                # outcome
TIME   <- "year"                          # calendar time
UNIT   <- "agency.id"                     # unit id (string)
UNITN  <- "agency.num"                    # numeric unit id (needed by did)
COHORT <- "year.changed"                  # first-treatment (drop) year; NA = never
# stacked file:
S_TREAT <- "no.req"                       # treatment indicator (requirement absent)
S_TIME  <- "year.cohort"                  # stack-specific time FE (year x cohort)
S_STACK <- "cohort"                       # sub-experiment / stack id

# =============================================================================
# Load + diagnostics
# =============================================================================
dta <- read.csv(file.path(DATA_DIR, "dta.csv"))
dta[[UNITN]] <- as.numeric(factor(dta[[UNIT]]))

# Enforce irreversible treatment: two agencies (memphis TN, portsmouth NH) both
# ADOPTED and later DROPPED a requirement, so year.changed is not unit-constant.
# The paper flags exactly these two (callaway_replication.R drops them), and
# did >= 2.2 now errors on reversible gname. Collapse each unit to its FIRST
# treatment year (standard staggered-adoption timing) -- this reproduces the
# behavior of did 2.1.2 for CS and applies the same definition to SA, so the two
# estimators use a consistent cohort assignment.
dta[[COHORT]] <- ave(dta[[COHORT]], dta[[UNIT]],
                     FUN = function(x) { v <- x[!is.na(x)]; if (length(v)) min(v) else NA_real_ })

# P&P coding of the cohort variable for CS/SA (from callaway_replication.R):
#   never-treated (always has a requirement)  -> 0   (clean controls)
#   always-no-requirement (treated pre-sample)-> 1987 (dropped by att_gt: no pre)
dta[[COHORT]] <- ifelse(is.na(dta[[COHORT]]) & dta[["no.req"]] == 0, 0,
                 ifelse(is.na(dta[[COHORT]]) & dta[["no.req"]] == 1, 1987,
                        dta[[COHORT]]))

cat("== dta.csv ==  rows:", nrow(dta),
    " units:", length(unique(dta[[UNIT]])),
    " years:", paste(range(dta[[TIME]]), collapse = "-"), "\n")
cat("   outcome non-missing years:",
    paste(range(dta[[TIME]][!is.na(dta[[Y]])]), collapse = "-"), "\n")

# =============================================================================
# Step 2 -- CS atoms   (att_gt auto-drops already-treated + missing-outcome rows)
# =============================================================================
cs <- extract_cs_atoms(dta, yname = Y, tname = TIME, idname = UNITN,
                       gname = COHORT, xformla = ~1,
                       control_group = "nevertreated", clustervars = UNITN)

# =============================================================================
# Step 3 -- SA atoms
# Comparable sample: never-treated controls + cohorts estimable within the
# outcome window (first-treated >= 2001); drop units treated before the sample
# (1987..2000), exactly the set att_gt drops as "already treated in first period".
# =============================================================================
d_sa <- dta[!is.na(dta[[Y]]) &
              (dta[[COHORT]] == 0 | dta[[COHORT]] >= 2001), , drop = FALSE]
sa <- extract_sa_atoms(d_sa, yname = Y, idname = UNITN, tname = TIME,
                       gname = COHORT, never_value = 0, ref_c = 10000)

# =============================================================================
# Step 4 -- Stacked atoms  (their pre-built clean-control stack)
# =============================================================================
stacked <- read.csv(file.path(DATA_DIR, "stacked_fatal.csv"))   # P&P full panel = reconstruct_stacked.R output
cat("== stacked_fatal.csv ==  rows:", nrow(stacked),
    " stacks:", length(unique(stacked[[S_STACK]])), "\n")

# LEAKAGE-FREE (correct) stacked FE: interact the unit FE with the stack, so a
# control agency (which appears in 11 stacks) is demeaned separately within each
# stack -- the fully-interacted Cengiz FE (agency x stack + year x cohort). The
# paper's m3/m4 instead use a shared agency.id FE, which leaks levels across
# stacks and artificially homogenizes the per-cohort effects (see fe_diagnostic.R
# and RECOMBINATION_NOTES.md). Clustering stays on the real agency.id.
STACKUNIT <- "agency_stack"
stacked[[STACKUNIT]] <- paste(stacked[[UNIT]], stacked[[S_STACK]], sep = "__")

st <- extract_stacked_atoms(stacked, yname = Y, idname = STACKUNIT,
                            stacktime = S_TIME, stackid = S_STACK,
                            treatname = S_TREAT, cluster = UNIT,
                            resid_on = "estimation")

# also compute the naive "residualize on all rows" version for the diagnostic
st_naive <- extract_stacked_atoms(stacked, yname = Y, idname = STACKUNIT,
                                  stacktime = S_TIME, stackid = S_STACK,
                                  treatname = S_TREAT, cluster = UNIT,
                                  resid_on = "all")

# =============================================================================
# Step 5 -- Recombination check
# =============================================================================
recomb_row <- function(df, label = attr(df, "estimator_label")) {
  recon  <- sum(df$weight * df$estimate)
  pooled <- attr(df, "pooled_estimate")
  data.frame(
    estimator            = df$estimator[1],
    pooled_label         = attr(df, "pooled_label"),
    reconstructed_pooled = recon,
    package_pooled       = pooled,
    abs_diff             = abs(recon - pooled),
    rel_diff             = abs(recon - pooled) / abs(pooled),
    stringsAsFactors = FALSE
  )
}

check <- rbind(recomb_row(cs), recomb_row(sa), recomb_row(st))

# extra diagnostic line: naive stacked residualization (all rows)
naive_recon <- sum(st_naive$weight * st_naive$estimate)
check_naive <- data.frame(
  estimator            = "stacked_naive_residALLrows",
  pooled_label         = attr(st, "pooled_label"),
  reconstructed_pooled = naive_recon,
  package_pooled       = attr(st, "pooled_estimate"),
  abs_diff             = abs(naive_recon - attr(st, "pooled_estimate")),
  rel_diff             = abs(naive_recon - attr(st, "pooled_estimate")) /
                           abs(attr(st, "pooled_estimate")),
  stringsAsFactors = FALSE)
check_full <- rbind(check, check_naive)

cat("\n================ RECOMBINATION CHECK ================\n")
print(check_full, row.names = FALSE, digits = 8)

# =============================================================================
# Step 6 -- Output
# =============================================================================
atoms_long <- rbind(cs, sa, st)
write.csv(atoms_long, file.path(OUT_DIR, "atoms_long.csv"), row.names = FALSE)
write.csv(check_full, file.path(OUT_DIR, "recombination_check.csv"), row.names = FALSE)
cat("\nWrote:", file.path(OUT_DIR, "atoms_long.csv"), "(", nrow(atoms_long), "atoms )\n")
cat("Wrote:", file.path(OUT_DIR, "recombination_check.csv"), "\n")
```

### Reconstruction rule (`R/reconstruct_stacked.R`, core loop)

```r
cohorts  <- sort(unique(dta$year.changed[!is.na(dta$year.changed) & dta$year.changed >= 2000]))
nochange <- dta[dta$change.type == "No Change", ]          # 741 never-changers
parts <- list()
for (g in cohorts) {
  tr <- dta[!is.na(dta$year.changed) & dta$year.changed == g, ]   # treated segment (row-wise)
  tr$cohort <- g; tr$treat <- 1L
  parts[[length(parts)+1L]] <- tr
  if (g - 4 >= 2000) {                                            # control eligibility
    ct <- nochange; ct$cohort <- g; ct$treat <- 0L
    parts[[length(parts)+1L]] <- ct
  }
}
recon <- do.call(rbind, parts)
recon$scaled.year <- ifelse(abs(recon$year - recon$cohort) <= 4, recon$year - recon$cohort, NA_real_)
recon$year.cohort <- as.numeric(paste0(recon$year, recon$cohort))
```

### Cengiz windowed build (`R/cengiz_notyet.R`, core; `k` = half-width)

```r
REVERSIBLE <- c("memphis tennessee","portsmouth new hampshire")
dta <- dta[!dta$agency.id %in% REVERSIBLE, ]
cohorts  <- allcoh[allcoh - k >= 2000 & allcoh + k <= 2020]      # balanced-window eligibility
nochange <- dta[dta$change.type == "No Change", ]
for (g in cohorts) {
  tr <- dta[!is.na(dta$year.changed) & dta$year.changed == g &
              dta$year >= g - k & dta$year <= g + k, ]           # 2k+1 rows/unit
  ct <- nochange
  if (control == "not-yet-treated")
    ct <- rbind(ct, dta[!is.na(dta$year.changed) & dta$year.changed > g + k, ])
  ct <- ct[ct$year >= g - k & ct$year <= g + k, ]
  # ... rbind(tr, ct); stackunit = paste(agency.id, cohort); stacktime = paste(year, cohort)
}
# decompose with idname = stackunit, stacktime = stacktime, cluster = agency.id
```

### FE diagnostic (`R/fe_diagnostic.R`, core)

```r
st$agency_stack <- paste(st$agency.id, st$cohort, sep = "__")
# shared FE (paper):     any.fatalities ~ no.req + <covs> | agency.id     + year.cohort | 0 | agency.id
# interacted FE (correct): any.fatalities ~ no.req + <covs> | agency_stack + year.cohort | 0 | agency.id
# covs = log.pop + log.med.inc + pct.white + pct.white.officers.imputed
```
