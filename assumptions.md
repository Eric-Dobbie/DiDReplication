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
- **These two are a subset of a larger issue:** 6 of 37 treated stack-members are
  *adopters* (`no.req` 1→0), not droppers, because stacks are keyed on the calendar
  year of change, not its direction. See **§18** — the coding handles the sign
  correctly (the regressor is the state `no.req`), but the adopters are load-bearing
  for the clean-control-free early cohorts.

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

**Robustness — stacked on never-treated controls only (`stacked_nevertreated.R`,
`output/stacked_nevertreated.csv`).** Re-running the full-panel stacked regression
with the control pool restricted to the **41 strictly never-treated** agencies
(dropping the 700 always-treated non-changers), corrected interactive FE:

| control pool | plain | m3 (covariates) | N |
|---|---|---|---|
| never-treated only (41) | **−0.0933** (0.0364) | **−0.1308** (0.0368) | 10,206 |
| full No-Change pool (741) | −0.0989 (0.0394) | −0.0971 (0.0394) | 171,444/171,906 |

- **The plain estimate barely moves (−0.0989 → −0.0933, Δ = +0.006, within the
  SE)** even though the 700 always-treated units are 94% of the control pool and
  ~161k of the ~172k rows. The negative stacked result is **not** an artifact of
  using always-treated units as controls — it survives on the 41 clean never-treated
  controls alone. The SE is essentially unchanged (identification comes from the
  treated switchers, not the control-pool size).
- **With covariates the pools diverge** (m3: −0.097 → −0.131): covariate adjustment
  interacts with the control composition, giving a more negative conditional estimate
  on the clean never-treated pool.
- **2002 is unchanged (+0.007)** — it has no controls in either spec, always
  identified by its within-treated drop-vs-adopter contrast (§18); same 11 identified
  cohorts.
- **Weighting not reported:** the entropy-balancing weights were built for the 741-
  control pool and do not validly rebalance the 41-control pool (would need `ebal`
  re-run); only unweighted plain and m3 are shown.

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

### Reproducing the paper's table (`spec_grid.R`, `output/spec_grid.{tex,csv}`)
The published fatal-encounters table (paper p.25; `police_residency_main.R` m1–m4,
covariates = **`log.pop + log.med.inc + pct.white + pct.white.officers.imputed`**)
reports for "Requirement Dropped": col 1 −0.097 (0.037), col 2 −0.093 (0.036),
**col 3 −0.103 (0.039)**, **col 4 −0.091 (0.037)**. Note **col 3's SE is 0.039, not
0.037** (an earlier note/prompt had 0.037 — corrected here against the PDF). A
2×2×2 grid (FE × covariates × weights, SE clustered on `agency.id`) reproduces both
stacked columns **exactly** in the **additive-`agency.id` FE + 4-controls** row:
- **col 3** = additive FE, 4 controls, **unweighted** → −0.1026 (0.0393), N 171,444.
- **col 4** = same cell, **weighted** → −0.0908 (0.0373), N 171,444.

Both require the **shared `agency.id` FE**; the leakage-free `agency×stack`
counterparts are −0.0971 (0.0400) and −0.0891 (0.0375) — close but distinct
(swapping FE moves col 3 by +0.005, col 4 by +0.002). No other cell matches either
target. **Implementation note:** felm's weighted `$se` is the *non-clustered*
analytic SE (~0.004); the clustered SE is `$cse` — `spec_grid.R` reports `$cse`.

### The full 15-cohort stack "as described", and a Table 5 / Table A.5 inconsistency
(`stacked_full15.R` → `output/stacked_full15.{tex,csv}`, `stacked_full15_regtable.tex`)

The paper *describes* (p.17) the stack as each change-cohort "along with **all pure
control cities**", no window. Building that — every change-cohort (≥2000) given all
741 controls over the full panel, so controls are **added to the 2000–2003 stacks**
the shipped file leaves treated-only — gives **N = 233,520 with covariates**, which
**exactly matches the paper's own Table A.5 col 3** (233,520). On this full stack,
unweighted, clustered on `agency.id`:

| FE | covariates | coef | cluster SE | N |
|---|---|---|---|---|
| interacted (agency×stack) | 4 P&P | **−0.092** | 0.039 | 233,520 |
| interacted (agency×stack) | none | −0.093 | 0.038 | 234,150 |
| additive (agency.id) | 4 P&P | −0.097 | 0.037 | 233,520 |
| additive (agency.id) | none | −0.099 | 0.037 | 234,150 |

The "as described" cell (interacted FE + controls on all cohorts) is **−0.092**,
vs the published Table 5 col 3 **−0.103**. The two differ on two axes, each pulling
toward zero: adding controls to 2000–2003 (additive FE: −0.103 → −0.097) and
switching to the interacted FE (−0.097 → −0.092). 13 of 15 cohorts identify (2001–
2003 now do; 2000 has no observed pre-period, 2016 is absorbed).

**The inconsistency (documented, not fixed):** the published fatal-encounters
results are not the design the paper describes.
- **FE:** the prose/label say "City and Cohort × Year FEs" and the design implies a
  Cengiz-style interacted unit FE, but the code (`police_residency_main.R` m3/m4)
  uses the **additive `agency.id`** FE (§9).
- **Sample:** **Table 5** (binary outcome) col 3 uses **N = 171,444** — the shipped
  stack, with controls only on 2004–2020, so 2000–2003 are treated-only. But the
  paper's **Table A.5** (count outcome) col 3 uses **N = 233,520** — the full stack
  with controls on *all* cohorts. So the paper already ran the "controls everywhere"
  construction for the count outcome but not for the binary Table 5; applying it to
  the binary outcome is what yields −0.092. (Both tables' col 4, weighted, are
  171,444, since the ebal weights exist only for the 11 control-bearing cohorts — §20.)

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
| `extract_atoms.R` | data-agnostic extractors `extract_cs_atoms`, `extract_sa_atoms`, `extract_stacked_atoms`, `fwl_stack_weights`, `build_stacks`. **Engine — embedded in §21.** |
| `run_extraction.R` | driver: hardcodes P&P column names, recodes, runs the three extractors + recombination, writes `atoms_long.csv`, `recombination_check.csv`. **Embedded in §21.** |
| `reconstruct_stacked.R` | rebuild `stacked_fatal.csv` from `dta.csv` and verify the exact match (§13). |
| `fe_diagnostic.R` | Step 1–4 FE-spec diagnostic (§9). |
| `cs_control_variants.R` | CS × {never/not-yet} × {no cov / paper cov}; coefficient plot. |
| `cengiz_stacked.R` | literal Cengiz ±4 stack (never-treated) + comparison to P&P. |
| `cengiz_notyet.R` | Cengiz ±4 with not-yet-treated controls. |
| `cengiz_window_sensitivity.R` | Cengiz coefficient vs window k = 2..10. |
| `stacked_dimensionality.R` | row counts by control group × window → `.tex`. |
| `stacked_loo.R` | leave-one-cohort-out on the shipped stack (plain + m4), corrected FE (§9). |
| `harmonize_atoms.R` | put CS/SA/stacked atoms on a common event-time grid → `atoms_harmonized.csv` (§19). |
| `stacked_pretrends.R` | per-sub-experiment pre-trend joint Wald tests, both stacks, corrected FE (§19). |
| `stacked_loo_both.R` | leave-one-cohort-out on both stacks (shipped m4, reconstructed m3), corrected FE (§19). |
| `fwl_decomp.R` | FWL variance-share decomposition of the R2 spec; identity checks unweighted/weighted, additive-FE contrast (§20). |
| `stack_composition.R` | per-stack treated (drop/adopt/absorbed) × control (never/always/not-yet) counts + identification source → `.tex` (§18). |
| `stacked_nevertreated.R` | stacked regression on never-treated controls only (41) vs full pool (741), corrected FE (§6). |
| `decomp_table.R` | LaTeX table of per-cohort β_s and weights w_s (unweighted + ebal-weighted) from the stacked decomposition → `.tex` (§10). |
| `spec_grid.R` | 2×2×2 FE×covariate×weight grid reproducing the paper's stacked cols 3–4 → `.tex` (§8). |
| `stacked_full15.R` | full 15-cohort stack (controls on all cohorts, no window) run "as described"; spec + standard regression tables → `.tex` (§8). |
| `data_description.R` | descriptive panel summary (years, cities, cohorts, N) → `.tex`. |
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
`cengiz_window_sensitivity.csv`, `stacked_dimensionality.{csv,tex}`,
`stacked_loo.{csv,tex}`, `atoms_harmonized.csv`, `stacked_pretrends.{csv,tex}`,
`stacked_loo_both.{csv,tex}`, `fwl_decomp_{unweighted,weighted,summary}.csv`,
`stack_composition.{csv,tex}`, `data_description.{csv,tex}`,
`stacked_nevertreated.csv`, `decomp_table.tex`, `spec_grid.{csv,tex}`,
`stacked_full15.{csv,tex}`, `stacked_full15_regtable.tex`.

### `figures/`
Nine figures from the plotting scripts above, each `.png` + `.pdf`:
`weight_vs_beta_decomposition`, `weight_vs_beta_smallmultiples`,
`cohort_level_estimates`, `cohort_level_variants`, `calendar_and_event_time`,
`cs_variants_comparison`, `cengiz_vs_pp_stacked`, `cengiz_control_groups`,
`cengiz_window_sensitivity`.

All nine use **`theme_minimal()` with no plot title, subtitle, or caption** (axis
labels, legends, reference lines and in-panel annotations retained);
`plot_calendar_event.R` composes its two panels with `cowplot::plot_grid` (the
env has no `patchwork`). There is no separate `*_minimal` variant — the canonical
figures are themselves minimal.

Not part of this pipeline: `figures/replication.png` and `figures/simulation/*.png`
are diagnostics produced by the separately-merged `stage_*.R` scripts (Stages
V/L/D/N from `main`), not by the fatal-stacked decomposition scripts inventoried
here.

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

## 17. Cohort 2002: present in the shipped stack, but clean-control-free (finding)

Reviewer Step 1 asked whether a 2002 sub-experiment exists in `stacked_fatal.csv`
at all, and whether the within-treated 2002 contrast we documented lives in the
shipped stack or only in the reconstructed 15-cohort stack. Verified against the
shipped data:

- **Cohort 2002 IS in the shipped `stacked_fatal.csv`** — 170 rows, but **zero
  clean controls**. Composition:

  | rows | agencies | direction |
  |---|---|---|
  | treated | 4 — elgin IL, covington KY, canton OH, green bay WI | `no.req` 0→1 (dropped) |
  | treated | 1 — revere MA | `no.req` 1→0 (adopted) |
  | "No Change" controls | **0** | — |

  This is the §4 boundary rule at work: 2002 needs `2002 − 4 = 1998 ≥ 2000` to
  attach controls, which fails, so 2002 is treated-only in the shipped weighted
  stack (as are 2000, 2001, 2003).

- **The within-treated 2002 contrast lives in the SHIPPED stack, not just the
  reconstruction.** Under the corrected agency×stack FE:

  ```
  2002 atom, full shipped stack (4 droppers + 1 adopter): +0.00658
  2002 atom, adopter (revere) removed (4 droppers only):  NA — UNIDENTIFIED (collinear)
  ```

  With no clean controls, the 2002 stacked coefficient is identified *only* by the
  opposing-direction contrast between the four droppers and the single reversing
  adopter (revere MA, whose `no.req` moves 1→0 in the same window). Remove revere
  and `no.req` becomes collinear with the agency×stack + year×cohort FE (all
  survivors move the same way at the same time) → the coefficient is not
  identified. **The entire shipped-stack 2002 estimate hinges on one reversing
  unit.**

- **CS's 2002 is a different object — do not conflate them.** CS identifies 2002
  off the **41 never-treated clean controls** on the full panel: mean ATT(2002,·) =
  **+0.222**, high CS weight, and it is the atom whose removal **flips CS's overall
  sign** (see the CS leave-one-out). So the sharpest statement of the cross-estimator
  contrast is: *2002 is present and clean-control-identified in CS (high-weight,
  sign-flipping), but the shipped stacked 2002 has no clean controls at all and
  survives only through the fragile within-treated dropper-vs-adopter contrast.*

  | | shipped **stacked** 2002 | **CS** 2002 |
  |---|---|---|
  | identified off | 4 droppers vs 1 adopter (within-treated) | 41 never-treated clean controls |
  | estimate | +0.007 | +0.222 |
  | depends on | revere (collinear without it) | clean controls |
  | pooled-sign role | negligible | sign-flipping |

- **The reconstructed 15-cohort *unweighted* stack is the contrast case:** there
  2002 *would* attach the 741 "No Change" controls (that build applies controls to
  all cohorts, not only `g−4 ≥ 2000`; §13–14, and Table A.5 col. 3 = 233,520 rows),
  giving 2002 a clean-control identification the shipped weighted stack denies it.
  This is why the shipped weighted stack (11 control-bearing cohorts, 171,444 N) and
  the reconstructed unweighted stack (15 cohorts, 233,520 N) treat 2002 differently.

Verification was done with scratch scripts against `DIDREP_DATA` (not committed,
per the "work from a copy" instruction); the facts above are the record.

---

## 18. Directional coding: adopters coded inside "drop" stacks (finding)

Reviewer Step 2 asked whether a reversal is being coded as an adoption inside a
stack (a bug) or whether the coding handles it (a sentence). **Both are true — the
sign is handled correctly, but the structure it exposes is load-bearing.**

**Fact.** Stacks are keyed on `year.changed` (the *calendar year* of the policy
change), not on its direction. So a unit lands in cohort `g` whether it dropped or
adopted a requirement in year `g`. Of the **37 treated stack-members**, **6 are
adopters** (`change.type == "Adopted"`, `no.req` moving **1→0**, the opposite of the
canonical "requirement dropped" 0→1 event):

| agency | cohort | within-stack `no.req` |
|---|---|---|
| revere massachusetts | 2002 | 1→0 (adopt) |
| portsmouth new hampshire | 2003 | 1→0 (adopt) |
| memphis tennessee | 2004 | 1→0 (adopt) |
| fall river massachusetts | 2012 | 1→0 (adopt) |
| ramapo new york | 2018 | 1→0 (adopt) |
| springfield massachusetts | 2018 | 1→0 (adopt) |

Plus **2 "absorbed"** members whose reversal falls exactly on a segment boundary so
they carry no within-stack variation (memphis 2009, portsmouth 2016 — their `no.req`
is constant within the second segment). The other 29 treated members are droppers
(0→1).

**The coding handles the sign (the sentence).** The regressor is **`no.req` — the
residency-requirement STATE** (1 = no requirement), *not* a directional treated×post
dummy. An adopter (1→0) and a dropper (0→1) therefore identify the *same* object —
the effect of the no-requirement state — symmetrically; no reversal is mis-signed.
Had the spec used a directional "post" dummy (treated-post = 1 regardless of
direction), the 6 adopters would be mis-coded: their requirement-ON post-period would
carry the same 1 as a dropper's requirement-OFF post-period, conflating opposite
states. It does not — `treat` is only a constant membership flag (absorbed by the
agency×stack FE); all identification runs through `no.req`. So the pooled sign and
magnitude are not mechanically distorted by the adopters.

**But the structure is load-bearing (the finding).** Under the corrected
agency×stack FE, identification of a control-free early cohort depends entirely on
whether it contains a *wrong-direction* member. Among the four cohorts with **no
clean controls** (`g − 4 < 2000`: 2000, 2001, 2002, 2003):

| cohort | composition | identified? | atom |
|---|---|---|---|
| 2000 | 5 droppers, no adopter | **dropped (NaN)** | — |
| 2001 | 1 dropper | **dropped (NaN)** | — |
| **2002** | **4 droppers + 1 adopter (revere)** | **YES** | +0.007 |
| 2003 | 1 adopter only (portsmouth) | **dropped (NaN)** | — |

2000/2001 collapse because every treated unit switches `no.req` in the same
direction in the same year → collinear with the year×cohort FE. 2003 collapses
because a single unit has no within-stack contrast. **2002 survives *only* because
the adopter revere moves `no.req` the opposite way to its four droppers**, breaking
the collinearity — exactly the Step-1 result (remove revere → 2002 becomes NaN).
This is why 2002 is the one control-free early cohort that enters the 11-atom
decomposition at all (§9), and why its estimate rests on a single reverse-direction
unit.

Adopters with clean controls (2012 = fall river only; 2018 = 1 dropper + 2 adopters;
2004 = 1 dropper + 1 adopter) are identified against the 741 controls, so their
survival does not hinge on the adopter — but their treated variation is still a mix
of drop and adopt transitions of `no.req`.

**Bottom line.** No reversal is mechanically mis-coded (the state-variable spec is
correct), but the "Requirement Dropped" label is imprecise for 6 of 37 treated
members, and for the clean-control-free cohorts the adopters are not a nuisance —
they are the *sole* source of identification (2002) or the reason a cohort would
otherwise vanish. The Cengiz windowed designs sidestep this entirely by dropping the
two multiply-reversing units (§5, §11); the shipped full-panel stack keeps them.

**Full per-stack composition** — treated members split by switch direction (drop /
adopt / absorbed) and controls by treatment history (never- / always- / not-yet-
treated), with the identification source of each stack — is in
`R/stack_composition.R` → `output/stack_composition.{tex,csv}`. It makes the 2002
case explicit at a glance: 2002 has 0 controls and 4 droppers + 1 adopter, so it is
identified purely by the within-treated drop-vs-adopt contrast, whereas its
control-free siblings 2000 (5 droppers, no adopter), 2001 and 2003 (single unit)
have no contrast and drop out. The control rows also show the §6/§12 pool
asymmetry: every control-bearing stack carries the same **41 never-treated + 700
always-treated + 0 not-yet-treated** controls.

---

## 19. Pre-trends and leave-one-out (findings)

Reviewer Step 3: harmonize the atoms onto a common event-time grid, run
per-sub-experiment pre-trend tests, and leave-one-cohort-out on both stacks.

### 19.1 Harmonized atoms (`harmonize_atoms.R` → `atoms_harmonized.csv`)

One row per `(estimator, cohort, event_time, estimate, se, is_pre, identified)`:
- **CS**: `ATT(g,t)`, `event_time = t − g` (t<g rows are the placebo pre-periods).
- **SA**: `CATT(g,e)`, `event_time = e`.
- **stacked**: per-cohort event study under the corrected interactive FE,
  reference `e = −1`, restricted to the **±4** balancing window. CS/SA rows span
  their full event-time range; the stacked rows are ±4 (the design's own window),
  so restrict CS/SA to ±4 before a like-for-like comparison. The CS/SA SEs here are
  **marginal** — a joint pre-trend test needs the full covariance and is computed
  by refitting (§19.2), not read off this file.

### 19.2 Per-sub-experiment pre-trend tests (`stacked_pretrends.R`)

For each cohort's sub-experiment, an event study `any.fatalities ~ i(etime, treat,
ref=−1) | agency.id + year` on the ±4 window (within one stack the interacted FE
reduce to agency + year), then a **joint Wald test of the pre-treatment leads**
`{e<0, e≠−1}=0`. Primary inference clusters on `agency.id`; a **hetero-robust**
p-value is reported alongside because these cohorts have **few treated clusters
(often 1–3)**, where cluster-robust inference degenerates and inflates the Wald
statistic.

A cohort is **testable only with clean controls AND an observed treated lead:**
- **Control-free early cohorts (2000–2003) in the shipped stack are UNTESTABLE** —
  no clean controls ⇒ the `treat×etime` leads are collinear. So the shipped stack's
  early cohorts have not only a fragile point estimate (§17–18) but an **untestable
  parallel-trends assumption**. Attaching controls (reconstructed-15) makes 2002 and
  2003 testable (2002 borderline, cluster/hetero p ≈ 0.06 on its single lead; 2003
  p ≈ 0.30).
- **2000, 2001, 2016 are untestable even with controls** — no observed treated
  pre-period (2000's treatment = the first outcome year; 2001 observes only `e=−1`;
  2016's treated unit, portsmouth, is only observed post).

**Results (control-bearing cohorts, ±4):** pre-trends are consistent with parallel
trends for **2004, 2008 (marginal, p≈0.058), 2009, 2012, 2013, 2014, 2020**. Two
late cohorts reject:
- **2019** (1 treated unit): lead coefficients ≈ **+1.0 on a binary outcome** — a
  single always-fatal agency; rejects under **both** cluster and hetero SEs. Genuine.
- **2017** (3 treated units): cluster p≈0 is a few-cluster artifact (lead SE 0.02),
  but hetero p≈0.007 with leads ≈ +0.7–1.0 — a real violation.
- **2018** (3 treated units): cluster p=0.045 but **hetero p=0.17** — the rejection
  is a few-cluster artifact; not a robust violation.

**Takeaway:** the well-populated cohorts pass; the failures are concentrated in
single-/few-treated-unit late cohorts dominated by one large agency, and the extreme
Wald magnitudes overstate them (few-cluster inference). The substantive caveat is
that the shipped stack cannot test pre-trends at all for its control-free early
cohorts — exactly the cohorts whose estimates are most fragile.

### 19.3 Leave-one-cohort-out, both stacks (`stacked_loo_both.R`)

Corrected FE; drop each cohort's whole sub-experiment and refit. `plain` (no
controls/weights) is the apples-to-apples target; the covariate column is **m4
(weighted, Table 3) on the shipped stack** and **m3 (covariates, unweighted) on the
reconstructed-15 stack** (ebal weights not reconstructible).

- **No single cohort flips the stacked sign** on either stack — plain stays in
  **[−0.115, −0.077]** (shipped) / **[−0.107, −0.077]** (reconstructed). The stacked
  negative estimate is **robust** to leave-one-out. *(Contrast: CS's positive overall
  hinges on the single clean-control 2002 atom — dropping 2002 flips CS's sign. The
  stacked estimate has no such single point of failure.)*
- **Dropping 2000/2001/2003/2016 moves the shipped plain coef by exactly 0** —
  direct confirmation they are unidentified and contribute nothing (§17–18).
- **2018 is the most influential single cohort:** shipped plain Δ=+0.017, and under
  **m4 dropping 2018 nearly halves the effect (−0.089 → −0.044, Δ=+0.045)** — the
  same cohort whose pre-trend rejection was a few-cluster artifact (a thin,
  influential late cohort). **2013** pulls the other way (Δ=−0.016/−0.014). **2009**
  (the high-weight cohort) shifts plain by +0.011.
- The reconstructed-15 baseline (plain −0.093, m3 −0.092) is close to the shipped
  −0.099, and its LOO pattern matches — attaching controls to the early cohorts
  barely moves the pooled estimate but makes 2002/2003 identified and testable.

---

## 20. FWL variance-share decomposition of the R2 spec (`fwl_decomp.R`)

Explicit-regression verification of the TWFE identity `pooled_beta = Σ_s w_s·β_s`
for the R2 stacked spec (**unit-by-stack = agency×cohort** + **time-by-stack =
year.cohort** FE). `v_hat` = residual of `no.req` partialled on both stack-
interacted FE sets by explicit regression (the stacks are unbalanced, so the
two-way demeaning closed form does not apply); `V_s = Σ_i v_hat_i²` within stack;
`w_s = V_s/ΣV_s`; `β_s` = per-stack regression of the outcome on `no.req` with
unit + time FE. Estimation sample = observed-outcome rows (**171,906**), the sample
the pooled beta is fit on. Reports faithfully — discrepancies flagged, not fixed.

- **Step 2 (within-stack mean of `v_hat`):** max abs = **3.6e-17** unweighted,
  **1.9e-16** weighted (both effectively zero).
- **Step 3 — `V_s = 0` stacks (no within-stack residual variation ⇒ β_s
  undefined, weight 0):** **2000, 2001, 2003, 2016** (unweighted). These are the
  cohorts where all treated units switch `no.req` in the same direction at the same
  time (collinear with the time-by-stack FE); 2002 survives (V_s = 5.79) despite
  having no controls because its adopter moves `no.req` the opposite way (§17–18).
  The high-weight stack is **2009 (w = 0.369)**; 2013 (0.153) and 2017 (0.100) next.
- **Steps 5–6 — identity holds to machine precision (unweighted):**
  ```
  pooled beta  = -9.891537325820261e-02
  Σ_s w_s·β_s  = -9.891537325820207e-02
  difference   =  5.41e-16
  ```
  Exact because both FE dimensions are stack-interacted ⇒ the pooled normal
  equations are block-diagonal across stacks up to the shared `no.req` coefficient.
- **Step 7 — weighted (entropy-balancing) identity also holds to machine
  precision:**
  ```
  pooled beta  = -1.083934962396384e-01
  Σ_s w_s·β_s  = -1.083934962396388e-01
  difference   = -4.16e-16
  ```
  Weighting reshuffles the pieces (2009's weight 0.369→0.170; 2018's 0.080→0.149;
  treated shares rise) but the decomposition stays exact.
  **Discrepancy (reported, not fixed):** weighted stack **2016** has V_s = 7.55e-44
  (numerically nonzero vs the exact 0 it takes unweighted) and a garbage
  β_s = −193.98 (near-collinear `no.req` over a ~0 residual) instead of the `NA` the
  unweighted per-stack fit returns; its weight is ~1.6e-45 so the contribution
  (~−3e-43) leaves the identity intact. A floating-point collinearity-handling
  artifact that differs between the weighted and unweighted per-stack fits.
- **Step 8 — additive (non-interacted) unit FE + time-by-stack FE: identity FAILS
  by +5.48e-03 (≈5.5%), as expected:**
  ```
  pooled beta (agency.id additive) = -1.043972170916846e-01
  Σ_s w_s·β_s (from R2 decomp)      = -9.891537325820207e-02
  difference                        =  5.481843833482533e-03
  ```
  With additive unit FE the design is no longer block-diagonal: the **743 of 776**
  agencies that appear in >1 stack (all 741 never-changing controls + the 2
  reversibles) tie the stacks through a shared intercept, so the pooled coefficient
  is not the variance-weighted average of the within-stack β_s. Same FE-leakage gap
  as §9 (−0.10440 additive vs −0.09892 interacted).
- **Aside — sample choice:** residualizing `v_hat` over **all** 278,324 rows
  (including pre-2000 missing-outcome rows) then summing over observed-outcome rows
  gives Σ w_s'·β_s = −9.051e-02, off the pooled by **8.40e-03 (≈8.5%)** — the §10
  gap; the estimation sample is the correct one.

**What each step ran on (15 vs 11 cohorts — clarification).** All of steps 1–8 run
on the **one** file `stacked_fatal.csv`, estimation sample = **171,906 rows, 15
cohorts** (2000–2020); step 7 uses that file's **native `weights` column** — no
second dataset, no foreign weights. The "11" that appears elsewhere is *not* the
file's cohort count; the shipped stack is 15 cohorts. Two different 11-subsets exist
and they are **not** the same set:
- **11 control-bearing cohorts** (2004, 2008, 2009, 2012, 2013, 2014, 2016, 2017,
  2018, 2019, 2020) — the ones with clean controls, hence the ones the entropy
  weights actually balance. The **4 treated-only cohorts** (2000, 2001, 2002, 2003)
  carry **weight = 1 on every row** (nothing to balance).
- **11 identified cohorts** under the interacted FE (2002, 2004, 2008, 2009, 2012,
  2013, 2014, 2017, 2018, 2019, 2020) — those with V_s > 0. This set **includes 2002**
  (treated-only but identified via its adopter) and **excludes 2016** (absorbed
  reversal, V_s = 0). §9.

Consequently the weighted pooled coefficient depends on which "11" is meant:

| spec | cohorts | N | pooled β |
|---|---|---|---|
| unweighted, all 15 (**R2**) | 15 | 171,906 | −0.098915 |
| weighted, all 15 (**step 7**) | 15 | 171,906 | −0.108393 |
| weighted, **11 identified** (incl. 2002, excl. 2016) | 11 | 156,198 | −0.108393 |
| weighted, **11 control-bearing** (incl. 2016, excl. 2002) | 11 | 171,659 | −0.124043 |

- The step-7 figure (−0.1084) is **valid for the file as shipped** and equals the
  **11-identified** pooled to 15 digits — the 4 unidentified cohorts (2000, 2001,
  2003, 2016; all V_s = 0) contribute ~1e-16.
- Restricting instead to the **11 cohorts the weights were built for** gives
  **−0.1240 — matching neither R2 nor step 7.** The gap is cohort **2002**: treated-
  only (weight = 1) but identified via its reverse-direction adopter (revere,
  β₂₀₀₂ = +0.007); keeping it pulls the weighted pooled toward zero.
- So the weights were **not** applied to a stack they weren't constructed for (native
  column of the same file), **but** the shipped weighted estimate is a *mix* — 11
  entropy-balanced cohorts plus 2002 riding in unweighted — and thus inherits the
  same adopter-dependent 2002 fragility as §17–18. Fully-balanced (drop the treated-
  only cohorts) the weighted pooled is −0.1240.

Outputs: `output/fwl_decomp_{unweighted,weighted,summary}.csv`.

---

## 21. Core code (embedded verbatim)

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
