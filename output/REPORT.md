# Control-pool verification & harmonized ladder — summary report

Payson & Parinandi (2024), *Residency Blues* — fatal-encounters outcome
(`any.fatalities`, Table 3). Stages **V → L → D → N**. Numbers and flags only;
no substantive interpretation (per prompt rule 5).

**Reproducibility.** conda-forge R 4.5.3; `did` 2.5.1, `fixest` 0.13.2, `lfe`
3.1.1; `set.seed(0)`. The environment reproduces every documented baseline to
machine precision: CS/SA overall +0.044765, stacked interacted-FE −0.098915,
naive-residualization −0.082281. Scripts: `R/stage_v_control_pool.R`,
`R/stage_l_ladder.R`, `R/stage_dn.R` (reuse `R/extract_atoms.R`). Full log +
`sessionInfo()` in `output/run.log`.

---

## Stage V — Control-pool verification (GATE): **PASSED / hypothesis CONFIRMED**

The stacked design's 741 `change.type=="No Change"` controls split **41 / 700**,
exactly as hypothesized.

| | count |
|---|---:|
| No-Change agencies, always `no.req==0` | **41** |
| No-Change agencies, always `no.req==1` | **700** |
| No-Change agencies, mixed | 0 |

- **V1.** Split is 41/700. Cross-tab against the `callaway_replication.R` gname
  recode agrees **exactly** (always-0 ↔ gname 0; always-1 ↔ gname 1987): **0
  disagreements**. The 46 "mixed" agencies (of all 787) are exactly the changers.
- **V2.** Of 741 control agencies (and 8,151 agency×stack units), **0** have any
  within-unit `no.req` variation; the 700 always-1 controls vary in **0**. They
  contribute only fixed effects and supply **161,700 / 171,906 = 94.1%** of the
  rows fitting `year.cohort`.
- **V3.** 2000–2020 trend of `any.fatalities`: always-1 slope 0.00769 (SE
  0.00063), always-0 slope 0.00535 (SE 0.00245); both upward; equal-slope test
  **p = 0.380** (not rejected); fitted 2000→2020 change +0.154 vs +0.107.
- **V4.** **700 agencies** are simultaneously CS-already-treated (`gname==1987`,
  `callaway_replication.R:17` recode → `att_gt` drops them, lines 52–58) **and**
  stacked clean controls (`treat==0`). Exactly the 700 always-1 non-changers.

---

## Stage L — Harmonized ladder

Support rule (mechanical): supports P&P iff estimate < 0 **and** significant at 5%.

| id | label | pool | n_ctrl | estimate | se | supports P&P |
|----|-------|------|-------:|---------:|---:|:---:|
| R0 | Table 3 m4 (shared FE, cov+ebal wts) | 741 | 741 | −0.09081 | 0.00433 | **yes** |
| R1 | R0 plain, shared FE | 741 | 741 | −0.10440 | 0.03827 | **yes** |
| R2 | R1 + interacted agency×stack FE | 741 | 741 | −0.09892 | 0.03944 | **yes** |
| R3 | Cengiz ±4, never-treated | 741 | 741 | −0.06151 | 0.06924 | no |
| R4 | Cengiz ±4, not-yet-treated | 741+nyt | — | −0.06163 | 0.06921 | no |
| R5 | CS=SA never-treated, no cov, overall | 41 | 41 | **+0.04477** | 0.06503 | no |
| **R2h** | **stacked interacted FE, harmonized** | **41** | **41** | **−0.09335** | **0.03635** | **yes** |
| R3h | Cengiz ±4, harmonized | 41 | 41 | −0.09878 | 0.06549 | no |
| R5d | CS full 741 pool *(DIAGNOSTIC ONLY)* | 741 | 741 | +0.03028 | 0.05826 | no |
| R5b_R2 | CS-simple wts × R2 atoms (741) | 741 | 741 | −0.06114 | — | no (no SE) |
| R5b_R2h | CS-simple wts × R2h atoms (41) | 41 | 41 | −0.06148 | — | no (no SE) |
| R6 | CS ATT(g)/ATT(e) disaggregated (41) | 41 | 41 | +0.04477 | 0.06503 | no |

**Key isolations / flags**
- **R2h − R2 = +0.005566 — sign does NOT flip.** Restricting the stacked pool
  741→41 (nothing else changed) moves the estimate only from −0.09892 to −0.09335,
  within SE; R2h stays **negative and significant** (t = −2.57); all 11 cohorts
  retained, none lost. → The 741-vs-41 control pool is **not** what separates the
  stacked estimate from CS.
- **R5b flag = AGGREGATION-driven (both move).** CS weights on R2 atoms → −0.06114;
  on R2h atoms → −0.06148 (CS target +0.04477). Both move ≈ +0.037 toward R5 and
  land at ≈ the same value regardless of pool; the residual gap to +0.045 is the
  atoms/estimand, not the pool.
- **R5d (diagnostic) = +0.03028** vs R5 +0.04477: forcing the 700 always-1 units
  into CS as never-treated pulls CS down only −0.0145 and it **stays positive**.
- **R6:** 13 cohort ATT(g); **53.8% negative**, **38.5% individually significant**
  at 5%. Strongly heterogeneous (2004 −0.260, 2017 +0.722, 2013 +0.302 sig).

---

## Stage D — Diagnostics

- **D1 effective cohorts** (1/Σw²): R2 5.23, R2h 5.68, R5 5.64 (of 11–13).
- **D2 Wald, H₀ cohort effects equal — all four REJECT:** stacked 741 χ²=1676
  (p≈0), stacked 41 χ²=167.7 (p=8e-31), CS 41 χ²=72.1 (p=1.3e-10), CS 741 χ²=457
  (p=3e-90). Dispersion exceeds sampling variation everywhere.
- **D3 leave-one-cohort-out:** the **only** sign flip is **CS dropping cohort
  2002 → −0.02753** (full +0.04477; cohort-2002 CS weight 0.290). No stacked
  cohort omission flips R2 or R2h. Cohort-2009 stacked weight **0.369** (R2),
  0.343 (R2h) — confirms `assumptions.md` §9 ≈0.37.
- **D4 weight concentration:** R2 top1 0.369 / top3 0.623 (2009); R5 top1 0.290 /
  top3 0.649 (2002).
- **D5 no clean controls:** cohorts 2000–2003 have zero controls (both pools);
  2000/2001/2003 unidentified, 2016 dropped (reversible collinear); **2002 is
  identified off within-treated variation only** (fragile, weight ≈0.06).

---

## Stage N — Reporting-number reconciliation

All counts match `assumptions.md` §14 **exactly**:

| quantity | value | kind |
|---|---:|---|
| file rows | 278,324 | file rows |
| outcome-observed rows | 171,906 | file rows (= 171,171 control + 735 treated) |
| felm N, plain | 171,906 | **estimation N** |
| felm N, covariate (m3/m4) | **171,444** | **estimation N** (462 dropped) |

**Labeling audit of `output/`: 0 collisions** — no file labels a windowed row
count and an estimation N the same way. `stacked_dimensionality` correctly labels
its 33k–48k as "Rows"; estimation Ns (171,444 / 171,906) appear only under
estimation-N labels. The "171,444 as Rows" confusion is in the paper's slides, not
in `output/`; `n_reconciliation.csv` is the single typed source of truth.

---

## Consolidated findings

1. **Control-pool hypothesis confirmed** (V): 741 = 41 + 700; the 700 always-1
   units are both dropped by CS and used as stacked controls.
2. **But the pool is not the driver of the sign disagreement** (L): harmonizing
   the stacked pool to CS's 41 units leaves the stacked estimate negative and
   significant (R2h −0.0933, t −2.57); forcing CS onto 741 leaves CS positive
   (R5d +0.030). The stacked/CS gap survives pool harmonization from both sides.
3. **Aggregation, not the pool, moves the reweighted estimate** (R5b): both pools
   reweight to ≈ −0.061.
4. **Cohort heterogeneity is real in every estimator/pool** (D2 all reject); the
   CS positive sign hinges on a single cohort (2002) that is itself fragile (D3,
   D5).
5. **All reporting numbers reconcile** (N); no mislabeling in `output/`.

## Cells not computed (with reason)
- **R5b SEs / CIs:** R5b is a point-estimate reweighting (CS weights × stacked
  atoms); no variance defined for the recombination → NA; support rule → "no".
- **n_units / n_cohorts** on some tabulated existing rows (R3, R4, R5) not carried
  by their source output files; see the `source` column in `ladder.csv`.

## Output index
`control_pool_verification.{csv,md}`, `v3_annual_means.csv`, `ladder.{csv,md}`,
`atoms_harmonized.csv`, `diagnostics.{csv,md}`, `n_reconciliation.{csv,md}`,
`n_label_audit.csv`, `run.log`.
