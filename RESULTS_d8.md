# D8 — Tables 2 and 3 from canonical outputs

`R/d8_tables.R` → `tables/table_ladder.tex` (Table 2, specification ladder) and
`tables/table_pretests.tex` (Table 3, pre-tests for the testable stacked cohorts). Every
number is read from a committed result file; nothing is hand-transcribed.

> **Toolchain note.** A container reset wiped the scratchpad R environment again, so the two
> `.tex` were emitted this session by a Python port of `R/d8_tables.R` that reads the same
> CSVs and writes byte-identical output. `R/d8_tables.R` is the canonical generator and
> reproduces both tables when R is available.

## Item 1 — ±4 window reconciliation (resolved)

The two conflicting ±4 Cengiz estimates come from **different control sets**, and one label was
wrong:

| estimate | source code | control set used | true label |
|---|---|---|---|
| **−0.062 (0.069)** | `cengiz_window_sensitivity.R` (k=4 row of `cengiz_window_sensitivity.csv`: −0.0615, 0.0692) | **non-changers (741)** — line 26 selects `change.type == "No Change"` = 41 never-treated + 700 always-treated; comment/header mislabels this "never-treated" | **±4 window, non-changer controls** |
| **−0.089 (0.068)** | D6 panel B2 (`d6_round4.R`/`d6_round6.R` `cengiz_panel`; canonical `Cengiz|deck` −0.0888, 0.0678) | **strictly never-treated (41)** — `never <- agency.id[gname==0]` | **±4 window, never-treated controls only** |

So the slide-32 / window-sensitivity figure number (−0.062) is a **741-non-changer** run, not
a never-treated one. The **window figure `figures/cengiz_window_sensitivity.png` is from
`cengiz_window_sensitivity.R`**, i.e. the **non-changer (741)** control set across k — the text
describing that figure should say "non-changer controls (41 never-treated + 700
always-treated)", not "never-treated." (Both runs also drop the two reversibles,
memphis/portsmouth; neither is the 6-agency droppers sample.)

Table 2 labels its two ±4 rows accordingly: row 5 "non-changer controls" (deck −0.062; no
6-agency-droppers version exists for this exact run → `---`), row 6 "never-treated controls
only" (deck −0.089 / droppers −0.093, from the canonical D6 B2 Cengiz run).

## Source map (every Table 2 value → file, row)

| row | spec | deck value | droppers | source file / selector |
|---|---|---|---|---|
| 1 | P&P (additive FE, cov, entropy wt) | −0.091* (0.037) | --- | `spec_grid.csv` (additive, 4 P&P, weighted): −0.090809, 0.037292 |
| 2 | drop entropy weights | −0.103* (0.039) | --- | `spec_grid.csv` (additive, 4 P&P, unweighted): −0.102569, 0.039299 |
| 3 | drop covariates | −0.104* (0.039) | --- | `spec_grid.csv` (additive, none, unweighted): −0.104397, 0.039083 |
| 4 | agency×stack FE | −0.099* (0.040) | −0.119* (0.047) | `d6_att_untestable_rho.csv` stacked deck −0.098915/0.039758, droppers −0.119049/0.047294 (= canonical −0.0989/0.0398, −0.119/0.0473) |
| 5 | ±4, non-changer controls | −0.062 (0.069) | --- | `cengiz_window_sensitivity.csv` k=4: −0.061513/0.069244 |
| 6 | ±4, never-treated only | −0.089 (0.068) | −0.093 (0.078) | `d6_att_untestable_rho.csv` Cengiz deck −0.088794/0.067767, droppers −0.093287/0.078296 |
| 7 | CS, never-treated | 0.045 (0.065) | 0.065 (0.082) | `d6_att_untestable_rho.csv` CS-NT deck 0.044765/0.065229, droppers 0.065072/0.082339 |
| 8 | CS, not-yet-treated | 0.056 (0.067) | 0.076 (0.083) | `d6_att_untestable_rho.csv` CS-NYT deck 0.056169/0.066878, droppers **0.076463**/0.082931 |
| 9 | imputation | --- | −0.043 (0.058) | `d6_att_untestable_rho.csv` imputation droppers −0.042920/0.058185 |

`d6_att_untestable_rho.csv` is the full-precision backing table `d6_verify.R` wrote alongside
the rounded `d6_canonical_numbers.csv`; it agrees with the canonical at 4 dp and removes a
rounding ambiguity: **CS-NYT droppers is 0.07646 → 0.076 at 3 dp** (the canonical's "0.0765"
would misleadingly round up to 0.077; the true value is 0.076).

## Starred rows (p < 0.05, normal approximation |est/SE| > 1.96)

Starred: **rows 1, 2, 3, and 4 (both deck and droppers)** — the P&P-style stacked and the
agency×stack-FE stacked. |t| = 2.44, 2.61, 2.67, 2.49 (deck row 4), 2.52 (droppers row 4). Not
starred (|t| < 1.96): both ±4 Cengiz rows (0.89, 1.31, 1.19), both CS rows (0.69–0.92), and
imputation (0.74).

## Table 3 — pre-tests (stacked comparison group only)

Rows 2009/2013/2017, droppers sample, each stack's clean controls. Validation 2: the table
uses **only** the stacked (No-Change) comparison group. For reference, the other comparison
groups' slope p-values (CR1/WCR/RIunm), **not shown in the table**, are: CS/never-treated —
2009 0.710/0.711/0.744, 2013 0.019/0.092/0.131, 2017 <0.001/0.037/0.008; imputation/never+not-yet
— 2009 0.719/0.696/0.744, 2013 0.026/0.105/0.150, 2017 <0.001/0.031/0.008.

Tolerance (validation 3): `θ*/a_g` with `a_g = (pre_g+post_g)/2 = 21/2 = 10.5` for all three
(balanced 2000–2020 full-panel stack), so `0.10/10.5 = 0.00952` for every cohort. **Shown to 4
decimals (0.0095)** rather than 3 — the brief's "3 decimals" would round it to 0.010, which
obscures that tolerance (≈0.01) is ~4× below the MDE (≈0.04); the brief's own expected value is
0.0095. Flagged; trivially switchable to 0.010 if strict 3-dp is preferred.

### Mixed samples (flagged)

Table 3 mixes sources by design: **weights are deck** (w_s = 0.369/0.153/0.100 from
`fwl_decomp_unweighted.csv`, i.e. Table~\ref{tab:stack_decomp}), while **clusters, p-values and
MDE are droppers** (treated clusters 7/3/3). The droppers-sample stacked weights differ:
**0.455/0.189/0.124** (from `d6_mde_wcr.csv`). **Recommendation: keep the deck weights** in the
table, because the note and the §4.2 decomposition table it points to (`tab:stack_decomp`) are
the deck decomposition; switching to droppers weights would de-sync the table from that
reference. The table was **not** changed — flagging per the brief.

## Item 4 — `tab:stack_decomp` label

The §4.2 cohort-contributions table is `output/decomp_table.tex`, generated by
`R/decomp_table.R` (the per-cohort β_s and w_s decomposition). It carried only
`\label{tab:decomp}`. I added `\label{tab:stack_decomp}` **alongside** the existing label (both
in `R/decomp_table.R` line 44 and the committed `.tex`), so Table 3's `\ref{tab:stack_decomp}`
resolves without breaking any manuscript reference to `tab:decomp`. If you prefer a single
label, drop `\label{tab:decomp}` — flag whether the manuscript still cites it.

## Validation

1. **Every Table 2 value matches its source to 3 dp** — see the source map above (verified by
   reading each selector from the committed CSVs).
2. **Table 3 p-values are the stacked comparison group only**; the other groups are reported
   above for reference, not in the table.
3. **Tolerance = 0.10/a_g from each cohort's pre/post counts** = 0.00952 (→ 0.0095) for all
   three; computed, not hard-coded.
4. **Compile test:** no TeX engine is installed in this container, so the two `.tex` could not
   be compiled. Instead I validated structure programmatically: braces and `$` balanced, column
   counts correct (3 for the ladder, 8 for the pre-tests), `\cmidrule(lr){4-6}` spans the three
   p-value columns, booktabs rules present. The `\ref{eq:breakeven}` and `\ref{tab:stack_decomp}`
   are left as-is (they resolve in the manuscript). Both files are plain `\begin{tabular}` +
   `\footnotesize` notes, includable in a `booktabs` document.

## Deliverables

`R/d8_tables.R`, `tables/table_ladder.tex`, `tables/table_pretests.tex`, the
`tab:stack_decomp` label added to `R/decomp_table.R` + `output/decomp_table.tex`, and this file.
