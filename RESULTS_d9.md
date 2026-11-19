# D9 — unchecked weight under every aggregation choice

`R/d9_aggregation_U.R` → `tables/table_aggregation_U.tex`, `output/d9_aggregation_U.csv`,
`output/d9_eventstudy_by_horizon.csv`. Droppers sample (§4.3–4.5). No new estimation: per-cohort
geometry (cohort, $n_g$, pre$_g$, post$_g$) and MDEs are the round-6 values in `d6_round6_mde.csv`;
reference rows are the round-6 single-weighting values in `d6_unchecked_weight.csv`. Aggregation
weights are built from panel geometry to match `did::aggte` defaults.

> **Toolchain note.** The scratchpad R environment was wiped again this session (no R, no `did`
> package), so the outputs were produced by a byte-identical Python port of `R/d9_aggregation_U.R`
> that reads the same CSVs. `R/d9_aggregation_U.R` is the canonical generator.

## Headline

| estimator | aggregation | $U(\infty)$ | $U(\theta^*{=}0.05)$ | $U(\theta^*{=}0.10)$ |
|---|---|---|---|---|
| CS-NT / CS-NYT / SA / imputation | Simple | **0.243** | 1.000 | 1.000 |
| CS-NT / CS-NYT / SA / imputation | Cohort (group) | **0.292** | 1.000 | 1.000 |
| CS-NT / CS-NYT / SA / imputation | Event study, overall | **0.292** | 1.000 | 1.000 |
| CS-NT / CS-NYT / SA / imputation | Calendar, overall | **0.292** | 1.000 | 1.000 |
| Stacked (reference) | full panel | 0.232 | 1.000 | 1.000 |
| Cengiz ±4 (reference) | ±4 window | 0.250 | 1.000 | 1.000 |

**All four estimators give identical columns.** CS-NT, CS-NYT, SA and imputation share the same
droppers geometry (same cohorts, $n_g$, post$_g$), and the small differences in their per-cohort
MDEs (all testable MDEs ≈ 0.039–0.96) are immaterial at these $\theta^*$: every testable cohort is
already unchecked under every overall aggregation, so only the *weights* would move the numbers, and
the weights are the same across these four estimators. SA ≡ CS-NT numerically (never-treated
controls ⇒ same geometry **and** same MDEs).

## Which aggregation minimizes $U$, and by how much

- **On the break-even margin ($\theta^*=0.05,0.10$): all four aggregations tie at $U=1.000$.** No
  overall aggregation achieves *any* checkable weight. Once a cohort's effect is averaged over its
  whole post-window (or over all event-study horizons), its accumulation factor $a_g$ is large
  enough that the break-even slope $\bar\delta_g=\theta^*/a_g$ falls below the ≈0.04 MDE floor, so
  every testable cohort is unchecked. The diagnostic does **not** separate the overall aggregations
  here.
- **On untestability ($U(\infty)$): simple minimizes, at 0.243 vs 0.292** for the other three — i.e.
  cohort/event-study/calendar each carry **+0.048** more untestable weight than simple. So the only
  aggregation that lowers the diagnostic is **simple**, and it does so on the $U(\infty)$ margin, not
  the break-even margin.

This is the paper's point in miniature: *the aggregation choice moves the share of weight on
cohorts you cannot check, and here the movement is against the aggregations usually preferred for
interpretability (group/event-study).*

## Does cohort aggregation lower $U(\infty)$ as expected? — No, it raises it (flag)

The brief anticipated that stripping the post$_g$ factor (simple → cohort) would **lower**
$U(\infty)$. It does the **opposite**: $U(\infty)$ rises from 0.243 (simple) to **0.292 = 7/24**
(cohort), +0.048.

Reason, specific to this sample: in droppers, every untestable cohort is a **single-unit** cohort
($n_g=1$: 2001, 2004, 2008, 2014, 2018, 2019, 2020 — seven of them), while the weight-dominant
cohort 2009 is testable ($n_g=7$). Simple weighting ($\propto n_g\,\text{post}_g$) concentrates
weight on the well-powered 2009 ($84/259=0.324$) and gives the seven singletons only
$\sum\text{post}_g/259 = 63/259 = 0.243$. Cohort weighting ($\propto n_g$) spreads weight by
treated-unit share: 2009 drops to $7/24=0.292$ and the seven singletons collectively **rise** to
$7/24=0.292$. Event-study-overall and calendar-overall land at essentially the group value
(0.29162), since both anchor cohort weights on $P(G=g)$-type shares.

So cohort/event-study/calendar aggregation, by treating cohorts evenhandedly, hands *more* weight to
the numerous unpowered singleton cohorts. Simple aggregation's post-window reward happens to favor
the one well-powered cohort. **Whether cohort aggregation helps or hurts $U(\infty)$ is
sample-specific, not a general ordering** — worth stating plainly in the text.

## How quickly does $U_e$ approach 1 as $e$ grows (event study)

`d9_eventstudy_by_horizon.csv`. Under the $g-1$ baseline, $a_g = e+1$ at horizon $e$, so the
break-even slope $\bar\delta_g=\theta^*/(e+1)$ shrinks with the horizon. With the MDE floor ≈ 0.04:

| horizon $e$ | $U_e(\theta^*{=}0.05)$ | $U_e(\theta^*{=}0.10)$ |
|---|---|---|
| 0 | 0.458 | 0.458 |
| 1 | **1.000** | 0.435 |
| 2 | 1.000 | **1.000** |
| ≥2 | 1.000 | 1.000 |

**$U_e$ reaches 1 almost immediately** — by $e=1$ at $\theta^*=0.05$ and by $e=2$ at
$\theta^*=0.10$. Only the **contemporaneous** effect ATT($e{=}0$) is checkable at $\theta^*=0.05$
($\bar\delta=0.05 \ge$ MDE 0.041–0.043 for cohorts 2009/2013/2017; checked weight $0.542$), and only
$e\le 1$ is (partly) checkable at $\theta^*=0.10$. One horizon of accumulation ($a_g:1\to2$) halves
the tolerance and pushes the testable cohorts below detectability. This is the sharpest illustration
of the diagnostic: *a pre-trend test that "passes" for the on-impact effect is already uninformative
about the dynamic effects two years out.* (`U_e(\infty)` wanders between 0.18 and 0.43 as the
estimable set shrinks, then jumps to 1.0 at $e=19$ where only the untestable 2001 cohort remains.)

Marginality worth noting: cohort 2017 sits right on the break-even edge — MDE $=0.0413$ vs
$\bar\delta=\theta^*/a_g=0.10/2.5=0.0400$ under cohort aggregation — so it is unchecked by only
0.0013. At its own first horizon ($e=0$, $a_g=1$) it is comfortably checked. The overall
aggregations are "all unchecked" partly because this near-miss falls the wrong way.

## Validation

1. **Simple reproduces round 6 — PASS.** CS-NT and imputation simple $U(\infty)=0.243243$ (target
   0.243243 from `d6_unchecked_weight.csv`); $U$ at both $\theta^*=1.000$, matching the committed
   droppers rows exactly.
2. **Weights match `did` — hand-verified (did unavailable this session).** The round-6 **simple**
   weights I anchor on were themselves produced by `did::att_gt + aggte(type="simple")` (gradient
   method, `d6_weight_checkability.R:96–105`), and this run reproduces them to the digit. Checks run
   by hand against the formulas: (a) **group** weights $W_g=n_g/\sum n_g=P(G=g)$ (e.g. 2009:
   $7/24=0.2917$), matching `aggte(type="group")`'s cohort weighting; (b) **event-study** cell
   weights $w_{g,g+e}=\frac1E P(G=g\mid G+e\le T)$ sum to 1 at **every** horizon $e$ (0 violations
   across all estimators), matching `aggte(type="dynamic")`'s equal-weight-over-horizons overall;
   (c) **calendar** cell weights sum to 1 over post-periods. The package equality check itself could
   not be executed without R/`did`.
3. **Weights sum to 1 — PASS.** max$|{\sum_g W_g}-1| = 2.2\times10^{-16}$ over every estimator ×
   aggregation.
4. **Bounds — PASS.** $U(\theta^*)\ge U(\infty)$ and $U$ non-increasing in $\theta^*$ for every row.

All MDEs are the **optimistic iid/WCR** round-6 values, so every $U$ here is a **lower bound** on the
unchecked weight.

## Deliverables

`R/d9_aggregation_U.R`, `tables/table_aggregation_U.tex`, `output/d9_aggregation_U.csv`,
`output/d9_eventstudy_by_horizon.csv`, and this file.
