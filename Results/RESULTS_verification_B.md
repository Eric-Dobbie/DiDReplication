# Verification B — are the pre-trend rejections real? (Tasks B1–B6)

After the `V_analytical` `/n` fix the CS per-cohort pre-trend tests swung from p≈1
to p≈0 for nearly every cohort. This chain asks whether those rejections are real.
**Answer: no — they are over-rejection artifacts of joint Wald tests that impose far
more restrictions than the few treated clusters can support.** Run order B1→B6;
nothing was written from the pre-trend numbers until the test was shown to misbehave
under the null (B2). Environment: R 4.3.3, `did` 2.1.2, `fixest` 0.13.2.

---

## B1 — Restrictions vs clusters  (`taskB1_restrictions_vs_clusters.R`)
`output/taskB1_restrictions.csv`. **10 of 12 cohorts have restrictions ≥ treated
agencies**, and p-values sort with the ratio (Spearman(ratio, p) = **−0.74**).

| cohort | n_pre | n_treated | ratio | p (CS) |
|---|---|---|---|---|
| 2002 | 1 | 5 | 0.2 | 0.078 |
| 2009 | 8 | 7 | 1.14 | 1.2e-4 |
| 2004 | 3 | 2 | 1.5 | 2.8e-3 |
| 2013 | 12 | 3 | 4.0 | 2.3e-6 |
| 2008 | 7 | 1 | 7.0 | 1e-102 |
| 2012 | 11 | 1 | 11 | **1e-243** |
| 2019 | 18 | 1 | 18 | ~0 |
| 2020 | 19 | 1 | 19 | 2e-9 |

The only cohort with ratio < 1 (2002) is the only one that does not reject. The
astronomical p-values come from single-treated-agency cohorts with 11–19
restrictions — the over-rejection fingerprint. CS and SA show the identical pattern.

## B2 — Does it reject under the null? (decisive)  (`taskB2_permutation_null.R`)
`output/taskB2_null_pvalues.csv`, `taskB2_size_by_cohort.csv`,
`figures/taskB2_null_pvalue_hist.*`. Permutation test, B=500: each draw assigns
placebo adoption years (observed cohort-size structure and spacing) to 28 of the 741
non-changers and uses 41 as controls — no real effect by construction, 41 controls
matching the real design, all cohorts full rank.

**Overall null rejection rate at α=0.05 = 0.713** (nominal 0.05). KS vs Uniform[0,1]:
p = 0. The null p-value histogram spikes at 0 (≈3,900 of 5,500 cells in the first
bin vs 275 expected). Rejection rate rises monotonically with pre-periods:

| n_pre | 2 | 4 | 8 | 9 | 12 | 13 | 14 | 17 | 18 | 19 | 20 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| null reject rate | **0.06** | 0.26 | 0.74 | 0.46 | 0.87 | 0.75 | 0.93 | 0.88 | 0.93 | 0.98 | **0.99** |

**Not the "size holds for larger cohorts" case.** The driver is the restriction count,
not the treated count: cohort 2002 holds size (0.06) because it has only 2
pre-periods despite 5 treated units; cohort 2009 (7 treated) still over-rejects at
**0.46**; the single-unit high-`n_pre` cohorts reject 93–99% of the time with nothing
to find. Every cohort with `n_pre ≥ 4` over-rejects.

**Consequence for the real data:** a cohort's real-data p carries information only
where the test holds size — i.e. only cohort 2002, which does **not** reject
(p=0.078). Cohort 2009's real-data p=0.0001 sits against a 46% false-positive rate,
so it is uninformative; every high-`n_pre` rejection is an artifact.

## B3 — End-to-end validation  (`taskB3B6_validation.R`)
`output/taskB3B6_validation.txt`. The machinery is correct — the over-rejection is a
property of the test, not a bug:
- **`V_analytical` is exactly `crossprod(inffunc)/n`** (max abs diff 0.0e+00), so the
  `/n` correction is a uniform rescale of the whole matrix (off-diagonals included),
  not a diagonal patch.
- **Analytical vs multiplier-bootstrap correlation** agree to 0.045 on the 2009 block
  (−0.673 vs −0.665, −0.871 vs −0.866, …), validating the off-diagonal structure the
  Wald statistic depends on.
- **did's own pooled pre-test `Wpval` is NULL** here: the pooled pre-vcov is rank 31
  of 127 (singular). The package refuses the pooled test; the per-cohort subsets are
  individually full-rank but are the over-rejecting objects B2 exposes.

## B4 — Use the method's own inference  (`taskB4_bootstrap_vcov.R`)
`output/taskB4_bootstrap_vcov.csv`. Per-cohort Wald on the multiplier-bootstrap vcov
vs analytical, both base periods. **They agree on reject/not for 12 of 12 cohorts**,
and the joint test is base-period invariant (a full-rank linear reparameterization
leaves the Wald unchanged). Neither the bootstrap vcov nor the base-period choice
fixes the over-rejection — confirming it is restrictions-vs-clusters.

## B5 — Look at the data  (`taskB5_raw_means.R`)
`figures/taskB5_raw_means.*`, `output/taskB5_raw_means.csv`. Raw outcome means,
treated vs the 41 never-treated controls, no adjustment:
- **Cohort 2002** (5 treated): series track closely, pre-period gap **−0.01**.
- **Cohort 2009** (8 treated): treated run ~**0.14** higher but roughly **parallel**
  in the pre-period — a level offset, differenced out by DiD — with no visible
  divergence, yet the joint test returns p=0.0001. This is the "they track, the test
  says p≈0, so the test is broken" case.

## B6 — Resolve n = 71  (`taskB3B6_validation.R`)
`n = 71` is the legitimate att_gt estimation sample: **41 never-treated controls + 30
treated** (cohorts 2001–2020). The 716 dropped are the already-treated units (700
always-no-requirement coded `gname=1987` + pre-2000 changers), exactly as CS drops
them. Not an unintended restriction; the `/n` scaling is correct.

---

## Verdict
**The tests over-reject under the null → the pre-trend rejections are an artifact.**
The design-validity reading is dropped. This is the first branch of the brief's
decision, not the mixed case: size is controlled only for the single lowest-`n_pre`
cohort (2002), and that cohort does not reject on the real data; every rejection the
real data shows is in a cohort where the test rejects 46–99% of the time under the
null.

**This strengthens the weight-vs-checkability finding (Task A4).** The standard
pre-trend diagnostic is uninformative precisely in the cohorts carrying the most
aggregation weight: 2009 (weight 0.29) has a 46% false-positive rate, and the
high-`n_pre` cohorts are at 90–99%. Pre-trends here cannot be assessed with the joint
Wald test at all; the honest statement is that for these cohort sizes the parallel-
trends assumption is **untestable**, not that it is violated.
