# =============================================================================
# d8_tables.R  --  generate paper Table 2 (specification ladder) and Table 3
# (pre-tests for the testable stacked cohorts) as booktabs .tex, entirely from
# committed result files. No re-estimation; every number is read, not typed.
#
#   tables/table_ladder.tex    (Table 2)
#   tables/table_pretests.tex  (Table 3)
#
# Sources: spec_grid.csv (Table A.4, rows 1-3), d6_att_untestable_rho.csv (the
# full-precision backing table d6_verify wrote alongside the rounded
# d6_canonical_numbers.csv; ATTs/SEs for rows 4,6-9), cengiz_window_sensitivity
# .csv (row 5, k=4; NON-CHANGER controls -- see item 1 in RESULTS_d8.md),
# d6_pretests_all.csv + d6_mde_wcr.csv (Table 3 p-values / MDEs), and
# fwl_decomp_unweighted.csv (deck stacked weights w_s).
# =============================================================================
OUT <- Sys.getenv("DIDREP_OUT", unset="output"); TAB <- Sys.getenv("DIDREP_TAB", unset="tables")
dir.create(TAB, showWarnings=FALSE, recursive=TRUE)
z95 <- qnorm(0.975)                                  # 1.959964
rd  <- function(f) read.csv(file.path(OUT,f), stringsAsFactors=FALSE, check.names=FALSE)

# ---- ATTs (full-precision backing table; agrees with canonical at 4dp) -------
a <- rd("d6_att_untestable_rho.csv"); attp <- list()
for (i in seq_len(nrow(a))) attp[[paste(a$estimator[i], a$sample[i], sep="|")]] <- c(a$ATT[i], a$SE[i])
att <- function(key) attp[[key]]
sg  <- rd("spec_grid.csv")
sgk <- function(fe,cov,w){ r<-sg[sg$fe==fe & sg$covariates==cov & sg$weights==w,]; c(r$coef[1], r$cse[1]) }
win <- rd("cengiz_window_sensitivity.csv"); w4 <- win[win$k==4 & !is.na(win$estimate),]; W4 <- c(w4$estimate[1], w4$se[1])

cell <- function(v){ if (is.null(v) || all(is.na(v))) return("---"); est<-v[1]; se<-v[2]
  star <- if (!is.na(se) && se>0 && abs(est)/se > z95) "^{*}" else ""
  sprintf("$%.3f%s$ (%.3f)", est, star, se) }
trow <- function(lbl,d,dr) sprintf("%s & %s & %s \\\\", lbl, cell(d), cell(dr))

# ---- Table 2: specification ladder ------------------------------------------
A <- list(
  list("P\\&P specification (additive FE, covariates, entropy weights)", sgk("additive","4 P&P","weighted"), NULL),
  list("\\quad Drop entropy weights", sgk("additive","4 P&P","unweighted"), NULL),
  list("\\quad Drop covariates", sgk("additive","none","unweighted"), NULL),
  list("\\quad Agency $\\times$ stack FE", att("stacked|deck"), att("stacked|droppers")),
  list("\\quad $\\pm 4$ window, non-changer controls", W4, NULL),
  list("\\quad $\\pm 4$ window, never-treated controls only", att("Cengiz|deck"), att("Cengiz|droppers")))
B <- list(
  list("Callaway--Sant'Anna, never-treated controls", att("CS-NT|deck"), att("CS-NT|droppers")),
  list("Callaway--Sant'Anna, not-yet-treated controls", att("CS-NYT|deck"), att("CS-NYT|droppers")),
  list("Imputation", NULL, att("imputation|droppers")))
L <- c("\\begin{tabular}{lcc}","\\toprule","Specification & Deck sample & Droppers sample \\\\","\\midrule",
       "\\multicolumn{3}{@{}l}{\\emph{Panel A: Stacked regression}} \\\\")
for (x in A) L <- c(L, trow(x[[1]], x[[2]], x[[3]]))
L <- c(L, "\\addlinespace", "\\multicolumn{3}{@{}l}{\\emph{Panel B: Alternative estimators (simple aggregation)}} \\\\")
for (x in B) L <- c(L, trow(x[[1]], x[[2]], x[[3]]))
L <- c(L, "\\bottomrule","\\end{tabular}","")
notes2 <- paste0("\\noindent Outcome is an indicator for any fatal encounter. Standard errors clustered by agency ",
"in parentheses; imputation standard error from a cluster bootstrap. $^{*}p<0.05$. The droppers sample ",
"excludes the six agencies that adopted a requirement and the two absorbed reversals; imputation requires ",
"absorbing treatment and is estimated on this sample only. All rows after the first use no covariates or ",
"weights unless noted; Panel B uses simple aggregation.")
L <- c(L, "\\smallskip", paste0("{\\footnotesize ", notes2, "\\par}"))
writeLines(L, file.path(TAB,"table_ladder.tex"))

# ---- Table 3: pre-tests, stacked comparison group only, droppers ------------
pr <- rd("d6_pretests_all.csv"); pr <- pr[grepl("stacked", pr$comparison, ignore.case=TRUE),]
md <- rd("d6_mde_wcr.csv"); md <- md[md$estimator=="stacked",]
fw <- rd("fwl_decomp_unweighted.csv")                 # deck w_s
pfmt <- function(x){ x<-as.numeric(x); if (x<0.001) "$<0.001$" else sprintf("%.3f", x) }
T <- c("\\begin{tabular}{lccccccc}","\\toprule",
       "& & & \\multicolumn{3}{c}{Slope test $p$-value} & & \\\\","\\cmidrule(lr){4-6}",
       "Cohort & Treated clusters & Stacked weight & CR1 & WCR & RI & MDE (WCR) & Tolerance ($\\theta^{*}{=}0.10$) \\\\","\\midrule")
for (c in c(2009,2013,2017)){
  p <- pr[pr$cohort==c,]; m <- md[md$cohort==c,]
  pre_g <- m$pre_g[1]; n_g <- m$n_g[1]; post_g <- 2020-c+1; a_g <- (pre_g+post_g)/2; tol <- 0.10/a_g
  ws <- fw$w_s[fw$stack==c]
  T <- c(T, sprintf("%d & %d & %.3f & %s & %s & %s & %.3f & %.4f \\\\", c, n_g, ws,
    pfmt(p$slope_p_CR1[1]), pfmt(p$slope_p_WCR[1]), pfmt(p$slope_p_RIunm[1]), m$mde_wcr[1], tol)) }
T <- c(T, "\\bottomrule","\\end{tabular}","")
notes3 <- paste0("\\noindent Linear pre-trend tests against each stack's clean controls, droppers sample. ",
"CR1: cluster-robust; WCR: restricted wild cluster bootstrap with Webb weights; RI: randomization ",
"inference over never-treated placebo cohorts. MDE is the minimum detectable trend per year at ",
"$\\alpha = 0.05$ and 80 percent power. Tolerance is the break-even slope $\\theta^{*}/a_g$ ",
"(Equation~\\ref{eq:breakeven}). Stacked weights are from the deck sample (Table~\\ref{tab:stack_decomp}).")
T <- c(T, "\\smallskip", paste0("{\\footnotesize ", notes3, "\\par}"))
writeLines(T, file.path(TAB,"table_pretests.tex"))
cat("Wrote tables/table_ladder.tex, tables/table_pretests.tex\n")
