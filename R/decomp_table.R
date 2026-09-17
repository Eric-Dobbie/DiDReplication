# =============================================================================
# decomp_table.R
#
# LaTeX table of the stacked-regression decomposition: per-cohort effect
# coefficients (beta_s) and aggregation weights (w_s) under the corrected
# interactive FE (agency x stack + year x cohort), for both the unweighted R2
# spec and the entropy-balancing-weighted spec. Reads the committed FWL outputs
# (fwl_decomp_{unweighted,weighted}.csv from fwl_decomp.R) and reports the
# reconstruction total sum_s w_s * beta_s next to each package pooled estimate.
#
# Cohorts with no within-stack residual variation (V_s ~ 0) are unidentified
# (beta_s undefined, weight 0) and shown as "--".
# =============================================================================
OUT_DIR <- Sys.getenv("DIDREP_OUT", unset = "output")
uw <- read.csv(file.path(OUT_DIR, "fwl_decomp_unweighted.csv"))
wt <- read.csv(file.path(OUT_DIR, "fwl_decomp_weighted.csv"))

TOL <- 1e-10                                    # treat w_s below this as unidentified
merge_cols <- function(d) {
  ident <- d$w_s > TOL
  data.frame(stack = d$stack,
             beta = ifelse(ident, d$beta_s, NA_real_),
             w    = ifelse(ident, d$w_s,   NA_real_))
}
U <- merge_cols(uw); W <- merge_cols(wt)
stopifnot(identical(U$stack, W$stack))

pooled_uw <- -0.098915; pooled_wt <- -0.108393                       # package pooled (R2 / step 7)
recon_uw  <- sum(U$beta * U$w, na.rm = TRUE)
recon_wt  <- sum(W$beta * W$w, na.rm = TRUE)

cell <- function(x, dig) ifelse(is.na(x), "--", sprintf(paste0("%.", dig, "f"), x))
L <- c(
  "\\begin{table}[t]\\centering",
  "\\caption{Stacked-regression decomposition: per-cohort effect coefficients",
  "$\\beta_s$ and variance-share weights $w_s$ under the corrected interactive",
  "fixed effects ($\\text{agency}\\times\\text{stack}+\\text{year}\\times\\text{cohort}$),",
  "for the unweighted R2 spec and the entropy-balancing-weighted spec. Weights are",
  "Frisch--Waugh--Lovell variance shares (residualized treatment). Cohorts with no",
  "within-stack residual variation are unidentified ($\\beta_s$ undefined, $w_s=0$),",
  "shown as ``--''; 2000/2001/2003 have no clean controls and no reverse-direction",
  "contrast, and 2016's lone treated unit is absorbed. The pooled coefficient equals",
  "$\\sum_s w_s\\beta_s$ to machine precision (see \\texttt{fwl\\_decomp.R}).}",
  "\\label{tab:decomp}",
  "\\begin{tabular}{lrrrr}",
  "\\toprule",
  "& \\multicolumn{2}{c}{Unweighted (R2)} & \\multicolumn{2}{c}{Weighted (ebal)} \\\\",
  "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}",
  "Cohort & $\\beta_s$ & $w_s$ & $\\beta_s$ & $w_s$ \\\\",
  "\\midrule")
for (i in seq_len(nrow(U))) {
  L <- c(L, sprintf("%d & %s & %s & %s & %s \\\\",
                    U$stack[i], cell(U$beta[i], 4), cell(U$w[i], 4),
                    cell(W$beta[i], 4), cell(W$w[i], 4)))
}
L <- c(L, "\\midrule",
  sprintf("Pooled $\\sum_s w_s\\beta_s$ & \\multicolumn{2}{c}{%.4f} & \\multicolumn{2}{c}{%.4f} \\\\",
          recon_uw, recon_wt),
  "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(L, file.path(OUT_DIR, "decomp_table.tex"))

# console echo
out <- data.frame(cohort = U$stack,
                  beta_unwt = round(U$beta, 4), w_unwt = round(U$w, 4),
                  beta_wtd  = round(W$beta, 4), w_wtd  = round(W$w, 4))
print(out, row.names = FALSE)
cat(sprintf("\nPooled sum w*beta: unweighted %.6f (pkg %.6f) | weighted %.6f (pkg %.6f)\n",
            recon_uw, pooled_uw, recon_wt, pooled_wt))
cat("Wrote", file.path(OUT_DIR, "decomp_table.tex"), "\n")
