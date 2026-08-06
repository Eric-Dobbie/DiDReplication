# =============================================================================
# stacked_dimensionality.R
#
# Summarize the dimensionality (row count) of the stacked dataset under
# different window (+/-k) and clean-control-group assumptions. Balanced Cengiz
# construction: each treated cohort g keeps years [g-k, g+k] (2k+1 rows/unit),
# eligible only if the window is fully observed (g-k >= 2000 & g+k <= 2020);
# the two reversible agencies are dropped. Writes a LaTeX table.
# =============================================================================

suppressMessages(library(dplyr))
DATA_DIR <- Sys.getenv("DIDREP_DATA", unset = "data")
OUT_DIR  <- Sys.getenv("DIDREP_OUT",  unset = "output")
FIRST_OUT <- 2000; LAST_OUT <- 2020
REVERSIBLE <- c("memphis tennessee", "portsmouth new hampshire")

dta <- read.csv(file.path(DATA_DIR, "dta.csv"))
dta <- dta[!dta$agency.id %in% REVERSIBLE, ]
nochange <- dta[dta$change.type == "No Change", ]
allcoh   <- sort(unique(dta$year.changed[!is.na(dta$year.changed) & dta$year.changed >= FIRST_OUT]))

n_rows <- function(control, k) {
  cohorts <- allcoh[allcoh - k >= FIRST_OUT & allcoh + k <= LAST_OUT]
  rows <- 0L; ncoh <- 0L
  for (g in cohorts) {
    tr <- dta[!is.na(dta$year.changed) & dta$year.changed == g &
                dta$year >= g - k & dta$year <= g + k, ]
    if (nrow(tr) == 0) next
    ncoh <- ncoh + 1L
    ct <- nochange
    if (control == "not-yet-treated") {
      future <- dta[!is.na(dta$year.changed) & dta$year.changed > g + k, ]
      ct <- rbind(ct, future)
    }
    ct <- ct[ct$year >= g - k & ct$year <= g + k, ]
    rows <- rows + nrow(tr) + nrow(ct)
  }
  data.frame(control = control, k = k, cohorts = ncoh, rows = rows)
}

grid <- expand.grid(control = c("never-treated", "not-yet-treated"),
                    k = 2:10, stringsAsFactors = FALSE)
res <- do.call(rbind, Map(n_rows, grid$control, grid$k))
res <- res[order(res$control, res$k), ]
print(res, row.names = FALSE)
res <- res[res$rows > 0, ]                       # drop empty windows (k >= 10)
cap1 <- function(s) paste0(toupper(substr(s, 1, 1)), substr(s, 2, nchar(s)))

# ---- write LaTeX table -----------------------------------------------------
fmt <- function(x) formatC(x, format = "d", big.mark = ",")
lines <- c(
  "\\begin{table}[t]",
  "\\centering",
  "\\caption{Dimensionality of the stacked dataset under alternative window and",
  "clean-control-group assumptions. Each treated cohort $g$ keeps a balanced",
  "$\\pm k$ event window ($2k+1$ rows per unit), and is included only if the",
  "window is fully observed ($g-k \\ge 2000$ and $g+k \\le 2020$); the two",
  "reversible agencies are dropped. Windows wider than $\\pm 9$ leave no",
  "fully-observed cohort. For reference, the Payson \\& Parinandi full-panel",
  "dataset (never-treated, no window) has 278{,}324 rows.}",
  "\\label{tab:stacked-dim}",
  "\\begin{tabular}{llr}",
  "\\toprule",
  "Control group & Window (yrs before / after) & Rows \\\\",
  "\\midrule")
prev <- ""
for (i in seq_len(nrow(res))) {
  cg <- if (res$control[i] != prev) res$control[i] else ""
  prev <- res$control[i]
  if (cg == "" && i > 1 && res$control[i] != res$control[i-1]) NULL
  wcell <- sprintf("%d / %d", res$k[i], res$k[i])
  lines <- c(lines, sprintf("%s & %s & %s \\\\",
                            cap1(cg), wcell, fmt(res$rows[i])))
  if (i < nrow(res) && res$control[i] != res$control[i+1]) lines <- c(lines, "\\midrule")
}
lines <- c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(lines, file.path(OUT_DIR, "stacked_dimensionality.tex"))
write.csv(res, file.path(OUT_DIR, "stacked_dimensionality.csv"), row.names = FALSE)
cat("\nWrote", file.path(OUT_DIR, "stacked_dimensionality.tex"), "and .csv\n")
