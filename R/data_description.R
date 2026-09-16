# =============================================================================
# data_description.R
#
# Descriptive summary of the Payson & Parinandi (2024) fatal-encounters panel
# (dta.csv): years covered, cities, cohorts, and sample size. Writes a LaTeX
# table (and a CSV) to DIDREP_OUT. All figures are computed from the data.
# =============================================================================
DATA_DIR <- Sys.getenv("DIDREP_DATA", unset = "data")
OUT_DIR  <- Sys.getenv("DIDREP_OUT",  unset = "output")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

dta <- read.csv(file.path(DATA_DIR, "dta.csv"))
st  <- tryCatch(read.csv(file.path(DATA_DIR, "stacked_fatal.csv")), error = function(e) NULL)

yrs    <- range(dta$year)
nyr    <- length(unique(dta$year))
ncity  <- length(unique(dta$agency.id))
nrow_t <- nrow(dta)
obs    <- !is.na(dta$any.fatalities)
nobs   <- sum(obs)
oyrs   <- range(dta$year[obs])
first_type <- dta$change.type[!duplicated(dta$agency.id)]
n_drop <- sum(first_type == "Dropped")
n_adopt<- sum(first_type == "Adopted")
n_none <- sum(first_type == "No Change")
n_change <- n_drop + n_adopt
coh_all <- sort(unique(dta$year.changed[!is.na(dta$year.changed)]))
coh_obs <- coh_all[coh_all >= 2000]
st_rows <- if (!is.null(st)) nrow(st) else NA
st_obs  <- if (!is.null(st)) sum(!is.na(st$any.fatalities)) else NA

# ---- CSV -------------------------------------------------------------------
tab <- data.frame(
  quantity = c("Years covered", "Outcome observed", "Cities (agencies)",
               "Panel rows (total)", "Panel rows (observed outcome)",
               "Cohorts (change years, total)", "Cohorts (observable, >=2000)",
               "Cities that dropped a requirement", "Cities that adopted a requirement",
               "Cities never changing", "Stacked-file rows", "Stacked-file rows (observed)"),
  value = c(sprintf("%d-%d (%d yrs)", yrs[1], yrs[2], nyr),
            sprintf("%d-%d (%d yrs)", oyrs[1], oyrs[2], length(oyrs[1]:oyrs[2])),
            ncity, nrow_t, nobs, length(coh_all), length(coh_obs),
            n_drop, n_adopt, n_none, st_rows, st_obs),
  stringsAsFactors = FALSE)
print(tab, row.names = FALSE)
write.csv(tab, file.path(OUT_DIR, "data_description.csv"), row.names = FALSE)

# ---- LaTeX -----------------------------------------------------------------
fmt <- function(x) formatC(x, format = "d", big.mark = ",")
L <- c(
  "\\begin{table}[t]\\centering",
  "\\caption{Description of the Payson \\& Parinandi (2024) fatal-encounters panel",
  "(\\texttt{dta.csv}). The outcome, an indicator for any fatal civilian encounter,",
  "is observed only from 2000; the panel is balanced across all cities and years.",
  "Cohorts are policy-change years; only the 15 with a change year $\\ge 2000$ are",
  "observable in the outcome window. Two cities (Memphis, Portsmouth) both adopted",
  "and dropped a requirement, so they contribute to two cohorts each.}",
  "\\label{tab:data-description}",
  "\\begin{tabular}{lr}",
  "\\toprule",
  "\\multicolumn{2}{l}{\\textit{Coverage}} \\\\",
  sprintf("\\quad Years covered & %d--%d (%d years) \\\\", yrs[1], yrs[2], nyr),
  sprintf("\\quad Outcome observed & %d--%d (%d years) \\\\", oyrs[1], oyrs[2], length(oyrs[1]:oyrs[2])),
  sprintf("\\quad Cities (agencies) & %s \\\\", fmt(ncity)),
  "\\midrule",
  "\\multicolumn{2}{l}{\\textit{Sample size}} \\\\",
  sprintf("\\quad Panel rows (total, balanced) & %s \\\\", fmt(nrow_t)),
  sprintf("\\quad Panel rows (observed outcome) & %s \\\\", fmt(nobs)),
  sprintf("\\quad Stacked-file rows (total) & %s \\\\", fmt(st_rows)),
  sprintf("\\quad Stacked-file rows (observed outcome) & %s \\\\", fmt(st_obs)),
  "\\midrule",
  "\\multicolumn{2}{l}{\\textit{Cohorts and treatment}} \\\\",
  sprintf("\\quad Cohorts (change years, total) & %d \\\\", length(coh_all)),
  sprintf("\\quad Cohorts (observable, $\\ge 2000$) & %d \\\\", length(coh_obs)),
  sprintf("\\quad Cities that dropped a requirement & %d \\\\", n_drop),
  sprintf("\\quad Cities that adopted a requirement & %d \\\\", n_adopt),
  sprintf("\\quad Cities never changing & %s \\\\", fmt(n_none)),
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}")
writeLines(L, file.path(OUT_DIR, "data_description.tex"))
cat("\nWrote", file.path(OUT_DIR, "data_description.tex"), "and .csv\n")
