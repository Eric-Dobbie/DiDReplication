# =============================================================================
# stack_composition.R
#
# Composition of each stack in the shipped stacked_fatal.csv: the treated cohort
# members split by which way they switch no.req, and the controls split by
# treatment history. Makes the identification of each stack -- especially the
# clean-control-free cohort 2002 -- explicit. Writes a LaTeX table + CSV.
#
#   treated (treat==1):
#     drop     : no.req 0 -> 1  (the canonical "requirement dropped" treatment)
#     adopt    : no.req 1 -> 0  (a changer switching THE OTHER WAY)
#     absorbed : reversal on a segment boundary -> no within-stack switch
#   controls (treat==0):
#     never  : no.req == 0 always (always has a requirement -- clean DiD control)
#     always : no.req == 1 always (always no requirement -- always in treated state)
#     notyet : no.req switches inside the stack (a later dropper used early)
# =============================================================================
DATA_DIR <- Sys.getenv("DIDREP_DATA", unset = "data")
OUT_DIR  <- Sys.getenv("DIDREP_OUT",  unset = "output")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
st <- read.csv(file.path(DATA_DIR, "stacked_fatal.csv"))

dir_treated <- function(v) {
  if (sum(diff(v) != 0) == 0) "absorbed"
  else if (v[1] == 0 & v[length(v)] == 1) "drop"
  else if (v[1] == 1 & v[length(v)] == 0) "adopt"
  else "other"
}
class_ctrl <- function(v) if (all(v == 0)) "never" else if (all(v == 1)) "always" else "notyet"

rows <- list()
for (g in sort(unique(st$cohort))) {
  o  <- st[st$cohort == g, ]; o <- o[order(o$agency.id, o$year), ]
  tr <- o[o$treat == 1, ]; ct <- o[o$treat == 0, ]
  td <- if (nrow(tr)) tapply(tr$no.req, tr$agency.id, dir_treated) else character(0)
  cc <- if (nrow(ct)) tapply(ct$no.req, ct$agency.id, class_ctrl) else character(0)
  d <- sum(td=="drop"); a <- sum(td=="adopt"); ab <- sum(td=="absorbed")
  nv <- sum(cc=="never"); al <- sum(cc=="always"); ny <- sum(cc=="notyet")
  # identification note
  has_ctrl <- (nv + al + ny) > 0
  n_tr_sw  <- d + a                       # treated units with a within-stack switch
  ident <- if (has_ctrl && n_tr_sw > 0) "vs.\\ controls"
           else if (has_ctrl && n_tr_sw == 0) "-- (treated unit absorbed)"
           else if (!has_ctrl && d > 0 && a > 0) "within-treated: drop vs.\\ adopt"
           else "-- (no contrast)"
  rows[[length(rows)+1]] <- data.frame(stack=g, drop=d, adopt=a, absorbed=ab,
    never=nv, always=al, notyet=ny, identified=ident, stringsAsFactors=FALSE)
}
tab <- do.call(rbind, rows)
print(tab, row.names = FALSE)
write.csv(tab, file.path(OUT_DIR, "stack_composition.csv"), row.names = FALSE)

# ---- LaTeX -----------------------------------------------------------------
L <- c(
  "\\begin{table}[t]\\centering",
  "\\caption{Composition of each stack in the shipped \\texttt{stacked\\_fatal.csv}.",
  "Treated cohort members are split by the direction of their \\texttt{no.req}",
  "switch: \\emph{drop} ($0\\!\\to\\!1$, the canonical treatment), \\emph{adopt}",
  "($1\\!\\to\\!0$, a changer switching the other way), and \\emph{absorbed}",
  "(reversal on a segment boundary, no within-stack switch). Controls are split by",
  "treatment history: never-treated (always a requirement), always-treated (always",
  "no requirement), and not-yet-treated. The four early cohorts (2000--2003) have no",
  "controls; among them only 2002 is identified, purely through the within-treated",
  "contrast between its four droppers and its single adopter (Revere). Cohort 2016",
  "has controls but its lone treated unit is absorbed, so it too is unidentified.}",
  "\\label{tab:stack-composition}",
  "\\begin{tabular}{lrrrrrrl}",
  "\\toprule",
  "& \\multicolumn{3}{c}{Treated} & \\multicolumn{3}{c}{Controls} & \\\\",
  "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}",
  "Stack & Drop & Adopt & Absorbed & Never & Always & Not-yet & Identified via \\\\",
  "\\midrule")
for (i in seq_len(nrow(tab))) {
  bold <- tab$stack[i] == 2002
  fmt <- if (bold) "\\textbf{%d} & \\textbf{%d} & \\textbf{%d} & %d & %d & %d & %d & %s \\\\"
         else      "%d & %d & %d & %d & %d & %d & %d & %s \\\\"
  L <- c(L, sprintf(fmt, tab$stack[i], tab$drop[i], tab$adopt[i], tab$absorbed[i],
                    tab$never[i], tab$always[i], tab$notyet[i], tab$identified[i]))
}
L <- c(L, "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(L, file.path(OUT_DIR, "stack_composition.tex"))
cat("\nWrote output/stack_composition.tex and .csv\n")
