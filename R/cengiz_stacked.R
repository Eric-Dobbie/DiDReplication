# =============================================================================
# cengiz_stacked.R
#
# A literal Cengiz, Dube, Lindner & Zipperer (2019) stacked design, contrasted
# with Payson & Parinandi's full-panel stack.
#
# Cengiz features implemented here:
#   * SYMMETRIC EVENT WINDOW: each sub-experiment keeps only years [g-4, g+4],
#     i.e. exactly 9 rows per unit (4 pre, treatment year, 4 post). Balanced.
#   * Only cohorts whose full window is observed: g-4 >= 2000 AND g+4 <= 2020
#     (outcome any.fatalities exists 2000-2020) -> cohorts 2004..2016.
#   * CLEAN never-treated controls (the 741 "No Change" agencies), windowed.
#   * STACK-SPECIFIC unit AND time fixed effects: agency x cohort + year x cohort
#     (id^stack + time^stack) -- the defining Cengiz FE structure. (P&P instead
#     used a shared agency FE + year.cohort.)
#   * The two reversible agencies (memphis, portsmouth) are dropped: their two
#     events sit < 8 years apart so their +-4 windows overlap, which would reuse
#     the same rows across stacks. (The paper drops these two as well.)
#
# We then run the pooled stacked coefficient, decompose it into cohort atoms +
# FWL weights (extract_stacked_atoms), check the recombination, and compare to
# the P&P full-panel stack and to CS/SA.
# =============================================================================

suppressMessages({library(lfe); library(dplyr); library(ggplot2)})
.here <- tryCatch({ a <- commandArgs(FALSE); f <- sub("^--file=","",a[grep("^--file=",a)])
  if (length(f)) dirname(normalizePath(f)) else "R" }, error=function(e) "R")
source(file.path(.here, "extract_atoms.R"))

DATA_DIR <- Sys.getenv("DIDREP_DATA", unset="data")
OUT_DIR  <- Sys.getenv("DIDREP_OUT",  unset="output")
FIG_DIR  <- Sys.getenv("DIDREP_FIG",  unset="figures")
dir.create(OUT_DIR, showWarnings=FALSE, recursive=TRUE); dir.create(FIG_DIR, showWarnings=FALSE, recursive=TRUE)

K <- 4; FIRST_OUT <- 2000; LAST_OUT <- 2020
REVERSIBLE <- c("memphis tennessee","portsmouth new hampshire")

dta <- read.csv(file.path(DATA_DIR, "dta.csv"))
dta <- dta[!dta$agency.id %in% REVERSIBLE, ]

cohorts <- sort(unique(dta$year.changed[!is.na(dta$year.changed)]))
cohorts <- cohorts[cohorts - K >= FIRST_OUT & cohorts + K <= LAST_OUT]     # 2004..2016
nochange <- dta[dta$change.type == "No Change", ]

# ---- build the windowed stack ----------------------------------------------
parts <- list()
for (g in cohorts) {
  tr <- dta[!is.na(dta$year.changed) & dta$year.changed == g &
              dta$year >= g - K & dta$year <= g + K, ]
  if (nrow(tr) == 0) next
  tr$cohort <- g; tr$treat <- 1L
  ct <- nochange[nochange$year >= g - K & nochange$year <= g + K, ]
  ct$cohort <- g; ct$treat <- 0L
  parts[[length(parts) + 1L]] <- rbind(tr, ct)
}
cz <- do.call(rbind, parts)
cz$stackunit <- paste0(cz$agency.id, "__", cz$cohort)    # id^stack
cz$stacktime <- paste0(cz$year, "__", cz$cohort)         # time^stack

cat("== Cengiz stack ==\n")
cat("cohorts:", paste(cohorts, collapse=","), "\n")
cat("rows:", nrow(cz), " | rows per (stack,unit):",
    paste(range(as.integer(table(cz$stackunit))), collapse="-"), "\n")
cat("treated units:", sum(cz$treat==1)/9, " control appearances:", sum(cz$treat==0)/9, "\n\n")

# ---- decompose (stack-specific FE) -----------------------------------------
atoms <- extract_stacked_atoms(cz, yname="any.fatalities", idname="stackunit",
                               stacktime="stacktime", stackid="cohort",
                               treatname="no.req", cluster="agency.id",
                               resid_on="estimation")
recon  <- sum(atoms$weight * atoms$estimate)
pooled <- attr(atoms, "pooled_estimate")
cat(sprintf("pooled Cengiz stacked coef (no.req) : %+.5f  (se %.5f)\n", pooled, attr(atoms,"pooled_se")))
cat(sprintf("recombination sum(w*beta)           : %+.5f  | abs diff %.2e\n", recon, abs(recon-pooled)))

cat("\ncohort | estimate | weight\n")
for (i in order(atoms$group)) cat(sprintf("  %d | %+.4f | %.4f\n", atoms$group[i], atoms$estimate[i], atoms$weight[i]))

# ---- reference: P&P full-panel stack & shared-FE windowed -------------------
pp <- read.csv(file.path(DATA_DIR, "stacked_fatal.csv"))
b_pp <- coef(felm(any.fatalities ~ no.req | agency.id + year.cohort | 0 | agency.id, data=pp))["no.req"]
# same windowed data but P&P-style shared agency FE (isolates the FE choice)
b_shared <- coef(felm(any.fatalities ~ no.req | agency.id + stacktime | 0 | agency.id, data=cz))["no.req"]
cat(sprintf("\nreference: P&P full-panel stacked    : %+.5f\n", b_pp))
cat(sprintf("reference: windowed, shared agency FE: %+.5f\n", b_shared))

# ---- save atoms + summary --------------------------------------------------
atoms$design <- "cengiz_window9_stackFE"
write.csv(atoms, file.path(OUT_DIR, "cengiz_stacked_atoms.csv"), row.names=FALSE)
summ <- data.frame(
  quantity=c("Cengiz windowed, stack-specific FE","windowed, shared agency FE (P&P-style FE)","P&P full-panel stacked"),
  estimate=c(pooled, b_shared, b_pp))
write.csv(summ, file.path(OUT_DIR, "cengiz_stacked_summary.csv"), row.names=FALSE)

# ---- figure: cohort atoms, Cengiz vs P&P full-panel ------------------------
ppw <- read.csv(file.path(OUT_DIR, "atoms_long.csv")) %>%
  filter(estimator=="stacked", weight>0) %>%
  transmute(group, estimate, weight, design="P&P full-panel")
plotdat <- bind_rows(transmute(atoms, group, estimate, weight, design="Cengiz +/-4 window"), ppw)
ov <- data.frame(design=c("Cengiz +/-4 window","P&P full-panel"), x=c(pooled, b_pp))
pal <- c("Cengiz +/-4 window"="#CC79A7","P&P full-panel"="#009E73")

p <- ggplot(plotdat, aes(estimate, weight, color=design)) +
  geom_vline(xintercept=0, linetype="dashed", color="grey60", linewidth=0.35) +
  geom_vline(data=ov, aes(xintercept=x, color=design), linetype="dotted", linewidth=0.7, show.legend=FALSE) +
  geom_point(size=3, alpha=0.85) +
  ggrepel::geom_text_repel(aes(label=group), size=2.7, fontface="bold",
                           min.segment.length=0, box.padding=0.3, max.overlaps=20,
                           seed=6, show.legend=FALSE) +
  scale_color_manual(values=pal, name=NULL) +
  labs(title="Stacked cohort atoms: literal Cengiz +/-4 window vs. P&P full-panel",
       subtitle="One point per treated cohort (FWL variance-share weight vs. cohort estimate). Dotted lines = pooled stacked coefficient.",
       x="Cohort-level stacked estimate (effect on P(fatal encounter))",
       y="FWL cohort weight") +
  theme_bw(base_size=12) +
  theme(plot.title=element_text(face="bold", size=12.5),
        plot.subtitle=element_text(size=8.7, color="grey30"),
        legend.position="top", panel.grid.minor=element_blank())
ggsave(file.path(FIG_DIR, "cengiz_vs_pp_stacked.pdf"), p, width=9, height=6, device=cairo_pdf)
ggsave(file.path(FIG_DIR, "cengiz_vs_pp_stacked.png"), p, width=9, height=6, dpi=200)
cat("\nWrote figures/cengiz_vs_pp_stacked.{pdf,png} and output/cengiz_stacked_{atoms,summary}.csv\n")
