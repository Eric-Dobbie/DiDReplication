# =============================================================================
# cengiz_notyet.R
#
# Cengiz +/-4 stacked design with NOT-YET-TREATED clean controls, compared to
# the never-treated Cengiz stack and to the P&P full-panel stack.
#
# For each cohort g (window [g-4, g+4], 9 rows/unit), the control pool is every
# unit that does NOT change treatment inside the window:
#     never-treated ("No Change")  UNION  future-treated (year.changed > g+4).
# Units that change during the window (year.changed in [g-4, g+4], g excepted)
# are excluded, as are the two reversible agencies. This is the canonical
# Cengiz "clean = not-yet-treated" control group (never-treated is the limiting
# case). Same balanced window, same stack-specific unit+time FE as cengiz_stacked.R.
# =============================================================================

suppressMessages({library(lfe); library(dplyr); library(ggplot2)})
.here <- tryCatch({ a <- commandArgs(FALSE); f <- sub("^--file=","",a[grep("^--file=",a)])
  if (length(f)) dirname(normalizePath(f)) else "R" }, error=function(e) "R")
source(file.path(.here, "extract_atoms.R"))

DATA_DIR <- Sys.getenv("DIDREP_DATA", unset="data")
OUT_DIR  <- Sys.getenv("DIDREP_OUT",  unset="output")
FIG_DIR  <- Sys.getenv("DIDREP_FIG",  unset="figures")

K <- 4; FIRST_OUT <- 2000; LAST_OUT <- 2020
REVERSIBLE <- c("memphis tennessee","portsmouth new hampshire")

dta <- read.csv(file.path(DATA_DIR, "dta.csv"))
dta <- dta[!dta$agency.id %in% REVERSIBLE, ]
cohorts  <- sort(unique(dta$year.changed[!is.na(dta$year.changed)]))
cohorts  <- cohorts[cohorts - K >= FIRST_OUT & cohorts + K <= LAST_OUT]     # 2004..2016
nochange <- dta[dta$change.type == "No Change", ]

# ---- build a Cengiz stack under a chosen control group ---------------------
build_cengiz <- function(control = c("never","notyet")) {
  control <- match.arg(control)
  parts <- list()
  for (g in cohorts) {
    tr <- dta[!is.na(dta$year.changed) & dta$year.changed == g &
                dta$year >= g - K & dta$year <= g + K, ]
    if (nrow(tr) == 0) next
    tr$cohort <- g; tr$treat <- 1L
    ct <- nochange
    if (control == "notyet") {
      future <- dta[!is.na(dta$year.changed) & dta$year.changed > g + K, ]  # not-yet-treated
      ct <- rbind(ct, future)
    }
    ct <- ct[ct$year >= g - K & ct$year <= g + K, ]
    ct$cohort <- g; ct$treat <- 0L
    parts[[length(parts) + 1L]] <- rbind(tr, ct)
  }
  cz <- do.call(rbind, parts)
  cz$stackunit <- paste0(cz$agency.id, "__", cz$cohort)
  cz$stacktime <- paste0(cz$year, "__", cz$cohort)
  cz
}

run_one <- function(control) {
  cz <- build_cengiz(control)
  a  <- extract_stacked_atoms(cz, yname="any.fatalities", idname="stackunit",
                              stacktime="stacktime", stackid="cohort",
                              treatname="no.req", cluster="agency.id", resid_on="estimation")
  list(atoms=a, pooled=attr(a,"pooled_estimate"), se=attr(a,"pooled_se"),
       n_ctrl=sum(cz$treat==0)/9, recon=sum(a$weight*a$estimate))
}

never  <- run_one("never")
notyet <- run_one("notyet")

cat("== Cengiz +/-4, control-group comparison ==\n")
cat(sprintf("never-treated  : %+.5f (se %.5f)  controls/stack~%.0f  recomb-diff %.1e\n",
            never$pooled, never$se, never$n_ctrl/length(cohorts), abs(never$recon-never$pooled)))
cat(sprintf("not-yet-treated: %+.5f (se %.5f)  controls/stack~%.0f  recomb-diff %.1e\n",
            notyet$pooled, notyet$se, notyet$n_ctrl/length(cohorts), abs(notyet$recon-notyet$pooled)))

cat("\ncohort | never | not-yet\n")
m <- merge(never$atoms[c("group","estimate")], notyet$atoms[c("group","estimate")],
           by="group", suffixes=c("_never","_notyet"))
for (i in order(m$group)) cat(sprintf("  %d | %+.4f | %+.4f\n", m$group[i], m$estimate_never[i], m$estimate_notyet[i]))

# ---- references + coefficient plot -----------------------------------------
pp <- read.csv(file.path(DATA_DIR, "stacked_fatal.csv"))
m4 <- felm(any.fatalities ~ no.req + log.pop + log.med.inc + pct.white + pct.white.officers.imputed |
             agency.id + year.cohort | 0 | agency.id, data=pp, weights=pp$weights)
pl <- felm(any.fatalities ~ no.req | agency.id + year.cohort | 0 | agency.id, data=pp)

pd <- data.frame(
  label = c("Stacked: P&P full-panel (Table 3, +cov+wts)",
            "Stacked: P&P full-panel (plain)",
            "Cengiz +/-4: never-treated controls",
            "Cengiz +/-4: not-yet-treated controls"),
  est = c(coef(m4)["no.req"], coef(pl)["no.req"], never$pooled, notyet$pooled),
  se  = c(m4$se["no.req"], pl$se["no.req"], never$se, notyet$se),
  grp = c("P&P full-panel","P&P full-panel","Cengiz window","Cengiz window"))
pd$lo <- pd$est - 1.96*pd$se; pd$hi <- pd$est + 1.96*pd$se
pd$label <- factor(pd$label, levels=rev(pd$label))
pal <- c("P&P full-panel"="#009E73","Cengiz window"="#CC79A7")

p <- ggplot(pd, aes(est, label, color=grp)) +
  geom_vline(xintercept=0, linetype="dashed", color="grey55", linewidth=0.4) +
  geom_errorbarh(aes(xmin=lo, xmax=hi), height=0.22, linewidth=0.7) +
  geom_point(size=3.4) +
  geom_text(aes(label=sprintf("%+.3f", est)), vjust=-1, size=3, fontface="bold", show.legend=FALSE) +
  scale_color_manual(values=pal, name=NULL) +
  labs(title="Stacked estimate: window and clean-control choice",
       subtitle="Pooled coefficient on 'requirement dropped' (P(fatal encounter)), 95% CI. Cengiz = balanced +/-4 window, 9 rows/unit.",
       x="Pooled stacked estimate", y=NULL) +
  theme_bw(base_size=12) +
  theme(plot.title=element_text(face="bold", size=13),
        plot.subtitle=element_text(size=8.5, color="grey30"),
        legend.position="top", panel.grid.minor=element_blank(),
        axis.text.y=element_text(size=9.5))
ggsave(file.path(FIG_DIR, "cengiz_control_groups.pdf"), p, width=9.5, height=4.6, device=cairo_pdf)
ggsave(file.path(FIG_DIR, "cengiz_control_groups.png"), p, width=9.5, height=4.6, dpi=200)

notyet$atoms$design <- "cengiz_notyet"
write.csv(notyet$atoms, file.path(OUT_DIR, "cengiz_notyet_atoms.csv"), row.names=FALSE)
write.csv(pd[c("label","est","se")], file.path(OUT_DIR, "cengiz_control_groups_summary.csv"), row.names=FALSE)
cat("\nWrote figures/cengiz_control_groups.{pdf,png} and output/cengiz_notyet_atoms.csv, cengiz_control_groups_summary.csv\n")
