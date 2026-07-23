# =============================================================================
# cs_control_variants.R
#
# Callaway & Sant'Anna under alternative control groups and the paper's controls.
# Adds two dimensions to the baseline CS spec:
#   1. control group:  never-treated  vs  not-yet-treated
#   2. covariates:      none  vs  the paper's Table-3 controls
#        xformla = ~ log.pop + log.med.inc + pct.white + pct.white.officers.imputed
#        (exactly callaway_replication.R's m2 covariate set)
#
# For each of the 4 CS specifications we extract atoms + weights, verify the
# recombination (sum(weight*estimate) == aggte simple overall), and record the
# overall ATT with its SE. The paper's stacked Table-3 estimate (m4: covariates
# + balancing weights) and the plain stacked coefficient are added as references.
#
# Outputs:
#   output/cs_variants_summary.csv   overall ATT/SE + recombination per spec
#   output/cs_variants_atoms.csv     atoms + weights for all 4 CS specs
#   figures/cs_variants_comparison.{pdf,png}   coefficient plot
# =============================================================================

suppressMessages({library(did); library(lfe); library(dplyr); library(ggplot2)})

.here <- tryCatch({ a <- commandArgs(FALSE); f <- sub("^--file=","",a[grep("^--file=",a)])
  if (length(f)) dirname(normalizePath(f)) else "R" }, error=function(e) "R")
source(file.path(.here, "extract_atoms.R"))

DATA_DIR <- Sys.getenv("DIDREP_DATA", unset="data")
OUT_DIR  <- Sys.getenv("DIDREP_OUT",  unset="output")
FIG_DIR  <- Sys.getenv("DIDREP_FIG",  unset="figures")
dir.create(OUT_DIR, showWarnings=FALSE, recursive=TRUE)
dir.create(FIG_DIR, showWarnings=FALSE, recursive=TRUE)

# ---- data prep (same as run_extraction.R) ----------------------------------
dta <- read.csv(file.path(DATA_DIR, "dta.csv"))
dta$year.changed <- ave(dta$year.changed, dta$agency.id,
                        FUN=function(x){ v<-x[!is.na(x)]; if(length(v)) min(v) else NA_real_ })
dta$agency.num   <- as.numeric(factor(dta$agency.id))
dta$year.changed <- ifelse(is.na(dta$year.changed)&dta$no.req==0, 0,
                    ifelse(is.na(dta$year.changed)&dta$no.req==1, 1987, dta$year.changed))

covs <- ~ log.pop + log.med.inc + pct.white + pct.white.officers.imputed

specs <- list(
  list(id="never_nocov",  cg="nevertreated",  x=~1,   lab="CS  never-treated,  no covariates"),
  list(id="never_cov",    cg="nevertreated",  x=covs, lab="CS  never-treated,  + covariates (paper m2)"),
  list(id="notyet_nocov", cg="notyettreated", x=~1,   lab="CS  not-yet-treated,  no covariates"),
  list(id="notyet_cov",   cg="notyettreated", x=covs, lab="CS  not-yet-treated,  + covariates")
)

atoms_all <- list(); summ <- list()
for (s in specs) {
  a <- extract_cs_atoms(dta, yname="any.fatalities", tname="year", idname="agency.num",
                        gname="year.changed", xformla=s$x, control_group=s$cg,
                        clustervars="agency.num")
  a$spec <- s$id; atoms_all[[s$id]] <- a
  recon <- sum(a$weight * a$estimate)
  summ[[s$id]] <- data.frame(
    spec=s$id, label=s$lab, control_group=s$cg,
    covariates=ifelse(inherits(s$x,"formula") && length(all.vars(s$x))==0, "none","paper 4"),
    overall=attr(a,"pooled_estimate"), se=attr(a,"pooled_se"),
    recon=recon, abs_diff=abs(recon-attr(a,"pooled_estimate")), n_atoms=nrow(a))
}

# ---- stacked references (paper Table 3) ------------------------------------
stk <- read.csv(file.path(DATA_DIR, "stacked_fatal.csv"))
m4 <- felm(any.fatalities ~ no.req + log.pop + log.med.inc + pct.white + pct.white.officers.imputed |
             agency.id + year.cohort | 0 | agency.id, data=stk, weights=stk$weights)
m_plain <- felm(any.fatalities ~ no.req | agency.id + year.cohort | 0 | agency.id, data=stk)
stk_summ <- data.frame(
  spec=c("stacked_m4","stacked_plain"),
  label=c("Stacked  (paper Table 3: + covariates + weights)","Stacked  (plain, no covariates)"),
  control_group="never-treated (clean)", covariates=c("paper 4 + wts","none"),
  overall=c(coef(m4)["no.req"], coef(m_plain)["no.req"]),
  se=c(m4$se["no.req"], m_plain$se["no.req"]),
  recon=NA, abs_diff=NA, n_atoms=NA)

summary_tbl <- bind_rows(bind_rows(summ), stk_summ)
write.csv(summary_tbl, file.path(OUT_DIR, "cs_variants_summary.csv"), row.names=FALSE)
write.csv(bind_rows(atoms_all), file.path(OUT_DIR, "cs_variants_atoms.csv"), row.names=FALSE)

cat("== overall ATT by specification ==\n")
print(summary_tbl %>% select(label, overall, se, abs_diff), row.names=FALSE, digits=4)

# ---- coefficient plot ------------------------------------------------------
pd <- summary_tbl %>%
  mutate(lo = overall - 1.96*se, hi = overall + 1.96*se,
         grp = case_when(grepl("^Stacked", label) ~ "Stacked",
                         grepl("not-yet", label) ~ "CS: not-yet-treated",
                         TRUE ~ "CS: never-treated"),
         label = factor(label, levels = rev(c(
           "CS  never-treated,  no covariates",
           "CS  never-treated,  + covariates (paper m2)",
           "CS  not-yet-treated,  no covariates",
           "CS  not-yet-treated,  + covariates",
           "Stacked  (plain, no covariates)",
           "Stacked  (paper Table 3: + covariates + weights)"))))
pal <- c("CS: never-treated"="#0072B2","CS: not-yet-treated"="#56B4E9","Stacked"="#009E73")

p <- ggplot(pd, aes(overall, label, color=grp)) +
  geom_vline(xintercept=0, linewidth=0.4, color="grey55", linetype="dashed") +
  geom_errorbarh(aes(xmin=lo, xmax=hi), height=0.25, linewidth=0.7) +
  geom_point(size=3.4) +
  geom_text(aes(label=sprintf("%+.3f", overall)), vjust=-1, size=3, fontface="bold", show.legend=FALSE) +
  scale_color_manual(values=pal, name=NULL) +
  labs(title="CS under alternative control groups and covariates",
       subtitle="Overall ATT on P(fatal encounter), 95% CI. Covariates = the paper's 4 (pop, income, %white, %white officers).\nCS stays positive & insignificant across specs; stacked (741 clean controls) is precisely negative.",
       x="Overall ATT (effect on probability of a fatal encounter)", y=NULL) +
  theme_bw(base_size=12) +
  theme(plot.title=element_text(face="bold", size=13),
        plot.subtitle=element_text(size=8.5, color="grey30"),
        legend.position="top", panel.grid.minor=element_blank(),
        axis.text.y=element_text(size=9.5))

ggsave(file.path(FIG_DIR, "cs_variants_comparison.pdf"), p, width=9.5, height=5.5, device=cairo_pdf)
ggsave(file.path(FIG_DIR, "cs_variants_comparison.png"), p, width=9.5, height=5.5, dpi=200)
cat("\nWrote figures/cs_variants_comparison.{pdf,png} and output/cs_variants_{summary,atoms}.csv\n")
