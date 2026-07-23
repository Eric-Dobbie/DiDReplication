# =============================================================================
# cohort_level_variants.R
#
# Cohort-level view (one point per treated cohort) across control-group and
# covariate variants, for BOTH Callaway & Sant'Anna and Sun & Abraham.
#   - CS  never-treated   controls  (with / without the paper's 4 covariates)
#   - CS  not-yet-treated controls  (with / without covariates)
#   - SA  (with / without covariates)
# Atoms are rolled up to cohort by their own weights:
#   estimate_g = sum_i(w_i*beta_i)/sum_i(w_i),  weight_g = sum_i(w_i).
# Faceted by covariates; the stacked overall is drawn as a per-panel reference.
# =============================================================================

suppressMessages({library(did); library(fixest); library(lfe); library(dplyr); library(ggplot2); library(ggrepel)})

.here <- tryCatch({ a <- commandArgs(FALSE); f <- sub("^--file=","",a[grep("^--file=",a)])
  if (length(f)) dirname(normalizePath(f)) else "R" }, error=function(e) "R")
source(file.path(.here, "extract_atoms.R"))

DATA_DIR <- Sys.getenv("DIDREP_DATA", unset="data")
OUT_DIR  <- Sys.getenv("DIDREP_OUT",  unset="output")
FIG_DIR  <- Sys.getenv("DIDREP_FIG",  unset="figures")
dir.create(FIG_DIR, showWarnings=FALSE, recursive=TRUE)

# ---- data prep -------------------------------------------------------------
dta <- read.csv(file.path(DATA_DIR, "dta.csv"))
dta$year.changed <- ave(dta$year.changed, dta$agency.id,
                        FUN=function(x){ v<-x[!is.na(x)]; if(length(v)) min(v) else NA_real_ })
dta$agency.num   <- as.numeric(factor(dta$agency.id))
dta$year.changed <- ifelse(is.na(dta$year.changed)&dta$no.req==0, 0,
                    ifelse(is.na(dta$year.changed)&dta$no.req==1, 1987, dta$year.changed))
covs   <- c("log.pop","log.med.inc","pct.white","pct.white.officers.imputed")
covfml <- reformulate(covs)

# SA comparable sample (never-treated controls + estimable cohorts)
d_sa <- dta %>% filter(!is.na(any.fatalities)) %>% filter(year.changed==0 | year.changed>=2001)

to_cohort <- function(atoms) atoms %>% filter(weight>0) %>% group_by(group) %>%
  summarize(estimate=sum(weight*estimate)/sum(weight), weight=sum(weight), .groups="drop")

# ---- compute cohort-level estimates for each spec --------------------------
build <- function(cov) {
  xcs <- if (cov) covfml else ~1
  xsa <- if (cov) covs  else NULL
  cs_n <- to_cohort(extract_cs_atoms(dta, "any.fatalities","year","agency.num","year.changed",
                                     xformla=xcs, control_group="nevertreated", clustervars="agency.num"))
  cs_y <- to_cohort(extract_cs_atoms(dta, "any.fatalities","year","agency.num","year.changed",
                                     xformla=xcs, control_group="notyettreated", clustervars="agency.num"))
  sa   <- to_cohort(extract_sa_atoms(d_sa, "any.fatalities","agency.num","year","year.changed",
                                     covars=xsa))
  bind_rows(
    mutate(cs_n, spec="CS: never-treated"),
    mutate(cs_y, spec="CS: not-yet-treated"),
    mutate(sa,   spec="Sun & Abraham")
  ) %>% mutate(cov = if (cov) "+ paper covariates" else "no covariates")
}
dat <- bind_rows(build(FALSE), build(TRUE))
dat$cov  <- factor(dat$cov, levels=c("no covariates","+ paper covariates"))
dat$spec <- factor(dat$spec, levels=c("CS: never-treated","CS: not-yet-treated","Sun & Abraham"))

# ---- stacked reference (per panel): plain vs covariates+weights -------------
stk <- read.csv(file.path(DATA_DIR, "stacked_fatal.csv"))
b_plain <- coef(felm(any.fatalities ~ no.req | agency.id + year.cohort | 0 | agency.id, data=stk))["no.req"]
b_m4    <- coef(felm(any.fatalities ~ no.req + log.pop + log.med.inc + pct.white + pct.white.officers.imputed |
                       agency.id + year.cohort | 0 | agency.id, data=stk, weights=stk$weights))["no.req"]
ref <- data.frame(cov=factor(c("no covariates","+ paper covariates"),
                             levels=levels(dat$cov)),
                  x=c(b_plain, b_m4),
                  lab=c(sprintf("stacked = %+.3f", b_plain), sprintf("stacked (Table 3) = %+.3f", b_m4)))

pal <- c("CS: never-treated"="#0072B2","CS: not-yet-treated"="#56B4E9","Sun & Abraham"="#D55E00")
shp <- c("CS: never-treated"=16,"CS: not-yet-treated"=17,"Sun & Abraham"=5)

# label a few high-weight cohorts (on the CS never-treated series) per panel
lab <- dat %>% filter(spec=="CS: never-treated" & group %in% c(2002,2004,2009,2013,2017))

p <- ggplot(dat, aes(estimate, weight, color=spec, shape=spec)) +
  geom_vline(xintercept=0, linewidth=0.35, color="grey60", linetype="dashed") +
  geom_vline(data=ref, aes(xintercept=x), color="#009E73", linewidth=0.6,
             linetype="dotted", inherit.aes=FALSE) +
  geom_text(data=ref, aes(x=x, label=lab), y=Inf, hjust=-0.03, vjust=1.6,
            color="#009E73", size=2.9, fontface="bold", inherit.aes=FALSE) +
  geom_point(size=2.8, stroke=0.9, alpha=0.9) +
  geom_text_repel(data=lab, aes(label=group), size=2.7, fontface="bold",
                  color="#0072B2", min.segment.length=0, box.padding=0.3,
                  segment.size=0.2, max.overlaps=20, seed=4, show.legend=FALSE) +
  facet_wrap(~cov, ncol=1, scales="free_x") +
  scale_color_manual(values=pal, name=NULL) +
  scale_shape_manual(values=shp, name=NULL) +
  scale_x_continuous(breaks=scales::pretty_breaks(9)) +
  labs(title="Cohort-level estimates by control group and covariates (CS and SA)",
       subtitle="One point per treated cohort; atoms rolled up to cohort by their weights. Green dotted line = stacked overall (per panel).\nCovariates = log.pop, log.med.inc, pct.white, pct.white.officers.imputed.",
       x="Cohort-level estimate (effect on probability of a fatal encounter)",
       y="Cohort aggregation weight") +
  theme_bw(base_size=12) +
  theme(plot.title=element_text(face="bold", size=13),
        plot.subtitle=element_text(size=8.5, color="grey30"),
        strip.text=element_text(face="bold", size=10.5),
        strip.background=element_rect(fill="grey93"),
        legend.position="top", panel.grid.minor=element_blank())

ggsave(file.path(FIG_DIR, "cohort_level_variants.pdf"), p, width=9, height=9, device=cairo_pdf)
ggsave(file.path(FIG_DIR, "cohort_level_variants.png"), p, width=9, height=9, dpi=200)
cat("Wrote figures/cohort_level_variants.{pdf,png}\n")

# quick console check: do CS-never and SA still coincide with covariates?
chk <- dat %>% filter(spec %in% c("CS: never-treated","Sun & Abraham")) %>%
  select(cov, spec, group, estimate) %>%
  tidyr::pivot_wider(names_from=spec, values_from=estimate)
cat("\nmax |CS(never) - SA| by covariate set:\n")
print(chk %>% group_by(cov) %>%
        summarize(max_abs_diff=max(abs(`CS: never-treated` - `Sun & Abraham`), na.rm=TRUE)),
      row.names=FALSE)
