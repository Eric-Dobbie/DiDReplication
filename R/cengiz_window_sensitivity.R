# =============================================================================
# cengiz_window_sensitivity.R
#
# Sensitivity of the Cengiz stacked coefficient to the event-window half-width k.
# For each k in 2..10 we build a balanced +/-k stack (2k+1 rows per unit),
# never-treated clean controls, stack-specific unit+time FE, dropping the two
# reversible agencies, and report the pooled coefficient on `no.req`.
#
# CAVEAT built into the design: the outcome (any.fatalities) is observed only
# 2000-2020, so a balanced +/-k window needs g-k >= 2000 AND g+k <= 2020. As k
# grows the eligible-cohort set shrinks -- by k=9 only cohort 2009 qualifies and
# k=10 leaves none. So larger k trades more periods-per-cohort for fewer cohorts;
# the coefficient movement is a mix of both and is NOT a clean "long-run" trace.
# =============================================================================

suppressMessages({library(lfe); library(dplyr); library(ggplot2)})

DATA_DIR <- Sys.getenv("DIDREP_DATA", unset="data")
OUT_DIR  <- Sys.getenv("DIDREP_OUT",  unset="output")
FIG_DIR  <- Sys.getenv("DIDREP_FIG",  unset="figures")
FIRST_OUT <- 2000; LAST_OUT <- 2020
REVERSIBLE <- c("memphis tennessee","portsmouth new hampshire")

dta <- read.csv(file.path(DATA_DIR, "dta.csv"))
dta <- dta[!dta$agency.id %in% REVERSIBLE, ]
nochange   <- dta[dta$change.type == "No Change", ]
allcohorts <- sort(unique(dta$year.changed[!is.na(dta$year.changed) & dta$year.changed >= FIRST_OUT]))

fit_k <- function(k) {
  cohorts <- allcohorts[allcohorts - k >= FIRST_OUT & allcohorts + k <= LAST_OUT]
  parts <- list()
  for (g in cohorts) {
    tr <- dta[!is.na(dta$year.changed) & dta$year.changed == g &
                dta$year >= g - k & dta$year <= g + k, ]
    if (nrow(tr) == 0) next
    tr$cohort <- g; tr$treat <- 1L
    ct <- nochange[nochange$year >= g - k & nochange$year <= g + k, ]
    ct$cohort <- g; ct$treat <- 0L
    parts[[length(parts) + 1L]] <- rbind(tr, ct)
  }
  used <- unique(unlist(lapply(parts, function(p) p$cohort[p$treat==1])))
  if (length(parts) == 0)
    return(data.frame(k=k, rows_per_unit=2*k+1, n_cohorts=0, n_treated=0,
                      estimate=NA, se=NA))
  cz <- do.call(rbind, parts)
  cz$stackunit <- paste0(cz$agency.id, "__", cz$cohort)
  cz$stacktime <- paste0(cz$year, "__", cz$cohort)
  m <- felm(any.fatalities ~ no.req | stackunit + stacktime | 0 | agency.id, data = cz)
  data.frame(k=k, rows_per_unit=2*k+1, n_cohorts=length(used),
             n_treated=sum(cz$treat==1)/(2*k+1),
             estimate=unname(coef(m)["no.req"]), se=unname(m$se["no.req"]))
}

res <- do.call(rbind, lapply(2:10, fit_k))
res$lo <- res$estimate - 1.96*res$se
res$hi <- res$estimate + 1.96*res$se

# P&P full-panel plain reference
pp <- read.csv(file.path(DATA_DIR, "stacked_fatal.csv"))
b_pp <- unname(coef(felm(any.fatalities ~ no.req | agency.id + year.cohort | 0 | agency.id, data=pp))["no.req"])

cat("== Cengiz window sensitivity (never-treated controls) ==\n")
print(transform(res, estimate=round(estimate,4), se=round(se,4))[c("k","rows_per_unit","n_cohorts","n_treated","estimate","se")],
      row.names=FALSE)
cat(sprintf("\nreference: P&P full-panel stacked = %+.4f\n", b_pp))
write.csv(res, file.path(OUT_DIR, "cengiz_window_sensitivity.csv"), row.names=FALSE)

# ---- plot: coefficient vs window half-width --------------------------------
rp <- res[!is.na(res$estimate), ]
p <- ggplot(rp, aes(k, estimate)) +
  geom_hline(yintercept=0, linetype="dashed", color="grey60", linewidth=0.35) +
  geom_hline(yintercept=b_pp, linetype="dotted", color="#009E73", linewidth=0.7) +
  annotate("text", x=max(rp$k), y=b_pp, label=sprintf("P&P full-panel = %+.3f", b_pp),
           vjust=-0.6, hjust=1, color="#009E73", size=3, fontface="bold") +
  geom_ribbon(aes(ymin=lo, ymax=hi), fill="#CC79A7", alpha=0.18) +
  geom_line(color="#CC79A7", linewidth=0.7) +
  geom_point(color="#CC79A7", size=3) +
  geom_text(aes(label=n_cohorts), vjust=-1.1, size=2.9, color="grey30") +
  scale_x_continuous(breaks=2:10) +
  labs(title="Cengiz stacked coefficient vs. event-window half-width (+/-k)",
       subtitle="Pooled effect on P(fatal encounter), 95% CI band. Number above each point = eligible cohorts.\nWider windows require more observed pre/post years, so fewer cohorts qualify (k>=9: only cohort 2009).",
       x="Window half-width k (rows per unit = 2k+1)", y="Pooled stacked estimate") +
  theme_bw(base_size=12) +
  theme(plot.title=element_text(face="bold", size=12.5),
        plot.subtitle=element_text(size=8.5, color="grey30"),
        panel.grid.minor=element_blank())
ggsave(file.path(FIG_DIR, "cengiz_window_sensitivity.pdf"), p, width=9, height=5.5, device=cairo_pdf)
ggsave(file.path(FIG_DIR, "cengiz_window_sensitivity.png"), p, width=9, height=5.5, dpi=200)
cat("Wrote figures/cengiz_window_sensitivity.{pdf,png} and output/cengiz_window_sensitivity.csv\n")
