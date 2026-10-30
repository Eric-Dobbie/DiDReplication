# =============================================================================
# cengiz_asym_window.R
#
# Cengiz stacked design with an ASYMMETRIC event window: the pre-treatment side
# is fixed at 4 periods, and the post-treatment side p grows from 4 upward. Each
# sub-experiment keeps years [g-4, g+p] (4 pre + treatment year + p post =
# p+5 rows per unit); the 741 No-Change clean controls (P&P's pool) are windowed
# the same way. The two reversible agencies are dropped.
#
# Eligibility: a cohort g is included only if its whole window is observed --
#   pre observable:  g - 4 >= 2000   (so g >= 2004; 2000-2003 never qualify)
#   post observable: g + p <= 2020
# As p grows, later cohorts fall out; we sweep p up until only ONE cohort
# remains eligible (2004, whose full observable post reaches 2020 at p=16).
#
# Stack-specific unit + time FE (agency x cohort + year x cohort); SE clustered
# on agency.id. Reports the pooled coefficient on `no.req` vs the post-window p.
# =============================================================================
suppressMessages({library(lfe); library(ggplot2)})
DATA_DIR <- Sys.getenv("DIDREP_DATA", unset = "data")
OUT_DIR  <- Sys.getenv("DIDREP_OUT",  unset = "output")
FIG_DIR  <- Sys.getenv("DIDREP_FIG",  unset = "figures")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
FIRST_OUT <- 2000; LAST_OUT <- 2020; PRE <- 4
REVERSIBLE <- c("memphis tennessee", "portsmouth new hampshire")

dta <- read.csv(file.path(DATA_DIR, "dta.csv"))
dta <- dta[!dta$agency.id %in% REVERSIBLE, ]
nochange   <- dta[dta$change.type == "No Change", ]              # 741 clean controls
allcohorts <- sort(unique(dta$year.changed[!is.na(dta$year.changed) & dta$year.changed >= FIRST_OUT]))

fit_post <- function(p) {
  # fixed 4 pre-periods, p post-periods
  cohorts <- allcohorts[allcohorts - PRE >= FIRST_OUT & allcohorts + p <= LAST_OUT]
  parts <- list()
  for (g in cohorts) {
    tr <- dta[!is.na(dta$year.changed) & dta$year.changed == g &
                dta$year >= g - PRE & dta$year <= g + p, ]
    if (nrow(tr) == 0) next
    tr$cohort <- g; tr$treat <- 1L
    ct <- nochange[nochange$year >= g - PRE & nochange$year <= g + p, ]
    ct$cohort <- g; ct$treat <- 0L
    parts[[length(parts) + 1L]] <- rbind(tr, ct)
  }
  used <- unique(unlist(lapply(parts, function(x) x$cohort[x$treat == 1])))
  if (length(parts) == 0)
    return(data.frame(post = p, rows_per_unit = PRE + 1 + p, n_cohorts = 0,
                      n_treated = 0, estimate = NA, se = NA))
  cz <- do.call(rbind, parts)
  cz$stackunit <- paste0(cz$agency.id, "__", cz$cohort)
  cz$stacktime <- paste0(cz$year, "__", cz$cohort)
  m <- felm(any.fatalities ~ no.req | stackunit + stacktime | 0 | agency.id, data = cz)
  data.frame(post = p, rows_per_unit = PRE + 1 + p, n_cohorts = length(used),
             n_treated = sum(cz$treat == 1) / (PRE + 1 + p),
             estimate = unname(coef(m)["no.req"]), se = unname(m$se["no.req"]))
}

# sweep post from 4 upward while >= 1 cohort qualifies
res <- list(); p <- PRE
repeat {
  r <- fit_post(p); res[[length(res) + 1L]] <- r
  if (r$n_cohorts <= 1) break        # stop once only one cohort remains
  p <- p + 1
}
res <- do.call(rbind, res)
res$lo <- res$estimate - 1.96 * res$se
res$hi <- res$estimate + 1.96 * res$se

# P&P full-panel plain reference (shared agency.id + year.cohort FE)
pp   <- read.csv(file.path(DATA_DIR, "stacked_fatal.csv"))
b_pp <- unname(coef(felm(any.fatalities ~ no.req | agency.id + year.cohort | 0 | agency.id,
                         data = pp))["no.req"])

cat("== Cengiz asymmetric window: 4 pre, post 4..max (741 controls) ==\n")
print(transform(res, estimate = round(estimate, 4), se = round(se, 4))[
        c("post", "rows_per_unit", "n_cohorts", "n_treated", "estimate", "se")], row.names = FALSE)
cat(sprintf("\nreference: P&P full-panel stacked = %+.4f\n", b_pp))
write.csv(res, file.path(OUT_DIR, "cengiz_asym_window.csv"), row.names = FALSE)

# ---- plot: coefficient vs post-window length (minimal style) ----------------
rp <- res[!is.na(res$estimate), ]
p_fig <- ggplot(rp, aes(post, estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.35) +
  geom_hline(yintercept = b_pp, linetype = "dotted", color = "#009E73", linewidth = 0.7) +
  annotate("text", x = max(rp$post), y = b_pp, label = sprintf("P&P full-panel = %+.3f", b_pp),
           vjust = -0.6, hjust = 1, color = "#009E73", size = 3, fontface = "bold") +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#0072B2", alpha = 0.15) +
  geom_line(color = "#0072B2", linewidth = 0.7) +
  geom_point(color = "#0072B2", size = 3) +
  geom_text(aes(label = n_cohorts), vjust = -1.1, size = 2.9, color = "grey30") +
  scale_x_continuous(breaks = min(rp$post):max(rp$post)) +
  labs(x = "Post-treatment periods (pre-periods fixed at 4; rows per unit = 5 + post)",
       y = "Pooled stacked estimate") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(FIG_DIR, "cengiz_asym_window.pdf"), p_fig, width = 9, height = 5.5, device = cairo_pdf)
ggsave(file.path(FIG_DIR, "cengiz_asym_window.png"), p_fig, width = 9, height = 5.5, dpi = 200)
cat("Wrote figures/cengiz_asym_window.{pdf,png} and output/cengiz_asym_window.csv\n")
