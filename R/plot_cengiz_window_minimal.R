# =============================================================================
# plot_cengiz_window_minimal.R
#
# Minimal reproduction of the Cengiz window-sensitivity figure: the pooled
# stacked coefficient vs the event-window half-width k, with its 95% CI band and
# the P&P full-panel reference line. No plot title, subtitle, or caption;
# theme_minimal(). Reads output/cengiz_window_sensitivity.csv (from
# cengiz_window_sensitivity.R); recomputes only the P&P reference line.
# =============================================================================
suppressMessages({library(lfe); library(ggplot2)})
DATA_DIR <- Sys.getenv("DIDREP_DATA", unset = "data")
OUT_DIR  <- Sys.getenv("DIDREP_OUT",  unset = "output")
FIG_DIR  <- Sys.getenv("DIDREP_FIG",  unset = "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

res <- read.csv(file.path(OUT_DIR, "cengiz_window_sensitivity.csv"))
rp  <- res[!is.na(res$estimate), ]

# P&P full-panel plain reference (shared agency.id + year.cohort FE)
pp   <- read.csv(file.path(DATA_DIR, "stacked_fatal.csv"))
b_pp <- unname(coef(felm(any.fatalities ~ no.req | agency.id + year.cohort | 0 | agency.id,
                         data = pp))["no.req"])

p <- ggplot(rp, aes(k, estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.35) +
  geom_hline(yintercept = b_pp, linetype = "dotted", color = "#009E73", linewidth = 0.7) +
  annotate("text", x = max(rp$k), y = b_pp, label = sprintf("P&P full-panel = %+.3f", b_pp),
           vjust = -0.6, hjust = 1, color = "#009E73", size = 3, fontface = "bold") +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#CC79A7", alpha = 0.18) +
  geom_line(color = "#CC79A7", linewidth = 0.7) +
  geom_point(color = "#CC79A7", size = 3) +
  geom_text(aes(label = n_cohorts), vjust = -1.1, size = 2.9, color = "grey30") +
  scale_x_continuous(breaks = 2:10) +
  labs(x = "Window half-width k (rows per unit = 2k+1)", y = "Pooled stacked estimate") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(FIG_DIR, "cengiz_window_sensitivity_minimal.pdf"), p,
       width = 9, height = 5.5, device = cairo_pdf)
ggsave(file.path(FIG_DIR, "cengiz_window_sensitivity_minimal.png"), p,
       width = 9, height = 5.5, dpi = 200)
cat("Wrote figures/cengiz_window_sensitivity_minimal.{pdf,png}\n")
