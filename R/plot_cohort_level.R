# =============================================================================
# plot_cohort_level.R  (Figure 1 of the follow-up pair)
#
# Aggregate every estimator to a COMMON level -- one point per treated cohort --
# so CS, SA, and stacked are directly comparable. Stacked is already cohort-
# level; CS (ATT(g,t)) and SA (CATT(g,e)) are rolled up within each cohort by
# their own atom weights:
#     estimate_g = sum_i(w_i * beta_i) / sum_i(w_i)
#     weight_g   = sum_i(w_i)
# an exact regrouping, so sum_g(weight_g * estimate_g) still equals each
# estimator's overall ATT.
#
# NOTE: for this balanced panel with a never-treated control group, CS and SA
# are algebraically equivalent, so their cohort points coincide to machine
# precision -- drawn with distinct shapes so the overlap is visible. The real
# contrast is with the stacked estimator, which uses a different (and larger)
# clean-control pool.
# =============================================================================

suppressMessages({library(ggplot2); library(dplyr); library(ggrepel)})

OUT_DIR <- Sys.getenv("DIDREP_OUT", unset = "output")
FIG_DIR <- Sys.getenv("DIDREP_FIG", unset = "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

atoms <- read.csv(file.path(OUT_DIR, "atoms_long.csv"))
chk   <- read.csv(file.path(OUT_DIR, "recombination_check.csv"))

lev  <- c("CS", "SA", "stacked")
labs <- c(CS = "Callaway & Sant'Anna", SA = "Sun & Abraham", stacked = "Stacked")
pal  <- c(CS = "#0072B2", SA = "#D55E00", stacked = "#009E73")
shp  <- c(CS = 16, SA = 5, stacked = 15)          # filled circle / open diamond / filled square
sz   <- c(CS = 2.8, SA = 4.4, stacked = 2.8)      # SA larger so its ring encloses the CS dot

coh <- atoms %>%
  filter(weight > 0) %>%
  group_by(estimator, group) %>%
  summarize(estimate = sum(weight * estimate) / sum(weight),
            weight   = sum(weight), .groups = "drop") %>%
  mutate(estimator = factor(estimator, lev))

# overall ATTs: CS and SA coincide, stacked separate -> vertical reference lines
overall <- chk %>% filter(estimator %in% lev) %>%
  transmute(estimator = factor(estimator, lev), x = package_pooled)

# label key cohorts on the CS(=SA) cloud and on the stacked column
key <- c(2002, 2004, 2009, 2013, 2017)
coh_lab <- coh %>% filter(group %in% key & estimator %in% c("CS", "stacked"))

p <- ggplot() +
  geom_vline(xintercept = 0, linewidth = 0.35, color = "grey60", linetype = "dashed") +
  geom_vline(data = overall, aes(xintercept = x, color = estimator),
             linewidth = 0.6, linetype = "dotted", alpha = 0.7, show.legend = FALSE) +
  geom_point(data = coh, aes(x = estimate, y = weight, color = estimator,
                             shape = estimator, size = estimator), stroke = 0.9) +
  geom_text_repel(data = coh_lab,
                  aes(x = estimate, y = weight, label = group, color = estimator),
                  size = 2.9, fontface = "bold", min.segment.length = 0,
                  box.padding = 0.4, segment.size = 0.2, max.overlaps = 20,
                  seed = 5, show.legend = FALSE) +
  scale_color_manual(values = pal, labels = labs, name = "Estimator") +
  scale_shape_manual(values = shp, labels = labs, name = "Estimator") +
  scale_size_manual(values = sz, labels = labs, name = "Estimator", guide = "none") +
  scale_x_continuous(breaks = scales::pretty_breaks(9)) +
  guides(color = guide_legend(override.aes = list(size = c(3, 4.4, 3)))) +
  labs(
    title    = "Cohort-level view: one point per treated cohort",
    subtitle = "Outcome: P(fatal civilian encounter). CS/SA atoms rolled up to cohort by their weights; stacked is already cohort-level.\nCS and SA coincide exactly here (balanced panel + never-treated controls) -- their points overlap. Dotted lines = overall ATT.",
    x = "Cohort-level estimate (effect on probability of a fatal encounter)",
    y = "Cohort aggregation weight"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 8.5, color = "grey30"),
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(FIG_DIR, "cohort_level_estimates.pdf"), p,
       width = 9, height = 6.5, device = cairo_pdf)
ggsave(file.path(FIG_DIR, "cohort_level_estimates.png"), p,
       width = 9, height = 6.5, dpi = 200)
cat("Wrote figures/cohort_level_estimates.{pdf,png}\n")
