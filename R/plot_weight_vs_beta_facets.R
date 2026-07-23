# =============================================================================
# plot_weight_vs_beta_facets.R
#
# Small-multiples version of the weight-vs-beta figure: one panel per estimator,
# with a FREE weight (y) axis so the CS and SA weight structure -- crushed in the
# combined plot by the stacked cohort-2009 atom at 0.37 -- becomes legible.
# Shared x (component estimate) keeps the estimate axis comparable across panels.
# Each panel adds a vertical line + diamond at that estimator's overall
# aggregated ATT, and labels the highest-weight cohorts.
# =============================================================================

suppressMessages({library(ggplot2); library(dplyr); library(ggrepel)})

OUT_DIR <- Sys.getenv("DIDREP_OUT", unset = "output")
FIG_DIR <- Sys.getenv("DIDREP_FIG", unset = "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

atoms <- read.csv(file.path(OUT_DIR, "atoms_long.csv"))
chk   <- read.csv(file.path(OUT_DIR, "recombination_check.csv"))

lev  <- c("CS", "SA", "stacked")
labs <- c(CS = "Callaway & Sant'Anna  —  ATT(g,t) atoms",
          SA = "Sun & Abraham  —  CATT(g,e) atoms",
          stacked = "Stacked regression  —  cohort-stack atoms")
pal  <- c(CS = "#0072B2", SA = "#D55E00", stacked = "#009E73")

atoms <- atoms %>% filter(weight > 0) %>%
  mutate(estimator = factor(estimator, lev))

# overall aggregated ATT per estimator, placed at the panel's mean atom weight
agg <- atoms %>% group_by(estimator) %>%
  summarize(y = mean(weight), ymax = max(weight), .groups = "drop")
pooled <- chk %>% filter(estimator %in% lev) %>%
  transmute(estimator = factor(estimator, lev), x = package_pooled)
agg <- left_join(agg, pooled, by = "estimator") %>%
  mutate(lab = sprintf("aggregate ATT = %+.3f", x))

# label the highest-weight cohorts in each panel (one label per cohort, at its
# max-weight atom); top 3 for CS/SA, top 5 for stacked, + the stacked 2002 pivot
top_lab <- atoms %>%
  group_by(estimator, group) %>%
  slice_max(weight, n = 1, with_ties = FALSE) %>%
  group_by(estimator) %>%
  mutate(rk = rank(-weight, ties.method = "first")) %>%
  filter(rk <= ifelse(estimator == "stacked", 5, 3) |
           (estimator == "stacked" & estimate > 0)) %>%   # keep 2002 pivot
  ungroup() %>%
  mutate(clab = ifelse(estimator == "stacked" & estimate > 0,
                       paste0(group, "\n(no controls)"), as.character(group)))

p <- ggplot(atoms, aes(x = estimate, y = weight, color = estimator)) +
  geom_vline(xintercept = 0, linewidth = 0.35, color = "grey65", linetype = "dashed") +
  # aggregate ATT: vertical line + point on the mean-weight level
  geom_vline(data = agg, aes(xintercept = x, color = estimator),
             linewidth = 0.6, linetype = "solid", alpha = 0.5, show.legend = FALSE) +
  geom_point(size = 2.1, alpha = 0.6, stroke = 0.2) +
  geom_point(data = agg, aes(x = x, y = y, fill = estimator),
             shape = 23, size = 4.6, color = "black", stroke = 0.7,
             inherit.aes = FALSE) +
  geom_text(data = agg, aes(label = lab, color = estimator),
            x = Inf, y = Inf, hjust = 1.03, vjust = 1.6, size = 3.1,
            fontface = "bold", inherit.aes = FALSE, show.legend = FALSE) +
  geom_text_repel(data = top_lab, aes(label = clab),
                  size = 2.9, fontface = "bold", lineheight = 0.85,
                  min.segment.length = 0, box.padding = 0.35, seed = 3,
                  segment.size = 0.25, show.legend = FALSE) +
  facet_wrap(~estimator, ncol = 1, scales = "free_y",
             labeller = labeller(estimator = labs)) +
  scale_color_manual(values = pal, guide = "none", aesthetics = c("color", "fill")) +
  scale_x_continuous(breaks = scales::pretty_breaks(9)) +
  labs(
    title    = "Component estimate vs. aggregation weight, by estimator",
    subtitle = "Outcome: P(fatal civilian encounter). Each panel has its own weight scale. CS and SA coincide exactly here\n(balanced panel), so their panels are identical. Vertical line + diamond = overall ATT; labels = highest-weight cohorts.",
    x = "Component estimate (effect on probability of a fatal encounter)",
    y = "Aggregation weight"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold", size = 13),
    plot.subtitle   = element_text(size = 8.7, color = "grey30"),
    strip.text      = element_text(face = "bold", size = 10.5),
    strip.background = element_rect(fill = "grey93"),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(FIG_DIR, "weight_vs_beta_smallmultiples.pdf"), p,
       width = 8.5, height = 9, device = cairo_pdf)
ggsave(file.path(FIG_DIR, "weight_vs_beta_smallmultiples.png"), p,
       width = 8.5, height = 9, dpi = 200)
cat("Wrote figures/weight_vs_beta_smallmultiples.{pdf,png}\n")
