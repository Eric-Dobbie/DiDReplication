# =============================================================================
# plot_weight_vs_beta.R
#
# Weight-vs-beta decomposition figure for the fatal-encounters DiD estimators.
# Reads the atoms + recombination outputs written by run_extraction.R and draws:
#   - each estimator's component atoms as points  (x = estimate, y = weight)
#   - one large diamond per estimator at its overall/aggregated estimate
#     (x = pooled ATT, y = mean weight of that estimator's weighted atoms)
#   - a shaded triangle = convex hull of the three aggregate points
#   - color encodes estimator throughout
#
# Colorblind-safe categorical palette (validated): CS #0072B2, SA #D55E00,
# Stacked #009E73.
# =============================================================================

suppressMessages({library(ggplot2); library(dplyr); library(ggrepel)})

OUT_DIR <- Sys.getenv("DIDREP_OUT",  unset = "output")
FIG_DIR <- Sys.getenv("DIDREP_FIG",  unset = "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

atoms <- read.csv(file.path(OUT_DIR, "atoms_long.csv"))
chk   <- read.csv(file.path(OUT_DIR, "recombination_check.csv"))

lev  <- c("CS", "SA", "stacked")
labs <- c(CS = "Callaway & Sant'Anna", SA = "Sun & Abraham", stacked = "Stacked")
pal  <- c(CS = "#0072B2", SA = "#D55E00", stacked = "#009E73")

# atoms that actually enter the aggregate (positive weight)
atoms <- atoms %>% filter(weight > 0) %>% mutate(estimator = factor(estimator, lev))

# aggregate points: x = package-reported pooled ATT, y = mean weighted-atom weight
agg <- atoms %>% group_by(estimator) %>%
  summarize(y = mean(weight), .groups = "drop")
pooled <- chk %>% filter(estimator %in% lev) %>%
  transmute(estimator = factor(estimator, lev), x = package_pooled)
agg <- left_join(agg, pooled, by = "estimator") %>%
  mutate(lab = sprintf("%s\nATT = %+.3f", labs[as.character(estimator)], x))

# convex hull of the three aggregate points (a triangle for 3 non-collinear pts)
hull <- agg[chull(agg$x, agg$y), c("x", "y")]

# label anchors so ggrepel pushes the CS/SA labels (near-coincident) apart
agg$lx <- c(CS = -0.35, SA = 0.45, stacked = -0.45)[as.character(agg$estimator)]
agg$ly <- c(CS = 0.045, SA = 0.045, stacked = 0.16)[as.character(agg$estimator)]

p <- ggplot() +
  # shaded convex hull between the aggregated estimates
  geom_polygon(data = hull, aes(x = x, y = y),
               fill = "grey40", alpha = 0.25, color = "grey35",
               linewidth = 0.5) +
  # null-effect reference
  geom_vline(xintercept = 0, linewidth = 0.35, color = "grey60", linetype = "dashed") +
  # component atoms
  geom_point(data = atoms, aes(x = estimate, y = weight, color = estimator),
             size = 2.1, alpha = 0.6, stroke = 0.2) +
  # aggregate markers
  geom_point(data = agg, aes(x = x, y = y, fill = estimator),
             shape = 23, size = 5.4, color = "black", stroke = 0.8) +
  # aggregate labels with leader lines
  geom_text_repel(data = agg,
                  aes(x = x, y = y, label = lab, color = estimator),
                  nudge_x = agg$lx - agg$x, nudge_y = agg$ly - agg$y,
                  size = 3.2, fontface = "bold", lineheight = 0.9,
                  segment.color = "grey45", segment.size = 0.3,
                  min.segment.length = 0, box.padding = 0.4,
                  seed = 1, show.legend = FALSE) +
  scale_color_manual(values = pal, labels = labs, name = "Estimator",
                     aesthetics = c("color", "fill")) +
  scale_x_continuous(breaks = scales::pretty_breaks(7)) +
  labs(
    title    = "Weight vs. component estimate across three DiD estimators",
    subtitle = "Outcome: P(fatal civilian encounter). Points = component atoms; diamonds = overall aggregated ATT;\nshaded region = convex hull of the three aggregates.",
    x = "Component estimate (effect on probability of a fatal encounter)",
    y = "Aggregation weight"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9, color = "grey30"),
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(FIG_DIR, "weight_vs_beta_decomposition.pdf"), p,
       width = 8.5, height = 6, device = cairo_pdf)
ggsave(file.path(FIG_DIR, "weight_vs_beta_decomposition.png"), p,
       width = 8.5, height = 6, dpi = 200)
cat("Wrote figures/weight_vs_beta_decomposition.{pdf,png}\n")
cat("Aggregate points:\n"); print(as.data.frame(agg[, c("estimator","x","y")]), row.names = FALSE)
