# =============================================================================
# plot_calendar_event.R  (Figure 2 of the follow-up pair)
#
# Show each estimator aggregated along its "natural" reporting dimension:
#   - CS  -> CALENDAR time:  ATT(t)  = did::aggte(type = "calendar")
#   - SA  -> EVENT time:     ATT(e)  = the interaction-weighted event-study
#                                       coefficients from fixest::sunab()
# Each panel plots the aggregated estimate with a 95% CI. This is the standard
# calendar-time / event-study presentation (estimate vs. time), a complement to
# the weight-vs-beta atom views.
#
# Data prep mirrors run_extraction.R (P&P column names, cohort recoding).
# =============================================================================

suppressMessages({library(did); library(fixest); library(dplyr); library(ggplot2); library(patchwork)})

DATA_DIR <- Sys.getenv("DIDREP_DATA", unset = "data")
FIG_DIR  <- Sys.getenv("DIDREP_FIG",  unset = "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

col_cs <- "#0072B2"; col_sa <- "#D55E00"

# ---- data prep (identical to the driver) -----------------------------------
dta <- read.csv(file.path(DATA_DIR, "dta.csv"))
dta$agency.num <- as.numeric(factor(dta$agency.id))
# irreversible treatment: collapse each unit to its first treatment year (2 units
# adopted+dropped; see run_extraction.R). Required by did >= 2.2.
dta$year.changed <- ave(dta$year.changed, dta$agency.id,
                        FUN = function(x) { v <- x[!is.na(x)]; if (length(v)) min(v) else NA_real_ })
dta$year.changed <- ifelse(is.na(dta$year.changed) & dta$no.req == 0, 0,
                    ifelse(is.na(dta$year.changed) & dta$no.req == 1, 1987,
                           dta$year.changed))

# ---- CS aggregated to CALENDAR time ----------------------------------------
m <- did::att_gt(yname = "any.fatalities", tname = "year", idname = "agency.num",
                 gname = "year.changed", control_group = "nevertreated",
                 clustervars = "agency.num", data = dta)
cal <- did::aggte(m, type = "calendar", na.rm = TRUE)
cs_df <- data.frame(time = cal$egt, estimate = cal$att.egt, se = cal$se.egt) %>%
  mutate(lo = estimate - 1.96 * se, hi = estimate + 1.96 * se)
cs_overall <- cal$overall.att

# ---- SA aggregated to EVENT time -------------------------------------------
d <- dta %>% filter(!is.na(any.fatalities)) %>% filter(year.changed == 0 | year.changed >= 2001)
d$coh <- ifelse(d$year.changed == 0, 10000, d$year.changed)
res <- fixest::feols(any.fatalities ~ sunab(coh, year) | agency.num + year,
                     data = d, cluster = ~agency.num)
 es <- summary(res)$coeftable            # interaction-weighted att(e)
cn <- rownames(es)
e  <- as.integer(vapply(regmatches(cn, regexec("::(-?[0-9]+)$", cn)),
                        function(z) z[2], character(1)))
sa_df <- data.frame(time = e, estimate = es[, 1], se = es[, 2]) %>%
  filter(!is.na(time)) %>%
  mutate(lo = estimate - 1.96 * se, hi = estimate + 1.96 * se) %>%
  arrange(time)
sa_overall <- fixest::coeftable(summary(res, agg = "att"))["ATT", 1]

base <- theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 8.7, color = "grey30"),
        panel.grid.minor = element_blank())

p_cs <- ggplot(cs_df, aes(time, estimate)) +
  geom_hline(yintercept = 0, color = "grey60", linetype = "dashed", linewidth = 0.35) +
  geom_hline(yintercept = cs_overall, color = col_cs, linetype = "dotted", linewidth = 0.5) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2, color = col_cs, alpha = 0.8) +
  geom_line(color = col_cs, linewidth = 0.4, alpha = 0.6) +
  geom_point(color = col_cs, size = 2.4) +
  annotate("text", x = max(cs_df$time), y = cs_overall, vjust = -0.6, hjust = 1,
           label = sprintf("overall ATT = %+.3f", cs_overall), color = col_cs,
           size = 3, fontface = "bold") +
  labs(title = "Callaway & Sant'Anna, aggregated to calendar time",
       subtitle = "ATT(t): average effect on P(fatal encounter) in each calendar year (95% CI). Dotted line = overall ATT.",
       x = "Calendar year", y = "ATT(t)") + base

p_sa <- ggplot(sa_df, aes(time, estimate)) +
  geom_hline(yintercept = 0, color = "grey60", linetype = "dashed", linewidth = 0.35) +
  geom_vline(xintercept = -0.5, color = "grey60", linetype = "dotted", linewidth = 0.4) +
  geom_hline(yintercept = sa_overall, color = col_sa, linetype = "dotted", linewidth = 0.5) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.3, color = col_sa, alpha = 0.8) +
  geom_line(color = col_sa, linewidth = 0.4, alpha = 0.6) +
  geom_point(color = col_sa, size = 2.4) +
  annotate("text", x = max(sa_df$time), y = sa_overall, vjust = -0.6, hjust = 1,
           label = sprintf("overall ATT = %+.3f", sa_overall), color = col_sa,
           size = 3, fontface = "bold") +
  labs(title = "Sun & Abraham, aggregated to event time",
       subtitle = "ATT(e): interaction-weighted effect by event time relative to treatment (95% CI). Vertical line = treatment.",
       x = "Event time (years relative to treatment)", y = "ATT(e)") + base

p <- p_cs / p_sa +
  plot_annotation(title = "Standard aggregations: CS by calendar time, SA by event time",
                  theme = theme(plot.title = element_text(face = "bold", size = 13.5)))

ggsave(file.path(FIG_DIR, "calendar_and_event_time.pdf"), p,
       width = 9, height = 9, device = cairo_pdf)
ggsave(file.path(FIG_DIR, "calendar_and_event_time.png"), p,
       width = 9, height = 9, dpi = 200)
cat("Wrote figures/calendar_and_event_time.{pdf,png}\n")
cat(sprintf("CS calendar overall = %+.4f | SA event overall = %+.4f\n", cs_overall, sa_overall))
