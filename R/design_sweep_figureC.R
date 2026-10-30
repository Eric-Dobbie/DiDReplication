# =============================================================================
# design_sweep_figureC.R  --  weight on UNCHECKABLE cohorts across the design space
#
# Deterministic enumeration (no outcomes, no DGP). At each design configuration
# we compute, per aggregation scheme, the share of total weight landing on
# cohorts whose joint pre-trend test is uncheckable (k/n_g > 1, from B1/B2). The
# random-placement arm samples DESIGNS (adoption-period placements), not data.
# Locates P&P in the distribution. Saves the full swept grid to CSV.
# =============================================================================
suppressMessages({library(ggplot2)})
.sourced_quiet <- TRUE
source("R/design_sweep_core.R")
OUT <- Sys.getenv("DIDREP_OUT", unset="output"); FIG <- Sys.getenv("DIDREP_FIG", unset="figures")
THRESH <- 1   # k/n_g > 1 is uncheckable (B1/B2)
N_C <- 40     # control pool common to all stacks (~ P&P's 41)

# ---- design generators ------------------------------------------------------
place <- function(T, G, type) {
  slots <- 2:T
  if (type=="early")  head(slots, G)
  else if (type=="late")   tail(slots, G)
  else if (type=="spread") unique(round(seq(2, T, length.out=G)))
  else if (type=="random") sort(sample(slots, G))
  else stop(type)
}
sizes <- function(G, type) {
  if (type=="equal")  rep(3L, G)
  else if (type=="skewed") { s <- rep(1L, G); s[1] <- 8L; if (G>=2) s[2] <- 5L; if (G>=3) s[3] <- 3L; s }  # a few big, rest single
  else stop(type)
}
uncheck_share <- function(dw) {
  unc <- !is.na(dw$k_over_n) & dw$k_over_n > THRESH
  sapply(c("w_simple","w_group","w_dynamic","w_calendar","w_stacked"), function(w){
    ww <- dw[[w]]; tot <- sum(ww, na.rm=TRUE); if (tot==0) NA_real_ else sum(ww[unc], na.rm=TRUE)/tot })
}

rows <- list(); add <- function(T,G,pl,sz,arm){
  gp <- place(T,G,pl); if (length(gp)<G) return(invisible())
  ng <- sizes(length(gp), sz); dw <- design_weights(T, gp, ng, N_C)
  us <- uncheck_share(dw)
  rows[[length(rows)+1L]] <<- data.frame(T=T, G=length(gp), placement=pl, size=sz, arm=arm,
    n_uncheck=sum(dw$k_over_n>THRESH), n_unident_stacked=sum(dw$id_stacked!="ok"),
    simple=us["w_simple"], group=us["w_group"], dynamic=us["w_dynamic"],
    calendar=us["w_calendar"], stacked=us["w_stacked"], row.names=NULL)
}
# deterministic grid
for (T in c(15,21,30)) for (G in c(3,5,8,11)) for (pl in c("early","late","spread"))
  for (sz in c("equal","skewed")) if (G <= T-1) add(T,G,pl,sz,"grid")
# random-placement arm (sampling designs, not data)
set.seed(11)
for (T in c(21)) for (G in c(5,8,11)) for (sz in c("equal","skewed")) for (r in 1:300) add(T,G,"random",sz,"random")
grid <- do.call(rbind, rows)
write.csv(grid, file.path(OUT,"design_sweep_grid.csv"), row.names=FALSE)
cat("swept configurations:", nrow(grid), "\n")

# ---- P&P position -----------------------------------------------------------
yrs <- c(2001,2002,2003,2004,2008,2009,2012,2013,2014,2017,2018,2019,2020)
pp <- design_weights(21, yrs-2000+1, c(1,5,1,2,1,7,1,3,1,3,3,1,1), N_C)
pp_us <- uncheck_share(pp)
cat("\nP&P uncheckable-weight share by scheme:\n"); print(round(pp_us,3))
cat("P&P cohorts with k/n_g>1 (uncheckable):", paste((yrs)[pp$k_over_n>THRESH], collapse=", "), "\n")

# ---- distribution figure ----------------------------------------------------
library(tidyr)
long <- pivot_longer(grid, cols=c(simple,group,dynamic,calendar,stacked),
                     names_to="scheme", values_to="share")
long$scheme <- factor(long$scheme, levels=c("simple","group","dynamic","calendar","stacked"))
ppl <- data.frame(
  scheme = factor(c("simple","group","dynamic","calendar","stacked"),
                  levels=c("simple","group","dynamic","calendar","stacked")),
  share  = as.numeric(pp_us[c("w_simple","w_group","w_dynamic","w_calendar","w_stacked")]))
g <- ggplot(long[!is.na(long$share),], aes(share)) +
  geom_histogram(breaks=seq(0,1,0.05), fill="#0072B2", color="white") +
  geom_vline(data=ppl, aes(xintercept=share), color="#D55E00", linewidth=0.9) +
  geom_text(data=ppl, aes(x=share, y=Inf, label=sprintf("P&P %.2f", share)),
            color="#D55E00", vjust=1.4, hjust=1.1, size=2.8) +
  facet_wrap(~scheme, ncol=5) +
  labs(x=sprintf("share of aggregation weight on uncheckable cohorts (k/n_g > %g)", THRESH), y="design configurations") +
  theme_minimal(base_size=11) + theme(panel.grid.minor=element_blank())
ggsave(file.path(FIG,"design_uncheckable_share.pdf"), g, width=12, height=3.6, device=cairo_pdf)
ggsave(file.path(FIG,"design_uncheckable_share.png"), g, width=12, height=3.6, dpi=200)

# ---- named-configuration table ----------------------------------------------
named <- grid[grid$T==21 & grid$G==11 & grid$arm=="grid" & grid$size=="skewed",
              c("placement","simple","group","dynamic","calendar","stacked")]
named <- rbind(named, data.frame(placement="P&P (actual)", simple=pp_us["w_simple"], group=pp_us["w_group"],
  dynamic=pp_us["w_dynamic"], calendar=pp_us["w_calendar"], stacked=pp_us["w_stacked"]))
cat("\nuncheckable share at named configs (T=21, G=11, skewed sizes):\n")
print(transform(named, simple=round(simple,3), group=round(group,3), dynamic=round(dynamic,3),
                calendar=round(calendar,3), stacked=round(stacked,3)), row.names=FALSE)
write.csv(named, file.path(OUT,"design_uncheckable_named.csv"), row.names=FALSE)
cat("\nWrote design_sweep_grid.csv, design_uncheckable_named.csv, figures/design_uncheckable_share.*\n")
