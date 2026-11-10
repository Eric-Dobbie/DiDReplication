# =============================================================================
# design_sweep_figureAB.R  --  deterministic; no outcomes.
# Figure A: each cohort's aggregation weight vs its pre-period count, one line
#           per scheme; panel 1 equal cohort sizes, panel 2 unequal (so size vs
#           timing contributions separate).
# Figure B: across swept configs, cor(w_simple, w_stacked) and L1 distance vs
#           cohort placement (centroid). Reports the weighting-divergence bound
#           L1*(b-a)/2 for an assumed cohort-effect range [a,b].
# =============================================================================
suppressMessages({library(ggplot2); library(tidyr)})
.sourced_quiet <- TRUE; source("R/design_sweep_core.R")
OUT <- Sys.getenv("DIDREP_OUT", unset="output"); FIG <- Sys.getenv("DIDREP_FIG", unset="figures")
N_C <- 40

# ---- Figure A ---------------------------------------------------------------
mkA <- function(T, G, sizetype) {
  gp <- unique(round(seq(2, T, length.out=G)))
  ng <- if (sizetype=="equal") rep(3L,length(gp)) else as.integer(round(seq(8,1,length.out=length(gp))))
  dw <- design_weights(T, gp, ng, N_C)
  long <- pivot_longer(dw[,c("pre","n_g","w_simple","w_group","w_dynamic","w_calendar","w_stacked")],
    cols=starts_with("w_"), names_to="scheme", values_to="weight")
  long$scheme <- sub("^w_","",long$scheme); long$panel <- sprintf("%s sizes", sizetype); long
}
dA <- rbind(mkA(30,10,"equal"), mkA(30,10,"unequal"))
dA$scheme <- factor(dA$scheme, levels=c("simple","group","dynamic","calendar","stacked"))
pal <- c(simple="#0072B2", group="#56B4E9", dynamic="#E69F00", calendar="#CC79A7", stacked="#009E73")
gA <- ggplot(dA, aes(pre, weight, color=scheme)) +
  geom_line(linewidth=0.8) + geom_point(size=1.8) +
  facet_wrap(~panel) + scale_color_manual(values=pal, name=NULL) +
  labs(x="cohort pre-period count (pre_g = gp - 1)", y="aggregation weight on cohort") +
  theme_minimal(base_size=12) + theme(legend.position="top", panel.grid.minor=element_blank())
ggsave(file.path(FIG,"design_weight_vs_pre.pdf"), gA, width=10, height=4.5, device=cairo_pdf)
ggsave(file.path(FIG,"design_weight_vs_pre.png"), gA, width=10, height=4.5, dpi=200)
cat("Figure A shapes (equal sizes): is CS simple linear-decreasing in pre? cor(pre, w_simple) =",
    round(cor(mkA(30,10,"equal")$pre[1:10], design_weights(30,unique(round(seq(2,30,length.out=10))),rep(3L,10),N_C)$w_simple),3),"\n")

# ---- Figure B ---------------------------------------------------------------
place <- function(T,G,type){ s<-2:T; if(type=="early") head(s,G) else if(type=="late") tail(s,G)
  else if(type=="spread") unique(round(seq(2,T,length.out=G))) else sort(sample(s,G)) }
rows <- list(); set.seed(7)
for (T in c(21,30)) for (G in c(4,6,8,10)) for (pl in c("early","late","spread","random"))
  for (sz in c("equal","skewed")) for (rep in if(pl=="random") 1:40 else 1) {
    gp <- place(T,G,pl); if (length(gp)<3) next
    ng <- if (sz=="equal") rep(3L,length(gp)) else { v<-rep(1L,length(gp)); v[1]<-8L; if(length(gp)>1) v[2]<-5L; v }
    dw <- design_weights(T, gp, ng, N_C)
    ok <- !is.na(dw$w_stacked)
    if (sum(ok)<3) next
    cc <- cor(dw$w_simple[ok], dw$w_stacked[ok]); L1 <- sum(abs(dw$w_simple[ok]-dw$w_stacked[ok]))
    rows[[length(rows)+1L]] <- data.frame(T=T,G=length(gp),placement=pl,size=sz,
      centroid=mean(gp)/T, cor=cc, L1=L1) }
B <- do.call(rbind, rows); write.csv(B, file.path(OUT,"design_sweep_corL1.csv"), row.names=FALSE)
gB <- ggplot(B, aes(centroid, cor, color=placement)) +
  geom_hline(yintercept=0, linetype="dashed", color="grey60") +
  geom_point(aes(size=L1), alpha=0.6) +
  scale_color_manual(values=c(early="#0072B2",late="#D55E00",spread="#009E73",random="grey55"), name=NULL) +
  scale_size_continuous(name="L1") +
  labs(x="cohort placement (mean adoption period / T)", y="cor(w_CS_simple, w_stacked)") +
  theme_minimal(base_size=12) + theme(legend.position="top", panel.grid.minor=element_blank())
ggsave(file.path(FIG,"design_weight_agreement.pdf"), gB, width=9, height=5, device=cairo_pdf)
ggsave(file.path(FIG,"design_weight_agreement.png"), gB, width=9, height=5, dpi=200)

cat(sprintf("\ncor(w_simple,w_stacked): range [%.2f, %.2f]; configs with cor<0: %d of %d\n",
    min(B$cor), max(B$cor), sum(B$cor<0), nrow(B)))
cat("negative-correlation configs are concentrated in placement:\n"); print(table(B$placement[B$cor<0]))
cat(sprintf("L1 range [%.2f, %.2f].  Weighting-divergence bound for cohort effects in [a,b]: L1*(b-a)/2.\n", min(B$L1), max(B$L1)))
cat(sprintf("  e.g. [a,b]=[-0.5,0.5] (b-a=1): max |CS_simple - stacked| attributable to weighting alone = L1*0.5, up to %.3f across the sweep.\n", max(B$L1)*0.5))
# P&P point
yrs <- c(2001,2002,2003,2004,2008,2009,2012,2013,2014,2017,2018,2019,2020)
ppw <- design_weights(21, yrs-2000+1, c(1,5,1,2,1,7,1,3,1,3,3,1,1), N_C); ok<-!is.na(ppw$w_stacked)
cat(sprintf("P&P: cor(w_simple,w_stacked)=%.3f, L1=%.3f, centroid=%.3f\n",
    cor(ppw$w_simple[ok],ppw$w_stacked[ok]), sum(abs(ppw$w_simple[ok]-ppw$w_stacked[ok])), mean(yrs-2000+1)/21))
cat("Wrote figures/design_weight_vs_pre.*, design_weight_agreement.*, output/design_sweep_corL1.csv\n")
