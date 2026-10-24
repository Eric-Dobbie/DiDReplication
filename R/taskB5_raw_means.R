# =============================================================================
# taskB5_raw_means.R
# Raw outcome means (no adjustment, no FE, no covariates): treated vs the 41
# never-treated controls, for cohorts 2002 and 2009 (62% of CS weight). If the
# series visibly diverge before treatment the pre-trend rejection is real; if
# they track each other and the joint test returns p~0, the test is broken.
# =============================================================================
suppressMessages({library(ggplot2)})
DATA_DIR <- Sys.getenv("DIDREP_DATA", unset="data"); OUT_DIR <- Sys.getenv("DIDREP_OUT", unset="output")
FIG_DIR <- Sys.getenv("DIDREP_FIG", unset="figures")
dta <- read.csv(file.path(DATA_DIR,"dta.csv"))
dta$yc <- ave(dta$year.changed, dta$agency.id, FUN=function(x){v<-x[!is.na(x)]; if(length(v)) min(v) else NA_real_})
never <- unique(dta$agency.id[is.na(dta$yc) & dta$no.req==0])

mk <- function(g) {
  tr <- dta[dta$yc==g & !is.na(dta$any.fatalities),]
  nt <- dta[dta$agency.id %in% never & !is.na(dta$any.fatalities),]
  rbind(
    aggregate(any.fatalities~year, tr, mean) |> transform(grp=sprintf("treated (n=%d)", length(unique(tr$agency.id)))),
    aggregate(any.fatalities~year, nt, mean) |> transform(grp=sprintf("never-treated (n=%d)", length(never))))|>
    transform(cohort=sprintf("Cohort %d (adopts %d)", g, g))
}
d <- rbind(cbind(mk(2002), gg=2002), cbind(mk(2009), gg=2009))
write.csv(d, file.path(OUT_DIR,"taskB5_raw_means.csv"), row.names=FALSE)

vl <- data.frame(gg=c(2002,2009), cohort=c("Cohort 2002 (adopts 2002)","Cohort 2009 (adopts 2009)"), x=c(2002,2009))
g <- ggplot(d, aes(year, any.fatalities, color=grp)) +
  geom_vline(data=vl, aes(xintercept=x), linetype="dotted", color="grey50") +
  geom_line(linewidth=0.7) + geom_point(size=1.6) +
  facet_wrap(~cohort, scales="free_x") +
  scale_color_manual(values=c("#0072B2","#0072B2","grey45","grey45"), guide="none") +
  scale_x_continuous(breaks=seq(2000,2020,4)) +
  labs(x="year", y="raw P(any fatal encounter)") +
  theme_minimal(base_size=12) + theme(panel.grid.minor=element_blank())
# label the two series per panel directly
lab <- do.call(rbind, lapply(unique(d$cohort), function(c){ s<-d[d$cohort==c & d$year==max(d$year),]; s }))
g <- g + ggrepel::geom_text_repel(data=lab, aes(label=grp), size=2.6, hjust=1, direction="y", show.legend=FALSE)
suppressMessages(library(ggrepel))
ggsave(file.path(FIG_DIR,"taskB5_raw_means.pdf"), g, width=10, height=4.5, device=cairo_pdf)
ggsave(file.path(FIG_DIR,"taskB5_raw_means.png"), g, width=10, height=4.5, dpi=200)

# quick numeric pre-trend read: correlation/divergence of treated vs control pre-period
for (g0 in c(2002,2009)) {
  s <- d[d$gg==g0 & d$year < g0,]
  tr <- s[grepl("treated \\(", s$grp),]; nt <- s[grepl("never", s$grp),]
  m <- merge(tr[,c("year","any.fatalities")], nt[,c("year","any.fatalities")], by="year")
  cat(sprintf("Cohort %d pre-period (%d yrs): mean treated %.3f vs control %.3f; mean gap %.3f, sd of gap %.3f\n",
      g0, nrow(m), mean(m$any.fatalities.x), mean(m$any.fatalities.y),
      mean(m$any.fatalities.x-m$any.fatalities.y), sd(m$any.fatalities.x-m$any.fatalities.y)))
}
cat("Wrote figures/taskB5_raw_means.{pdf,png} and output/taskB5_raw_means.csv\n")
