# =============================================================================
# d6_round5.R  --  D6 round-5 (final).
#
# NOTE: the Liu (2025) PDF (arXiv 2509.01829) was requested "in repo" but is not
# present in this checkout/remote, and arxiv.org is blocked by the egress proxy,
# so the paper could not be opened here. Item 1 uses the Lemma-1 simple-means
# formula EXACTLY AS STATED IN THE BRIEF:
#   Delta^Imp_{g,0} = (Ybar^T_{g-1} - Ybar^T_preavg) - (Ybar^C_{g-1} - Ybar^C_preavg)
#   C = initial control group = never-treated  union  {cohorts g' > g};
#   last pre-period = g-1 ; pre-avg = mean over observed pre-periods (2000..g-1);
#   unit-equal weighting; computed from simple means, independent of stage-1.
# This equals tau_imp(g,g) - tau_CS-NYT(g,g) when imputation reduces to the
# simple-means long-difference DiD (balanced panel); the FE-vs-simple-means
# discrepancy on the unbalanced sample is IMBALANCE, not inheritance.
#
# 1. gap_s1 vs Delta^Imp on (a) full droppers and (b) a balanced subpanel
#    (agencies with complete 2000-2020 outcomes). Relabel the gap figure:
#    at s=1 the remainder is imbalance.
# 2. Contribution table w_g*(tau_imp,g - tau_CS-NYT,g), cohort-average over post,
#    summing to the aggregate gap; reconcile with the 2013/2017/2018 localization.
# 3. testable_g = FALSE when n_treated_clusters == 1; recompute untestable weight
#    shares; redraw the main figure.
# 4. RI with placebo cohorts matched on pre-period outcome MEAN (alongside WCR).
# =============================================================================
suppressMessages({library(fixest); library(did); library(ggplot2); library(ggrepel); library(tidyr)})
DATA<-Sys.getenv("DIDREP_DATA",unset="data"); OUT<-Sys.getenv("DIDREP_OUT",unset="output"); FIG<-Sys.getenv("DIDREP_FIG",unset="figures")
FIRST_OUT<-2000; LAST<-2020
REV<-c("revere massachusetts","portsmouth new hampshire","memphis tennessee","fall river massachusetts","ramapo new york","springfield massachusetts")
code<-function(d){ d$agency.num<-as.numeric(factor(d$agency.id))
  d$yc<-ave(d$year.changed,d$agency.id,FUN=function(x){v<-x[!is.na(x)];if(length(v))min(v) else NA_real_})
  d$gname<-ifelse(is.na(d$yc)&d$no.req==0,0,ifelse(is.na(d$yc)&d$no.req==1,1987,d$yc)); d }
d<-code(read.csv(file.path(DATA,"dta.csv"))); d<-code(d[!(d$agency.id %in% REV),])
dw<-d[d$year>=FIRST_OUT & d$year<=LAST & !is.na(d$any.fatalities),]
cohorts<-sort(unique(dw$gname[dw$gname %in% 2001:2020]))

# ---- simple-means Delta^Imp (Lemma 1, brief statement), unit-equal ----------
simple_delta<-function(g, dsub){
  Tids<-unique(dsub$agency.id[dsub$gname==g]); Cids<-unique(dsub$agency.id[dsub$gname==0 | (dsub$gname>g & dsub$gname<3000)])
  pre<-FIRST_OUT:(g-1)
  arm<-function(ids){ s<-dsub[dsub$agency.id %in% ids & dsub$year %in% pre & !is.na(dsub$any.fatalities),]
    lp<-tapply(s$any.fatalities[s$year==g-1], s$agency.id[s$year==g-1], mean)
    pa<-tapply(s$any.fatalities, s$agency.id, mean); u<-intersect(names(lp),names(pa))
    if(!length(u)) return(NA_real_); mean(lp[u]-pa[u]) }
  arm(Tids)-arm(Cids) }

# ---- gap_s1 from the FE imputation vs CS-NYT, on a given subpanel -----------
gap_s1_fe<-function(dsub){
  tr<-dsub$gname %in% 2001:2020 & dsub$year>=dsub$gname
  un<-(dsub$gname==0)|(dsub$gname %in% 2001:2020 & dsub$year<dsub$gname)
  s1<-feols(any.fatalities~1|agency.num+year,dsub[un,],fixef.rm="none"); fe<-fixef(s1); A<-fe$agency.num; G<-fe$year
  y0<-function(ag,yr) unname(A[as.character(ag)]+G[as.character(yr)])
  set.seed(0); mp<-suppressWarnings(did::att_gt(yname="any.fatalities",tname="year",idname="agency.num",gname="gname",
      xformla=~1,control_group="notyettreated",base_period="varying",data=dsub))
  cs<-function(g){ i<-which(mp$group==g & mp$t==g); if(length(i)) mp$att[i] else NA_real_ }
  sapply(cohorts,function(g){ tg<-dsub[dsub$gname==g & dsub$year==g,]
    if(!nrow(tg)) return(NA_real_); mean(tg$any.fatalities-y0(tg$agency.num,tg$year))-cs(g) }) }

cat("=== ITEM 1: gap_s1 (FE) vs Delta^Imp (simple means) ===\n")
gap_full<-gap_s1_fe(dw); del_full<-sapply(cohorts,simple_delta,dsub=dw)
I1<-data.frame(cohort=cohorts, gap_s1=gap_full, delta_simple=del_full, imbalance=gap_full-del_full)
# (b) balanced subpanel: agencies observed every year 2000..2020
cnt<-tapply(dw$year,dw$agency.id,function(y) length(unique(y))); bal_ids<-names(cnt[cnt==(LAST-FIRST_OUT+1)])
db<-code(dw[dw$agency.id %in% bal_ids,])
cat(sprintf("balanced subpanel: %d of %d agencies have complete 2000-2020 outcomes\n", length(bal_ids), length(unique(dw$agency.id))))
bcoh<-sort(unique(db$gname[db$gname %in% 2001:2020]))
gv<-gap_s1_fe(db); names(gv)<-cohorts                          # compute once for balanced subpanel
I1b<-data.frame(cohort=bcoh, gap_s1_bal=gv[as.character(bcoh)],
  delta_bal=sapply(bcoh,simple_delta,dsub=db)); I1b$imbalance_bal<-I1b$gap_s1_bal-I1b$delta_bal
print(transform(I1, gap_s1=round(gap_s1,4), delta_simple=round(delta_simple,4), imbalance=round(imbalance,4)), row.names=FALSE)
cat(sprintf("\n(a) full droppers: max |imbalance remainder| = %.4f (mean |.| = %.4f)\n", max(abs(I1$imbalance),na.rm=TRUE), mean(abs(I1$imbalance),na.rm=TRUE)))
cat("(b) balanced subpanel:\n"); print(transform(I1b, gap_s1_bal=round(gap_s1_bal,5), delta_bal=round(delta_bal,5), imbalance_bal=signif(imbalance_bal,3)), row.names=FALSE)
cat(sprintf("    max |gap_s1_bal - delta_bal| = %.2e  (expect ~0 => simple means == FE on balanced)\n", max(abs(I1b$imbalance_bal),na.rm=TRUE)))
write.csv(I1, file.path(OUT,"d6_lemma1_full.csv"), row.names=FALSE); write.csv(I1b, file.path(OUT,"d6_lemma1_balanced.csv"), row.names=FALSE)

# ---- relabel gap figure: own = Delta^Imp (simple means), remainder = imbalance
long<-pivot_longer(data.frame(cohort=I1$cohort, block_bias=I1$delta_simple, imbalance=I1$imbalance),
  cols=c("block_bias","imbalance"), names_to="component", values_to="value")
long$component<-factor(long$component, levels=c("block_bias","imbalance"),
  labels=c("Delta^Imp block bias (Lemma 1, simple means)","imbalance (FE != simple means)"))
pg<-ggplot(long,aes(factor(cohort),value,fill=component))+geom_col()+
  geom_point(data=I1,aes(factor(cohort),gap_s1),inherit.aes=FALSE,shape=18,size=2.6)+
  geom_hline(yintercept=0,color="grey60")+
  scale_fill_manual(values=c("Delta^Imp block bias (Lemma 1, simple means)"="#0072B2","imbalance (FE != simple means)"="#E69F00"),name=NULL)+
  labs(x="cohort",y=expression(hat(Delta)^Imp[g*",0"]*"  (imputation - CS-NYT at event time 0)"),
    title="CS-imputation gap at s=1 = Lemma-1 block bias + imbalance",
    subtitle="diamond = total gap (FE); blue = simple-means Delta^Imp vs initial control group; orange = imbalance remainder (0 on a balanced panel)",
    caption="Lemma 1 (Liu 2025, arXiv:2509.01829): Delta^Imp_{g,0} = (Ybar^T_{g-1}-Ybar^T_preavg) - (Ybar^C_{g-1}-Ybar^C_preavg), C = never-treated union {g'>g}. On the balanced subpanel the orange term is ~0.")+
  theme_minimal(base_size=12)+theme(panel.grid.minor=element_blank(),legend.position="top",plot.caption=element_text(size=7,hjust=0,color="grey40"))
ggsave(file.path(FIG,"d6_gap_decomposition.pdf"),pg,width=9.5,height=5,device=cairo_pdf)
ggsave(file.path(FIG,"d6_gap_decomposition.png"),pg,width=9.5,height=5,dpi=200)
file.copy(file.path(FIG,"d6_gap_decomposition.png"),file.path(OUT,"d6_gap_decomposition.png"),overwrite=TRUE)

# ============================================================================
# ITEM 2: contribution table w_g * (tau_imp,g - tau_CS-NYT,g), cohort-avg
# ============================================================================
gid<-read.csv(file.path(OUT,"d6_gap_identity_s1.csv"))        # has gap_avg per cohort (imp - CS-NYT avg)
TAB<-read.csv(file.path(OUT,"weight_checkability_table.csv"))
wimp<-TAB[TAB$estimator=="imputation"&TAB$sample=="droppers",c("cohort","weight")]
C2<-merge(wimp, gid[,c("cohort","gap_avg")], by="cohort")
C2$contribution<-C2$weight*C2$gap_avg; C2<-C2[order(-abs(C2$contribution)),]
agg<-sum(C2$contribution)
write.csv(C2, file.path(OUT,"d6_gap_contributions.csv"), row.names=FALSE)
cat("\n=== ITEM 2: contribution w_g*(tau_imp - tau_CS-NYT) to the aggregate gap ===\n")
print(transform(C2, weight=round(weight,3), gap_avg=round(gap_avg,3), contribution=round(contribution,4)), row.names=FALSE)
cat(sprintf("aggregate gap Sum w_g*gap_avg_g = %.4f ; top-3 cohorts = %s contribute %.1f%%\n",
  agg, paste(C2$cohort[1:3],collapse="/"), 100*sum(abs(C2$contribution[1:3]))/sum(abs(C2$contribution))))

# ============================================================================
# ITEM 3: testable_g = FALSE when n_treated_clusters == 1
# ============================================================================
v2<-TAB[TAB$sample=="droppers"&TAB$estimator %in% c("CS","stacked","stacked_cengiz","imputation"),]
v2$testable_old<-v2$testable_g
v2$testable_g<-v2$testable_g & v2$n_g>1                         # <-- new rule
cat("\n=== ITEM 3: untestable weight share, old vs new (n_g==1 -> untestable) ===\n")
sh<-do.call(rbind,lapply(c("CS","stacked","stacked_cengiz","imputation"),function(e){ p<-v2[v2$estimator==e,]
  data.frame(estimator=e, untestable_old=sum(p$weight[!p$testable_old],na.rm=TRUE),
    untestable_new=sum(p$weight[!p$testable_g],na.rm=TRUE),
    n_single=sum(p$n_g==1 & p$testable_old)) }))
print(transform(sh, untestable_old=round(untestable_old,3), untestable_new=round(untestable_new,3)), row.names=FALSE)
write.csv(sh, file.path(OUT,"d6_untestable_ncluster.csv"), row.names=FALSE)

# redraw main figure (CR1 mde x-axis) with the new rule
mkfig<-function(sub,file){
  base<-c(CS="A: CS (simple)",stacked="B: stacked (R2 FE)",stacked_cengiz="B2: stacked Cengiz +/-4 (NT)",imputation="C: imputation (BJS)")
  ord<-c("CS","stacked","stacked_cengiz","imputation")
  sp<-function(w,pg) if(length(w)>=3) suppressWarnings(cor(w,pg,method="spearman")) else NA
  st<-sapply(ord,function(e){ p<-sub[sub$estimator==e,]
    sprintf("%s\nuntestable wt=%.2f  rho=%.2f",base[e],sum(p$weight[!p$testable_g],na.rm=TRUE),sp(p$weight,p$pre_g)) })
  sub$estimator<-factor(sub$estimator,levels=ord,labels=st)
  fin<-sub[sub$testable_g & is.finite(sub$mde_g),]; xmax<-max(fin$mde_g,na.rm=TRUE); x_un<-xmax*1.9
  sub$xpos<-ifelse(sub$testable_g & is.finite(sub$mde_g),sub$mde_g,x_un)
  g<-ggplot(sub,aes(xpos,weight))+
    annotate("rect",xmin=xmax*1.35,xmax=xmax*2.6,ymin=-Inf,ymax=Inf,fill="grey88",alpha=0.6)+
    geom_point(aes(size=n_g,fill=testable_g),shape=21,color="grey25",alpha=0.9)+
    ggrepel::geom_text_repel(aes(label=cohort),size=2.3,max.overlaps=30,seed=1)+
    annotate("text",x=sqrt(xmax*1.35*xmax*2.6),y=Inf,label="not testable",angle=90,vjust=1.1,hjust=1,size=2.6,color="grey35")+
    scale_fill_manual(values=c(`TRUE`="#0072B2",`FALSE`="#D55E00"),name="testable\n(n_g>1 & pre>=3)")+
    scale_size_area(max_size=8,name=expression(n[g]))+scale_x_log10()+facet_wrap(~estimator,scales="free_y",nrow=1)+
    labs(x="minimum detectable pre-trend slope  mde_g  (CR1; log; right strip = not testable)",y="cohort aggregation weight",
      title="Weight vs. checkability (droppers) -- testable requires >1 treated cluster",
      caption="testable_g=FALSE when n_treated_clusters=1 (no valid few-cluster inference) or pre_g<3 or no comparison group. Panel B2 (Cengiz): pre_g<4 excluded by window.")+
    theme_minimal(base_size=11)+theme(legend.position="right",panel.grid.minor=element_blank(),
      strip.text=element_text(size=8,lineheight=1.05),plot.caption=element_text(size=7,hjust=0,color="grey40"))
  ggsave(file.path(FIG,paste0(file,".pdf")),g,width=15,height=4.2,device=cairo_pdf)
  ggsave(file.path(FIG,paste0(file,".png")),g,width=15,height=4.2,dpi=200)
  file.copy(file.path(FIG,paste0(file,".png")),file.path(OUT,paste0(file,".png")),overwrite=TRUE) }
mkfig(v2,"fig_weight_checkability_d6_v2_ncluster")

# ============================================================================
# ITEM 4: RI matched on pre-period outcome MEAN
# ============================================================================
never_ids<-unique(dw$agency.id[dw$gname==0])
premean<-function(ids,g){ pre<-FIRST_OUT:(g-1); s<-dw[dw$year %in% pre & dw$agency.id %in% ids,]
  tapply(s$any.fatalities,s$agency.id,mean) }
ri_joint_meanmatch<-function(g,n_g,B=499,seed=4){
  pre<-FIRST_OUT:(g-1); base<-dw[dw$year %in% pre & dw$agency.id %in% never_ids,]
  tids<-unique(dw$agency.id[dw$gname==g]); tmean<-mean(premean(tids,g),na.rm=TRUE)
  nmean<-premean(never_ids,g); so<-rbind(base, dw[dw$year %in% pre & dw$agency.id %in% tids,]); so$D<-as.integer(so$agency.id %in% tids)
  Wobs<-tryCatch(unname(fixest::wald(feols(any.fatalities~i(year,D,ref=g-1)|agency.num+year,so,cluster=~agency.num),"year::",print=FALSE)$stat),error=function(e)NA)
  pool<-names(sort(abs(nmean-tmean))); set.seed(seed); Wb<-numeric(B)
  for(b in 1:B){ tr<-sample(head(pool,max(n_g*3,n_g+5)),n_g); dd<-base; dd$D<-as.integer(dd$agency.id %in% tr)
    mj<-tryCatch(feols(any.fatalities~i(year,D,ref=g-1)|agency.num+year,dd),error=function(e)NULL)
    Wb[b]<- if(is.null(mj)) NA else tryCatch(unname(fixest::wald(mj,"year::",print=FALSE)$stat),error=function(e)NA) }
  (1+sum(Wb>=Wobs,na.rm=TRUE))/(1+sum(!is.na(Wb))) }
rob<-read.csv(file.path(OUT,"d6_pretests_2009_2013_robust.csv")); szm<-read.csv(file.path(OUT,"d6_ri_sizematched.csv"))
cat("\n=== ITEM 4: joint RI matched on pre-period outcome mean (alongside WCR, size-match) ===\n")
I4<-do.call(rbind,lapply(c(2009,2013),function(g){ n_g<-length(unique(dw$agency.id[dw$gname==g]))
  data.frame(cohort=g, joint_p_WCR=round(mean(rob$joint_p_WCR[rob$cohort==g]),3),
    joint_RI_unmatched=szm$joint_p_RI_unmatched[szm$cohort==g], joint_RI_sizematched=szm$joint_p_RI_sizematched[szm$cohort==g],
    joint_RI_meanmatched=round(ri_joint_meanmatch(g,n_g),3)) }))
print(I4,row.names=FALSE); write.csv(I4,file.path(OUT,"d6_ri_meanmatched.csv"),row.names=FALSE)

cat("\nWrote d6_lemma1_{full,balanced}.csv, d6_gap_decomposition.*, d6_gap_contributions.csv,\n")
cat("  d6_untestable_ncluster.csv, fig_weight_checkability_d6_v2_ncluster.*, d6_ri_meanmatched.csv\n")
