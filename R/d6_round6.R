# =============================================================================
# d6_round6.R  --  checkability threshold pre_g>=3 -> pre_g>=2, and the
# unchecked-weight diagnostic U^k(deltabar) (paper section 3.2).
# No new estimators. Reuses committed weights/estimates + round-4 WCR machinery.
#
# MDE^k_g = Zc * SE_WCR(slope_g)  if pre_g>=2 & comparison exists & n_clusters>=2 ; else Inf.
# U^k(d) = sum_g W_g 1{MDE_g > d}.  Break-even deltabar_g = theta*/a_g,
#   a_g = mean(post event times) - mean(baseline event times) (bias of a unit slope).
# =============================================================================
suppressMessages({library(fixest); library(ggplot2); library(ggrepel)})
DATA<-Sys.getenv("DIDREP_DATA",unset="data"); OUT<-Sys.getenv("DIDREP_OUT",unset="output"); FIG<-Sys.getenv("DIDREP_FIG",unset="figures")
FIRST_OUT<-2000; LAST<-2020; WINDOW<-4; Zc<-qnorm(0.975)+qnorm(0.80)
WEBB<-c(-sqrt(1.5),-1,-sqrt(.5),sqrt(.5),1,sqrt(1.5))
REV<-c("revere massachusetts","portsmouth new hampshire","memphis tennessee","fall river massachusetts","ramapo new york","springfield massachusetts")
code<-function(d){ d$agency.num<-as.numeric(factor(d$agency.id))
  d$yc<-ave(d$year.changed,d$agency.id,FUN=function(x){v<-x[!is.na(x)];if(length(v))min(v) else NA_real_})
  d$gname<-ifelse(is.na(d$yc)&d$no.req==0,0,ifelse(is.na(d$yc)&d$no.req==1,1987,d$yc)); d }
D_DECK<-code(read.csv(file.path(DATA,"dta.csv"))); D_DROP<-code(D_DECK[!(D_DECK$agency.id %in% REV),])
dwf<-function(d) d[d$year>=FIRST_OUT & d$year<=LAST & !is.na(d$any.fatalities),]
SAMP<-list(deck=dwf(D_DECK), droppers=dwf(D_DROP))
never_ids<-function(d) unique(d$agency.id[d$gname==0])
nochange_ids<-function(d) unique(d$agency.id[d$change.type=="No Change"])
notyet_ids<-function(d,g) unique(d$agency.id[d$gname>g & d$gname<3000])

# ---- round-4 WCR slope SE (verbatim) + CR1 SE -------------------------------
cr1_vcov<-function(X,resid,cl,XtXinv){ Gs<-unique(cl);p<-ncol(X);meat<-matrix(0,p,p)
  for(g in Gs){ix<-cl==g;Xg<-X[ix,,drop=FALSE];sc<-crossprod(Xg,resid[ix]);meat<-meat+tcrossprod(sc)}
  nG<-length(Gs);n<-nrow(X);adj<-(nG/(nG-1))*((n-1)/(n-p));XtXinv%*%(adj*meat)%*%XtXinv }
wcr_slope_se<-function(s,B=999,seed=1){ X<-model.matrix(~year*D,s);y<-s$any.fatalities;cl<-s$agency.num;k<-which(colnames(X)=="year:D")
  if(length(unique(cl))<3||length(unique(s$D))<2) return(NA_real_)
  XtXinv<-tryCatch(solve(crossprod(X)),error=function(e)NULL); if(is.null(XtXinv)) return(NA_real_)
  Xr<-X[,-k,drop=FALSE];br<-solve(crossprod(Xr))%*%crossprod(Xr,y);yhat_r<-Xr%*%br;ur<-y-yhat_r
  Gs<-unique(cl);set.seed(seed);bs<-numeric(B)
  for(b in 1:B){ wg<-sample(WEBB,length(Gs),replace=TRUE);names(wg)<-Gs;ystar<-yhat_r+wg[as.character(cl)]*ur;bs[b]<-(XtXinv%*%crossprod(X,ystar))[k] }
  sd(bs) }
cr1_slope_se<-function(s){ m<-tryCatch(feols(any.fatalities~year*D,s,cluster=~agency.num),error=function(e)NULL)
  if(is.null(m)) return(NA_real_); ct<-coeftable(m); if(!("year:D"%in%rownames(ct))) return(NA_real_); ct["year:D","Std. Error"] }
pre_sub<-function(d,g,tids,cids){ pre<-FIRST_OUT:(g-1); s<-d[d$year%in%pre & d$agency.id%in%c(tids,cids),]; s$D<-as.integer(s$agency.id%in%tids); s }

# comparison group + has-control per estimator
comp_of<-function(est,d,g){ switch(est,
  `CS-NT`=never_ids(d), `CS-NYT`=union(never_ids(d),notyet_ids(d,g)),
  imputation=union(never_ids(d),notyet_ids(d,g)), stacked=nochange_ids(d), Cengiz=never_ids(d)) }
has_ctrl<-function(est,g){ switch(est,
  `CS-NT`=TRUE,`CS-NYT`=TRUE,imputation=TRUE,
  stacked=(g-WINDOW>=FIRST_OUT), Cengiz=(g-WINDOW>=FIRST_OUT & g+WINDOW<=LAST)) }

# a_g : mean(post event times) - mean(baseline event times)
a_generic<-function(est,g){ post<-0:(LAST-g)
  if(est=="Cengiz"){ post<-0:WINDOW; base<- -(WINDOW:1) }
  else if(est%in%c("imputation","stacked")) base<- -((g-FIRST_OUT):1)      # pre-average baseline
  else base<- -1                                                           # g-1 baseline (CS-NT/NYT)
  mean(post)-mean(base) }
a_closed<-function(est,g){ pre_g<-g-FIRST_OUT; post_g<-LAST-g+1
  if(est=="Cengiz") WINDOW+0.5 else if(est%in%c("imputation","stacked")) (pre_g+post_g)/2 else (post_g+1)/2 }

# ---- per-cohort weights/pre_g/post_g/n_g per estimator x sample -------------
TAB<-read.csv(file.path(OUT,"weight_checkability_table.csv"))
pull<-function(tabest,samp){ t<-TAB[TAB$estimator==tabest & TAB$sample==samp,c("cohort","n_g","pre_g","post_g","weight","mde_g")]; t[order(t$cohort),] }
PANELS<-list(
  `CS-NT|deck`=pull("CS","deck"), `CS-NT|droppers`=pull("CS","droppers"),
  `CS-NYT|deck`=pull("CS","deck"), `CS-NYT|droppers`=pull("CS","droppers"),   # CS-NYT weights = CS simple (control-independent; verified in d6_verify)
  `imputation|droppers`=pull("imputation","droppers"),
  `stacked|deck`=pull("stacked","deck"), `stacked|droppers`=pull("stacked","droppers"),
  `Cengiz|deck`=pull("stacked_cengiz","deck"), `Cengiz|droppers`=pull("stacked_cengiz","droppers"))

# ---- compute MDE (WCR) per cell ---------------------------------------------
rows<-list()
for(nm in names(PANELS)){ parts<-strsplit(nm,"\\|")[[1]]; est<-parts[1]; samp<-parts[2]; d<-SAMP[[samp]]; P<-PANELS[[nm]]
  for(i in seq_len(nrow(P))){ g<-P$cohort[i]; n_g<-P$n_g[i]; pre_g<-P$pre_g[i]; W<-P$weight[i]
    hc<-has_ctrl(est,g)
    mde_raw<-NA_real_
    if(pre_g>=2 && hc && n_g>=2){ tids<-unique(d$agency.id[d$gname==g]); s<-pre_sub(d,g,tids,comp_of(est,d,g))
      mde_raw<-Zc*wcr_slope_se(s) }
    # testability is RULE-based (pre_g, comparison group, n_clusters>=2) -- matches the
    # canonical untestable_share. MDE = Zc*SE_WCR where rule-eligible (NA if WCR fails);
    # Inf where rule-ineligible. A rule-eligible cohort with NA WCR is TESTABLE but its MDE
    # is uncomputable, so it is never "checked" in U(deltabar) (flagged below).
    rule_old<-(pre_g>=3 && hc && n_g>=2); rule_new<-(pre_g>=2 && hc && n_g>=2)
    MDE<-if(rule_new) mde_raw else Inf          # NA if rule-eligible but WCR failed; Inf otherwise
    ag<-a_generic(est,g); ac<-a_closed(est,g)
    rows[[length(rows)+1L]]<-data.frame(estimator=est,sample=samp,cohort=g,weight=W,pre_g=pre_g,post_g=P$post_g[i],
      n_g=n_g,has_ctrl=hc,mde_wcr_raw=mde_raw,MDE=MDE,a_g=ag,a_closed=ac,
      wcr_na=(rule_new && !is.finite(mde_raw)),
      testable_old=rule_old, testable_new=rule_new, row.names=NULL) } }
R6<-do.call(rbind,rows)
R6$a_flag<-ifelse(abs(R6$a_g-R6$a_closed)>1e-9,"MISMATCH","ok")

# =============================================================================
# VALIDATION 1: reproduce d6_mde_wcr.csv for already-testable (n_g>=2,pre_g>=3)
# =============================================================================
mw<-read.csv(file.path(OUT,"d6_mde_wcr.csv")); mw$est<-c(CS="CS-NT",stacked="stacked",stacked_cengiz="Cengiz",imputation="imputation")[mw$estimator]
chk<-merge(R6[R6$sample=="droppers" & R6$testable_old,c("estimator","cohort","mde_wcr_raw")],
           mw[,c("est","cohort","mde_wcr")], by.x=c("estimator","cohort"), by.y=c("est","cohort"))
cat(sprintf("[V1] reproduce d6_mde_wcr (testable droppers, n=%d): max|diff| = %.2e\n", nrow(chk), max(abs(chk$mde_wcr_raw-chk$mde_wcr),na.rm=TRUE)))
wna<-R6[R6$wcr_na,]
cat(sprintf("[WCR-NA] rule-testable cohorts whose WCR MDE is uncomputable (-> never checked in U): %s\n",
  if(nrow(wna)) paste(sprintf("%s|%s|%d(n_g=%d,w=%.3f)",wna$estimator,wna$sample,wna$cohort,wna$n_g,wna$weight),collapse=", ") else "none"))

# =============================================================================
# ITEM 1: new MDEs for pre_g==2 cohorts (should be 2002 only in droppers)
# =============================================================================
new2<-R6[R6$pre_g==2 & R6$n_g>=2,]
cat("\n[ITEM 1] pre_g==2 cohorts with >=2 treated clusters:\n")
print(transform(new2[,c("estimator","sample","cohort","n_g","has_ctrl","mde_wcr_raw","MDE")],
  mde_wcr_raw=round(mde_wcr_raw,4), MDE=round(MDE,4)), row.names=FALSE)
ntc_2002<-4; cat(sprintf("\n2002 treated clusters (droppers) = %d. Distinct Webb draws: total 6^(G) with G=#clusters (effectively unlimited);\n  but only %d TREATED clusters -> 6^%d = %d distinct treated-side weight patterns. FLAG: WCR UNRELIABLE with <~10 treated clusters (severe below 5).\n",
  ntc_2002, ntc_2002, ntc_2002, 6^ntc_2002))
cat(sprintf("stacked 2002: has_ctrl=%s -> MDE=%s (confirm Inf)\n",
  R6$has_ctrl[R6$estimator=="stacked"&R6$sample=="droppers"&R6$cohort==2002],
  R6$MDE[R6$estimator=="stacked"&R6$sample=="droppers"&R6$cohort==2002]))

# =============================================================================
# ITEM 2: U(inf) old (pre_g>=3) vs new (pre_g>=2)
# =============================================================================
Uinf<-do.call(rbind,lapply(names(PANELS),function(nm){ parts<-strsplit(nm,"\\|")[[1]]; p<-R6[R6$estimator==parts[1]&R6$sample==parts[2],]
  data.frame(estimator=parts[1],sample=parts[2],
    U_inf_old=sum(p$weight[!p$testable_old]), U_inf_new=sum(p$weight[!p$testable_new]),
    diff=sum(p$weight[!p$testable_old])-sum(p$weight[!p$testable_new]),
    newly_testable=paste(p$cohort[p$testable_new & !p$testable_old],collapse=";"), row.names=NULL) }))
cat("\n[ITEM 2] U(inf): old (pre_g>=3) vs new (pre_g>=2)\n")
print(transform(Uinf,U_inf_old=round(U_inf_old,3),U_inf_new=round(U_inf_new,3),diff=round(diff,4)),row.names=FALSE)

# =============================================================================
# ITEM 3: U at break-even deltabar_g = theta*/a_g, theta* in {0.05,0.10}
# =============================================================================
U_break<-list(); PCT<-list()
for(th in c(0.05,0.10)){ for(nm in names(PANELS)){ parts<-strsplit(nm,"\\|")[[1]]; p<-R6[R6$estimator==parts[1]&R6$sample==parts[2],]
  p$deltabar<-th/p$a_g; p$checked<-is.finite(p$MDE) & p$MDE<=p$deltabar
  U_break[[length(U_break)+1L]]<-data.frame(estimator=parts[1],sample=parts[2],theta_star=th,
    U_breakeven=sum(p$weight[!p$checked]), n_checked=sum(p$checked), row.names=NULL)
  PCT[[length(PCT)+1L]]<-transform(p[,c("estimator","sample","cohort","weight","pre_g","post_g","a_g","MDE","deltabar","checked")], theta_star=th) } }
UB<-do.call(rbind,U_break); PCTtab<-do.call(rbind,PCT)
cat("\n[ITEM 3] U at break-even tolerance:\n")
print(transform(UB,U_breakeven=round(U_breakeven,3)),row.names=FALSE)

# =============================================================================
# ITEM 4: U curve over deltabar grid [0,0.06] step 0.001
# =============================================================================
grid<-seq(0,0.06,by=0.001)
CURVE<-do.call(rbind,lapply(names(PANELS),function(nm){ parts<-strsplit(nm,"\\|")[[1]]; p<-R6[R6$estimator==parts[1]&R6$sample==parts[2],]
  data.frame(estimator=parts[1],sample=parts[2],deltabar=grid,
    U=sapply(grid,function(d) sum(p$weight[!(is.finite(p$MDE)&p$MDE<=d)])), row.names=NULL) }))
# validation 2: non-increasing and >= U(inf)
v2ok<-all(tapply(CURVE$U,paste(CURVE$estimator,CURVE$sample),function(u) all(diff(u)<=1e-12)))
v2ge<-all(mapply(function(e,s){ u<-CURVE$U[CURVE$estimator==e&CURVE$sample==s]; uinf<-Uinf$U_inf_new[Uinf$estimator==e&Uinf$sample==s]; all(u>=uinf-1e-12) },
  Uinf$estimator,Uinf$sample))
cat(sprintf("\n[V2] U curve non-increasing: %s ; U(d) >= U(inf_new) everywhere: %s\n",v2ok,v2ge))

write.csv(R6, file.path(OUT,"d6_round6_mde.csv"), row.names=FALSE)
write.csv(rbind(cbind(kind="U_inf",estimator=Uinf$estimator,sample=Uinf$sample,theta_star=NA,value=Uinf$U_inf_new,note=paste0("old=",round(Uinf$U_inf_old,4))),
                cbind(kind="U_breakeven",estimator=UB$estimator,sample=UB$sample,theta_star=UB$theta_star,value=round(UB$U_breakeven,4),note="")),
          file.path(OUT,"d6_unchecked_weight.csv"), row.names=FALSE)
write.csv(PCTtab, file.path(OUT,"d6_unchecked_percohort.csv"), row.names=FALSE)
write.csv(CURVE, file.path(OUT,"d6_unchecked_curve.csv"), row.names=FALSE)

# =============================================================================
# VALIDATION 3 & 4
# =============================================================================
cat(sprintf("\n[V3] new U(inf) <= old for all: %s ; diff == weight on newly-testable: %s\n",
  all(Uinf$U_inf_new<=Uinf$U_inf_old+1e-12),
  all(abs(Uinf$diff - mapply(function(e,s){ p<-R6[R6$estimator==e&R6$sample==s,]; sum(p$weight[p$testable_new&!p$testable_old]) },Uinf$estimator,Uinf$sample))<1e-9)))
sw<-tapply(R6$weight,paste(R6$estimator,R6$sample),sum)
cat(sprintf("[V4] weights sum to 1 (max|.-1|)=%.2e (CS-NYT reuses CS weights; stacked includes w_s=0 cohorts)\n", max(abs(sw-1))))

# =============================================================================
# ITEM 5: refresh fig_weight_checkability_d6_v2_ncluster with pre_g>=2
# (same CR1 mde_g axis; 2002 now plotted where testable_new)
# =============================================================================
fig<-TAB[TAB$sample=="droppers" & TAB$estimator%in%c("CS","stacked","stacked_cengiz","imputation"),]
fig$testable_g <- fig$n_g>1 & fig$pre_g>=2 &
  ifelse(fig$estimator=="stacked", fig$cohort>=2004, TRUE) &
  ifelse(fig$estimator=="stacked_cengiz", fig$cohort>=2004 & fig$cohort<=2016, TRUE)
# CR1 mde for 2002 (CS & imputation; never / never+notyet comparison) on droppers
d<-SAMP[["droppers"]]; t2<-unique(d$agency.id[d$gname==2002])
mde2_cs<-Zc*cr1_slope_se(pre_sub(d,2002,t2,never_ids(d)))
mde2_imp<-Zc*cr1_slope_se(pre_sub(d,2002,t2,union(never_ids(d),notyet_ids(d,2002))))
fig$mde_g[fig$estimator=="CS"&fig$cohort==2002]<-mde2_cs
fig$mde_g[fig$estimator=="imputation"&fig$cohort==2002]<-mde2_imp
cat(sprintf("\n[ITEM 5] 2002 CR1 mde: CS=%.4f, imputation=%.4f (now on-axis); stacked/Cengiz 2002 stay untestable\n",mde2_cs,mde2_imp))
base<-c(CS="A: CS (simple)",stacked="B: stacked (R2 FE)",stacked_cengiz="B2: stacked Cengiz +/-4 (NT)",imputation="C: imputation (BJS)")
ord<-c("CS","stacked","stacked_cengiz","imputation"); sp<-function(w,pg) if(length(w)>=3) suppressWarnings(cor(w,pg,method="spearman")) else NA
st<-sapply(ord,function(e){ p<-fig[fig$estimator==e,]; sprintf("%s\nuntestable wt=%.2f  rho=%.2f",base[e],sum(p$weight[!p$testable_g],na.rm=TRUE),sp(p$weight,p$pre_g)) })
fig$estimator<-factor(fig$estimator,levels=ord,labels=st)
fin<-fig[fig$testable_g & is.finite(fig$mde_g),]; xmax<-max(fin$mde_g,na.rm=TRUE); x_un<-xmax*1.9
fig$xpos<-ifelse(fig$testable_g & is.finite(fig$mde_g),fig$mde_g,x_un)
g<-ggplot(fig,aes(xpos,weight))+
  annotate("rect",xmin=xmax*1.35,xmax=xmax*2.6,ymin=-Inf,ymax=Inf,fill="grey88",alpha=0.6)+
  geom_point(aes(size=n_g,fill=testable_g),shape=21,color="grey25",alpha=0.9)+
  ggrepel::geom_text_repel(aes(label=cohort),size=2.3,max.overlaps=30,seed=1)+
  annotate("text",x=sqrt(xmax*1.35*xmax*2.6),y=Inf,label="not testable",angle=90,vjust=1.1,hjust=1,size=2.6,color="grey35")+
  scale_fill_manual(values=c(`TRUE`="#0072B2",`FALSE`="#D55E00"),name="testable\n(n_g>1 & pre_g>=2)")+
  scale_size_area(max_size=8,name=expression(n[g]))+scale_x_log10()+facet_wrap(~estimator,scales="free_y",nrow=1)+
  labs(x="minimum detectable pre-trend slope  mde_g  (CR1; log; right strip = not testable)",y="cohort aggregation weight",
    title="Weight vs. checkability (droppers) -- checkability rule pre_g>=2",
    caption="pre_g>=2 (was >=3): 2002 now testable for CS and imputation (4 treated clusters -- WCR MDE unreliable, see RESULTS). testable_g also needs n_g>1 and a comparison group.")+
  theme_minimal(base_size=11)+theme(legend.position="right",panel.grid.minor=element_blank(),strip.text=element_text(size=8,lineheight=1.05),plot.caption=element_text(size=7,hjust=0,color="grey40"))
ggsave(file.path(FIG,"fig_weight_checkability_d6_v2_ncluster.pdf"),g,width=15,height=4.2,device=cairo_pdf)
ggsave(file.path(FIG,"fig_weight_checkability_d6_v2_ncluster.png"),g,width=15,height=4.2,dpi=200)
file.copy(file.path(FIG,"fig_weight_checkability_d6_v2_ncluster.png"),file.path(OUT,"fig_weight_checkability_d6_v2_ncluster.png"),overwrite=TRUE)

# =============================================================================
# ITEM 6: canonical update + change log
# =============================================================================
CN<-read.csv(file.path(OUT,"d6_canonical_numbers.csv"),stringsAsFactors=FALSE)
changelog<-list(); logchg<-function(key,old,new,cause) changelog[[length(changelog)+1L]]<<-data.frame(key=key,old=old,new=new,cause=cause,row.names=NULL)
# update untestable_share (now U_inf_new), n_testable, testable_list
for(i in seq_len(nrow(Uinf))){ r<-Uinf[i,]; key<-sprintf("%s|%s",r$estimator,r$sample)
  # only CS-NT/CS-NYT/imputation/stacked/Cengiz are in canonical with these exact names
  cn_key<-key
  oldrow<-CN[CN$section=="untestable_share" & CN$key==cn_key,]
  if(nrow(oldrow)){ if(round(as.numeric(oldrow$value[1]),3)!=round(r$U_inf_new,3)) logchg(sprintf("untestable_share|%s",key),oldrow$value[1],round(r$U_inf_new,3),"rule pre_g>=3 -> pre_g>=2 (2002 now testable)")
    CN$value[CN$section=="untestable_share"&CN$key==cn_key]<-round(r$U_inf_new,3) }
  p<-R6[R6$estimator==r$estimator&R6$sample==r$sample,]; tl<-paste(sort(p$cohort[p$testable_new]),collapse=";"); nt<-sum(p$testable_new)
  ntrow<-CN[CN$section=="n_testable" & CN$key==cn_key,]
  if(nrow(ntrow)){ if(as.numeric(ntrow$value[1])!=nt) logchg(sprintf("n_testable|%s",key),ntrow$value[1],nt,"rule pre_g>=2: +2002 where testable")
    CN$value[CN$section=="n_testable"&CN$key==cn_key]<-nt; CN$note[CN$section=="n_testable"&CN$key==cn_key]<-tl } }
# append new sections
addrows<-list(); ADD<-function(section,key,value,note="") addrows[[length(addrows)+1L]]<<-data.frame(section=section,key=key,value=value,note=note,row.names=NULL)
for(i in seq_len(nrow(Uinf))){ r<-Uinf[i,]; ADD("U_inf_old",sprintf("%s|%s",r$estimator,r$sample),round(r$U_inf_old,4)); ADD("U_inf_new",sprintf("%s|%s",r$estimator,r$sample),round(r$U_inf_new,4),sprintf("newly=%s",r$newly_testable)) }
for(i in seq_len(nrow(UB))){ r<-UB[i,]; ADD("U_breakeven",sprintf("%s|%s|theta%.2f",r$estimator,r$sample,r$theta_star),round(r$U_breakeven,4)) }
for(i in seq_len(nrow(new2))){ r<-new2[i,]; ADD("MDE_2002",sprintf("%s|%s",r$estimator,r$sample),round(r$MDE,4),sprintf("n_g=%d; WCR unreliable (4 treated clusters)",r$n_g)) }
# a_g (report per estimator, generic==closed)
for(est in c("CS-NT","CS-NYT","imputation","stacked","Cengiz")){ p<-R6[R6$estimator==est & R6$sample=="droppers",]
  ADD("a_g",sprintf("%s|formula",est), paste0(ifelse(est=="Cengiz","k+1/2",ifelse(est%in%c("imputation","stacked"),"(pre+post)/2","(post+1)/2"))),
      sprintf("flag=%s; 2009=%.1f 2013=%.1f 2017=%.1f", unique(p$a_flag)[1], p$a_g[p$cohort==2009][1], p$a_g[p$cohort==2013][1], p$a_g[p$cohort==2017][1])) }
CN2<-rbind(CN,do.call(rbind,addrows)); write.csv(CN2,file.path(OUT,"d6_canonical_numbers.csv"),row.names=FALSE)
CL<-do.call(rbind,changelog); write.csv(CL,file.path(OUT,"d6_round6_changelog.csv"),row.names=FALSE)
cat("\n[ITEM 6] canonical change log:\n"); print(CL,row.names=FALSE)
cat(sprintf("\na_g flag summary: %s\n", paste(unique(R6$a_flag),collapse=",")))
cat("\nWrote d6_round6_mde.csv, d6_unchecked_weight.csv, d6_unchecked_percohort.csv, d6_unchecked_curve.csv,\n")
cat("  d6_round6_changelog.csv, updated d6_canonical_numbers.csv, refreshed fig_weight_checkability_d6_v2_ncluster.*\n")
