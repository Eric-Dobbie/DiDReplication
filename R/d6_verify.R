# =============================================================================
# d6_verify.R  --  D6 verification pass: canonical numbers of record.
# No new estimators/specs. Reuses committed per-cohort estimates; recomputes the
# pre-test battery for 2009/2013/2017 (same functions/seeds as rounds 3/5 so
# 2009/2013 reproduce and are cross-checked), ATTs+SE, CS-NYT, untestable shares,
# rho, and the s=1 gap identity. Writes d6_canonical_numbers.csv + RESULTS_d6_canonical.md.
# =============================================================================
suppressMessages({library(did); library(fixest); library(lfe)})
DATA<-Sys.getenv("DIDREP_DATA",unset="data"); OUT<-Sys.getenv("DIDREP_OUT",unset="output")
FIRST_OUT<-2000; LAST<-2020; WINDOW<-4; Zc<-qnorm(0.975)+qnorm(0.80)
WEBB<-c(-sqrt(1.5),-1,-sqrt(.5),sqrt(.5),1,sqrt(1.5))
REV<-c("revere massachusetts","portsmouth new hampshire","memphis tennessee","fall river massachusetts","ramapo new york","springfield massachusetts")
code<-function(d){ d$agency.num<-as.numeric(factor(d$agency.id))
  d$yc<-ave(d$year.changed,d$agency.id,FUN=function(x){v<-x[!is.na(x)];if(length(v))min(v) else NA_real_})
  d$gname<-ifelse(is.na(d$yc)&d$no.req==0,0,ifelse(is.na(d$yc)&d$no.req==1,1987,d$yc)); d }
dta<-code(read.csv(file.path(DATA,"dta.csv"))); dta_drop<-code(dta[!(dta$agency.id %in% REV),])
dw<-dta_drop[dta_drop$year>=FIRST_OUT & dta_drop$year<=LAST & !is.na(dta_drop$any.fatalities),]
never_ids<-function(d) unique(d$agency.id[d$gname==0])
nochange_ids<-function(d) unique(d$agency.id[d$change.type=="No Change"])
notyet_ids<-function(d,g) unique(d$agency.id[d$gname>g & d$gname<3000])

# ---- CS (gradient simple weights + per-cohort tau + overall ATT/SE) ---------
cs_panel<-function(d, cg){ set.seed(0)
  mp<-suppressWarnings(did::att_gt(yname="any.fatalities",tname="year",idname="agency.num",gname="gname",
      xformla=~1,control_group=cg,clustervars="agency.num",base_period="universal",data=d))
  ag<-aggte(mp,type="simple"); base<-ag$overall.att; n<-length(mp$att); w<-numeric(n); dd<-1e-4
  for(k in seq_len(n)){ m2<-mp; m2$att[k]<-m2$att[k]+dd; w[k]<-(aggte(m2,type="simple")$overall.att-base)/dd }
  gg<-mp$group; wg<-tapply(w,gg,sum); tau<-tapply(w*mp$att,gg,sum)/wg; coh<-as.integer(names(wg))
  ng<-sapply(coh,function(g) length(unique(d$agency.id[d$gname==g & !is.na(d$any.fatalities)])))
  list(tab=data.frame(cohort=coh,weight=as.numeric(wg),tau=as.numeric(tau),n_g=ng,pre_g=coh-FIRST_OUT,row.names=NULL),
       att=base, se=ag$overall.se) }

# ---- stacked (R2) + Cengiz +/-4 builders, pooled ATT/SE/weights -------------
build_stack<-function(d){ coh<-sort(unique(d$year.changed[!is.na(d$year.changed)&d$year.changed>=FIRST_OUT]))
  nc<-d[d$change.type=="No Change",]; parts<-list()
  for(g in coh){ tr<-d[!is.na(d$year.changed)&d$year.changed==g,]; tr$cohort<-g; tr$treat<-1L; parts[[length(parts)+1L]]<-tr
    if(g-WINDOW>=FIRST_OUT){ ct<-nc; ct$cohort<-g; ct$treat<-0L; parts[[length(parts)+1L]]<-ct } }
  r<-do.call(rbind,parts); r$year.cohort<-as.numeric(paste0(r$year,r$cohort)); r$agency_stack<-paste(r$agency.id,r$cohort,sep="__"); r }
build_cengiz<-function(d){ K<-WINDOW; nt<-d[d$agency.id %in% never_ids(d),]
  coh<-sort(unique(d$year.changed[!is.na(d$year.changed)&d$year.changed>=FIRST_OUT])); coh<-coh[coh-K>=FIRST_OUT & coh+K<=LAST]
  parts<-list(); for(g in coh){ tr<-d[!is.na(d$year.changed)&d$year.changed==g & d$year>=g-K & d$year<=g+K,]
    if(!nrow(tr)) next; tr$cohort<-g; tr$treat<-1L; ct<-nt[nt$year>=g-K & nt$year<=g+K,]; ct$cohort<-g; ct$treat<-0L; parts[[length(parts)+1L]]<-rbind(tr,ct) }
  r<-do.call(rbind,parts); r$year.cohort<-as.numeric(paste0(r$year,r$cohort)); r$agency_stack<-paste(r$agency.id,r$cohort,sep="__"); r }
stack_panel<-function(rec){ est<-rec[!is.na(rec$any.fatalities),]
  m<-felm(any.fatalities~no.req|agency_stack+year.cohort|0|agency.id,data=est)
  att<-unname(coef(m)["no.req"]); se<-unname(sqrt(diag(vcov(m)))["no.req"])
  r1<-felm(no.req~1|agency_stack+year.cohort,data=est); est$vhat<-as.numeric(r1$residuals[,1])
  stacks<-sort(unique(est$cohort)); Vs<-sapply(stacks,function(s) sum(est$vhat[est$cohort==s]^2)); ws<-Vs/sum(Vs)
  betas<-sapply(stacks,function(s){ ds<-est[est$cohort==s,]; mm<-tryCatch(felm(any.fatalities~no.req|agency.id+year,ds),error=function(e)NULL)
    if(is.null(mm)) NA_real_ else unname(coef(mm)["no.req"]) })
  ng<-sapply(stacks,function(s) length(unique(est$agency.id[est$cohort==s & est$treat==1])))
  ncg<-sapply(stacks,function(s) length(unique(est$agency.id[est$cohort==s & est$treat==0])))
  list(tab=data.frame(cohort=stacks,weight=ws,tau=betas,n_g=ng,n_c=ncg,pre_g=stacks-FIRST_OUT,row.names=NULL),att=att,se=se) }

# ---- imputation (BJS) ATT + per-cohort + cluster-bootstrap SE ---------------
imp_panel<-function(d,B=499,seed=7){ d<-d[d$year>=FIRST_OUT&d$year<=LAST&!is.na(d$any.fatalities),]
  fit<-function(dd){ tr<-dd$gname%in%2001:2020 & dd$year>=dd$gname; un<-(dd$gname==0)|(dd$gname%in%2001:2020 & dd$year<dd$gname)
    s1<-feols(any.fatalities~1|agency.num+year,dd[un,],fixef.rm="none"); t<-dd[tr,]; t$yhat0<-predict(s1,newdata=t); t<-t[!is.na(t$yhat0),]
    t$tau<-t$any.fatalities-t$yhat0; t }
  t<-fit(d); att<-mean(t$tau); N1<-nrow(t)
  coh<-sort(unique(t$gname)); wtab<-data.frame(cohort=coh,
    weight=sapply(coh,function(g) sum(t$gname==g)/N1), tau=sapply(coh,function(g) mean(t$tau[t$gname==g])),
    n_g=sapply(coh,function(g) length(unique(t$agency.id[t$gname==g]))), pre_g=coh-FIRST_OUT, row.names=NULL)
  # cluster bootstrap SE (resample agencies with replacement)
  ids<-unique(d$agency.id); set.seed(seed); bb<-numeric(B)
  for(b in 1:B){ samp<-sample(ids,length(ids),replace=TRUE)
    db<-do.call(rbind,lapply(seq_along(samp),function(i){ x<-d[d$agency.id==samp[i],]; x$agency.id<-paste0(x$agency.id,"_",i); x }))
    db<-code(db); tb<-tryCatch(fit(db),error=function(e)NULL); bb[b]<-if(is.null(tb)||!nrow(tb)) NA else mean(tb$tau) }
  list(tab=wtab, att=att, se=sd(bb,na.rm=TRUE)) }

# ---- pre-test battery (CR1 / WCR / RI unmatched / RI mean-matched) ----------
cr1_vcov<-function(X,resid,cl,XtXinv){ Gs<-unique(cl);p<-ncol(X);meat<-matrix(0,p,p)
  for(g in Gs){ix<-cl==g;Xg<-X[ix,,drop=FALSE];sc<-crossprod(Xg,resid[ix]);meat<-meat+tcrossprod(sc)}
  nG<-length(Gs);n<-nrow(X);adj<-(nG/(nG-1))*((n-1)/(n-p));XtXinv%*%(adj*meat)%*%XtXinv }
pre_sub<-function(g,tids,cids){ pre<-FIRST_OUT:(g-1); s<-dw[dw$year%in%pre & dw$agency.id%in%c(tids,cids),]; s$D<-as.integer(s$agency.id%in%tids); s }
wcr_slope<-function(s,B=999,seed=1,want="p"){ X<-model.matrix(~year*D,s);y<-s$any.fatalities;cl<-s$agency.num;k<-which(colnames(X)=="year:D")
  XtXinv<-solve(crossprod(X));bhat<-XtXinv%*%crossprod(X,y);V<-cr1_vcov(X,y-X%*%bhat,cl,XtXinv);t_obs<-bhat[k]/sqrt(V[k,k])
  Xr<-X[,-k,drop=FALSE];br<-solve(crossprod(Xr))%*%crossprod(Xr,y);yhat_r<-Xr%*%br;ur<-y-yhat_r;Gs<-unique(cl);set.seed(seed);tb<-numeric(B);bs<-numeric(B)
  for(b in 1:B){ wg<-sample(WEBB,length(Gs),replace=TRUE);names(wg)<-Gs;ystar<-yhat_r+wg[as.character(cl)]*ur
    bst<-XtXinv%*%crossprod(X,ystar);Vs<-cr1_vcov(X,ystar-X%*%bst,cl,XtXinv);tb[b]<-bst[k]/sqrt(Vs[k,k]);bs[b]<-bst[k] }
  if(want=="se") sd(bs) else (1+sum(abs(tb)>=abs(t_obs)))/(1+B) }
cr1_slope<-function(s){ m<-feols(any.fatalities~year*D,s,cluster=~agency.num);ct<-coeftable(m)
  if(!("year:D"%in%rownames(ct))) return(c(est=NA,p=NA)); c(est=ct["year:D","Estimate"],p=ct["year:D","Pr(>|t|)"]) }
joint_obs<-function(s,g){ m<-feols(any.fatalities~i(year,D,ref=g-1)|agency.num+year,s,cluster=~agency.num)
  w<-fixest::wald(m,"year::",print=FALSE); c(stat=unname(w$stat),p=unname(w$p)) }
wcr_joint<-function(s,g,B=499,seed=1){ W<-joint_obs(s,g)["stat"]
  mr<-feols(any.fatalities~1|agency.num+year,s);yr<-predict(mr);ur<-resid(mr);cl<-s$agency.num;Gs<-unique(cl);set.seed(seed);Wb<-numeric(B)
  for(b in 1:B){ wg<-sample(WEBB,length(Gs),replace=TRUE);names(wg)<-Gs;s$ystar<-yr+wg[as.character(cl)]*ur
    mb<-tryCatch(feols(ystar~i(year,D,ref=g-1)|agency.num+year,s,cluster=~agency.num),error=function(e)NULL)
    Wb[b]<-if(is.null(mb)) NA else tryCatch(unname(fixest::wald(mb,"year::",print=FALSE)$stat),error=function(e)NA) }
  (1+sum(Wb>=W,na.rm=TRUE))/(1+sum(!is.na(Wb))) }
# RI null over never-treated placebos (seed 2 unmatched = round 3); pools can be matched
ri_null<-function(g,n_g,B=999,seed=2,match=c("none","mean")){ match<-match.arg(match); pre<-FIRST_OUT:(g-1)
  base<-dw[dw$year%in%pre & dw$agency.id%in%never_ids(dw),];ids<-unique(base$agency.id)
  if(match=="mean"){ tids<-unique(dw$agency.id[dw$gname==g]); pm<-function(x) tapply(x$any.fatalities,x$agency.id,mean)
    tmean<-mean(pm(dw[dw$year%in%pre & dw$agency.id%in%tids,]),na.rm=TRUE); nmean<-pm(base); pool<-names(sort(abs(nmean-tmean))) }
  set.seed(seed);sl<-numeric(B);wd<-numeric(B)
  for(b in 1:B){ tr<-if(match=="mean") sample(head(pool,max(n_g*3,n_g+5)),n_g) else sample(ids,n_g)
    d2<-base;d2$D<-as.integer(d2$agency.id%in%tr)
    ml<-tryCatch(feols(any.fatalities~year*D,d2),error=function(e)NULL);sl[b]<-if(is.null(ml))NA else unname(coef(ml)["year:D"])
    mj<-tryCatch(feols(any.fatalities~i(year,D,ref=g-1)|agency.num+year,d2),error=function(e)NULL)
    wd[b]<-if(is.null(mj))NA else tryCatch(unname(fixest::wald(mj,"year::",print=FALSE)$stat),error=function(e)NA) }
  list(slope=sl,wald=wd) }

GRP<-list(`CS/Cengiz (never-treated)`=function(d,g) never_ids(d),
          `stacked (No-Change)`=function(d,g) nochange_ids(d),
          `imputation (never+not-yet)`=function(d,g) union(never_ids(d),notyet_ids(d,g)))
pretest_cohort<-function(g){ tids<-unique(dw$agency.id[dw$gname==g]); n_g<-length(tids)
  ri_u<-ri_null(g,n_g,match="none"); ri_m<-ri_null(g,n_g,match="mean",seed=4)
  do.call(rbind,lapply(names(GRP),function(gl){ s<-pre_sub(g,tids,GRP[[gl]](dw,g)); cs<-cr1_slope(s); jo<-joint_obs(s,g)
    data.frame(cohort=g,comparison=gl,slope=unname(cs["est"]),
      slope_p_CR1=unname(cs["p"]), slope_p_WCR=wcr_slope(s), slope_mde_wcr=Zc*wcr_slope(s,want="se"),
      slope_p_RIunm=(1+sum(abs(ri_u$slope)>=abs(cs["est"]),na.rm=TRUE))/(1+sum(!is.na(ri_u$slope))),
      slope_p_RImean=(1+sum(abs(ri_m$slope)>=abs(cs["est"]),na.rm=TRUE))/(1+sum(!is.na(ri_m$slope))),
      joint_p_CR1=unname(jo["p"]), joint_p_WCR=wcr_joint(s,g),
      joint_p_RIunm=(1+sum(ri_u$wald>=jo["stat"],na.rm=TRUE))/(1+sum(!is.na(ri_u$wald))),
      joint_p_RImean=(1+sum(ri_m$wald>=jo["stat"],na.rm=TRUE))/(1+sum(!is.na(ri_m$wald))), row.names=NULL) })) }

# ---- testability + untestable share + rho -----------------------------------
testable_share<-function(tab, has_comp){ # tab: cohort,weight,n_g,pre_g ; has_comp: fn(cohort)->logical
  tot<-sum(tab$weight)
  untest<-with(tab, pre_g<3 | n_g==1 | !sapply(cohort,has_comp))
  list(untestable=sum(tab$weight[untest])/tot, n_pos=sum(tab$weight>0),
       rho=suppressWarnings(cor(tab$weight,tab$pre_g,method="spearman")),
       testable_cohorts=tab$cohort[!untest]) }

# =============================================================================
# RUN
# =============================================================================
cat("computing CS-NT / CS-NYT (deck+droppers)...\n")
csnt_deck<-cs_panel(dta,"nevertreated");  csnt_drop<-cs_panel(dta_drop,"nevertreated")
csnyt_deck<-cs_panel(dta,"notyettreated"); csnyt_drop<-cs_panel(dta_drop,"notyettreated")
cat("stacked / Cengiz...\n")
st_deck<-stack_panel(build_stack(dta)); st_drop<-stack_panel(build_stack(dta_drop))
cz_deck<-stack_panel(build_cengiz(dta)); cz_drop<-stack_panel(build_cengiz(dta_drop))
cat("imputation (droppers) + bootstrap SE...\n")
imp_drop<-imp_panel(dta_drop)

# has-comparison predicates
hc_never<-function(g) TRUE                                  # never-treated always exists
hc_nyt<-function(g) length(notyet_ids(dw,g))>0              # CS-NYT needs a not-yet-treated unit
hc_stack<-function(g) g-WINDOW>=FIRST_OUT                   # stack has clean controls only for g>=2004
hc_cengiz<-function(g) TRUE
PAN<-list(
  `CS-NT|deck`=list(p=csnt_deck,hc=hc_never), `CS-NT|droppers`=list(p=csnt_drop,hc=hc_never),
  `CS-NYT|deck`=list(p=csnyt_deck,hc=hc_nyt), `CS-NYT|droppers`=list(p=csnyt_drop,hc=hc_nyt),
  `imputation|droppers`=list(p=imp_drop,hc=hc_never),
  `stacked|deck`=list(p=st_deck,hc=hc_stack), `stacked|droppers`=list(p=st_drop,hc=hc_stack),
  `Cengiz|deck`=list(p=cz_deck,hc=hc_cengiz), `Cengiz|droppers`=list(p=cz_drop,hc=hc_cengiz))

I12<-do.call(rbind,lapply(names(PAN),function(nm){ P<-PAN[[nm]]; ts<-testable_share(P$p$tab,P$hc)
  es<-strsplit(nm,"\\|")[[1]]
  data.frame(estimator=es[1],sample=es[2],ATT=P$p$att,SE=P$p$se,
    untestable_share=ts$untestable, rho_weight_preg=ts$rho, n_pos=ts$n_pos,
    n_testable=length(ts$testable_cohorts), testable_list=paste(ts$testable_cohorts,collapse=";"),
    sum_w=sum(P$p$tab$weight), recon_att=sum(P$p$tab$weight*P$p$tab$tau), row.names=NULL) }))
cat("\n=== ITEMS 1-2: ATT/SE, untestable share, rho, testable cohorts ===\n")
print(transform(I12[,c("estimator","sample","ATT","SE","untestable_share","rho_weight_preg","n_pos","n_testable","testable_list")],
  ATT=round(ATT,4),SE=round(SE,4),untestable_share=round(untestable_share,3),rho_weight_preg=round(rho_weight_preg,2)),row.names=FALSE)

cat("\n=== ITEM 3: pre-tests 2009/2013/2017 (full battery) ===\n")
PT<-do.call(rbind,lapply(c(2009,2013,2017),pretest_cohort))
print(transform(PT, slope=round(slope,4), slope_p_CR1=round(slope_p_CR1,3), slope_p_WCR=round(slope_p_WCR,3),
  slope_mde_wcr=round(slope_mde_wcr,4), slope_p_RIunm=round(slope_p_RIunm,3), slope_p_RImean=round(slope_p_RImean,3),
  joint_p_CR1=round(joint_p_CR1,3), joint_p_WCR=round(joint_p_WCR,3), joint_p_RIunm=round(joint_p_RIunm,3),
  joint_p_RImean=round(joint_p_RImean,3)),row.names=FALSE)

# ---- gap numbers (reuse committed) ------------------------------------------
gid<-read.csv(file.path(OUT,"d6_gap_identity_s1.csv")); con<-read.csv(file.path(OUT,"d6_gap_contributions.csv"))
agg_gap_imp_minus_csnyt<-sum(con$contribution); share_three<-sum(abs(con$contribution[con$cohort%in%c(2002,2017,2013)]))/sum(abs(con$contribution))

# =============================================================================
# VALIDATION
# =============================================================================
cat("\n=== VALIDATION ===\n")
cat(sprintf("[2] weights sum to 1 (max |sum_w-1|) = %.2e ; max |recon_att-ATT| = %.2e\n",
  max(abs(I12$sum_w-1)), max(abs(I12$recon_att-I12$ATT),na.rm=TRUE)))
cat(sprintf("[3] s=1 gap identity: 2001 gap=%.2e, 2002 gap=%.4f ; verify (reuse d6_gap_identity_s1)\n",
  gid$gap_0[gid$cohort==2001], gid$gap_0[gid$cohort==2002]))
rob<-read.csv(file.path(OUT,"d6_pretests_2009_2013_robust.csv")); mm<-read.csv(file.path(OUT,"d6_ri_meanmatched.csv"))
chk<-merge(PT[PT$cohort%in%c(2009,2013),], rob, by=c("cohort","comparison"), suffixes=c("",".old"))
cat("[1] reproduction vs committed (max abs diff, 2009/2013):\n")
cat(sprintf("    slope_p_CR1 %.3f WCR %.3f RIunm %.3f ; joint_p_CR1 %.3f WCR %.3f RIunm %.3f\n",
  max(abs(chk$slope_p_CR1-chk$slope_p_CR1.old)), max(abs(chk$slope_p_WCR-chk$slope_p_WCR.old)),
  max(abs(chk$slope_p_RIunm-chk$slope_p_RI)), max(abs(chk$joint_p_CR1-chk$joint_p_CR1.old)),
  max(abs(chk$joint_p_WCR-chk$joint_p_WCR.old)), max(abs(chk$joint_p_RIunm-chk$joint_p_RI))))

# =============================================================================
# CANONICAL CSV
# =============================================================================
rows<-list()
add<-function(section,key,value,note="") rows[[length(rows)+1L]]<<-data.frame(section=section,key=key,value=value,note=note,row.names=NULL)
for(i in seq_len(nrow(I12))){ r<-I12[i,]
  add("ATT",sprintf("%s|%s|ATT",r$estimator,r$sample),round(r$ATT,4))
  add("ATT",sprintf("%s|%s|SE",r$estimator,r$sample),round(r$SE,4))
  add("untestable_share",sprintf("%s|%s",r$estimator,r$sample),round(r$untestable_share,3))
  add("rho_weight_preg",sprintf("%s|%s",r$estimator,r$sample),round(r$rho_weight_preg,3),sprintf("n_pos=%d",r$n_pos))
  add("n_testable",sprintf("%s|%s",r$estimator,r$sample),r$n_testable,r$testable_list) }
for(i in seq_len(nrow(gid))) add("gap_s1",sprintf("cohort_%d",gid$cohort[i]),round(gid$gap_0[i],4))
add("gap_aggregate","imputation_minus_CSNYT",round(agg_gap_imp_minus_csnyt,4))
add("gap_aggregate","CSNYT_minus_imputation",round(-agg_gap_imp_minus_csnyt,4))
add("gap_aggregate","share_2002_2017_2013",round(share_three,3))
for(i in seq_len(nrow(con))) add("gap_contribution",sprintf("cohort_%d",con$cohort[i]),round(con$contribution[i],4),sprintf("w=%.3f gap=%.3f",con$weight[i],con$gap_avg[i]))
for(i in seq_len(nrow(PT))){ r<-PT[i,]; pfx<-sprintf("%d|%s",r$cohort,r$comparison)
  for(cc in c("slope","slope_p_CR1","slope_p_WCR","slope_mde_wcr","slope_p_RIunm","slope_p_RImean","joint_p_CR1","joint_p_WCR","joint_p_RIunm","joint_p_RImean"))
    add("pretest",paste0(pfx,"|",cc),round(r[[cc]],4)) }
CN<-do.call(rbind,rows); write.csv(CN,file.path(OUT,"d6_canonical_numbers.csv"),row.names=FALSE)
write.csv(PT,file.path(OUT,"d6_pretests_all.csv"),row.names=FALSE)
write.csv(I12,file.path(OUT,"d6_att_untestable_rho.csv"),row.names=FALSE)
cat(sprintf("\nWrote d6_canonical_numbers.csv (%d rows), d6_pretests_all.csv, d6_att_untestable_rho.csv\n",nrow(CN)))
cat(sprintf("aggregate gap: imputation-CSNYT=%.4f, CSNYT-imputation=%.4f, 2002+2017+2013 share=%.1f%%\n",
  agg_gap_imp_minus_csnyt,-agg_gap_imp_minus_csnyt,100*share_three))
