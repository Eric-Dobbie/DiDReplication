# =============================================================================
# d7_geometry_sweep.R  --  analytic geometry sweep (paper section 3.3 / appendix).
# DETERMINISTIC: weights are closed-form functions of geometry, MDE is the iid
# closed form. No outcomes, no estimation, no random draws. Scope is the fixed
# 12-geometry grid below -- nothing added.
#
# Setup: t=1..T, sigma=1, N_c=40 never-treated, N_T=40 treated over K=4 cohorts,
# each unit a cluster, absorbing treatment, no covariates.
# =============================================================================
suppressMessages({library(ggplot2)})
OUT<-Sys.getenv("DIDREP_OUT",unset="output"); FIG<-Sys.getenv("DIDREP_FIG",unset="figures"); DATA<-Sys.getenv("DIDREP_DATA",unset="data")
Zc<-1.96+0.84                               # as written in the brief (approx z_.975+z_.80)
N_C<-40; N_T<-40; K<-4; THETA<-c(0.10,0.25,0.50)

# ---- adoption times: round(seq(a,b,len=4)), dedup by shifting up -------------
adopt_times<-function(T,timing){ half<-ceiling(T/2)
  win<-switch(timing, early=c(2,half), spread=c(2,T), late=c(half,T))
  g<-round(seq(win[1],win[2],length.out=4))
  for(i in 2:4) if(g[i]<=g[i-1]) g[i]<-g[i-1]+1      # shift up until 4 distinct
  g }
cohort_sizes<-function(size){ if(size=="equal") rep(N_T/K,K) else { v<-rep(5L,K); v[2]<-25L; v } }  # concentrated: 2nd-earliest=25

# ---- closed-form pieces -----------------------------------------------------
mde_fn<-function(pre,n_g,n_c){ if(pre<2 || n_g<2 || n_c==0) return(Inf)
  Zc*sqrt( (1/n_g + 1/n_c) / (pre*(pre^2-1)/12) ) }
stacked_w<-function(g,n,T){ Ns<-n+N_C; p<-n/Ns; q<-(T-g+1)/T; Vs<-Ns*T*p*(1-p)*q*(1-q); if(sum(Vs)>0) Vs/sum(Vs) else rep(NA_real_,length(Vs)) }

# ---- one geometry -> per-cohort rows for CS / imputation / stacked -----------
geom_rows<-function(T,timing,size){
  g<-adopt_times(T,timing); n<-cohort_sizes(size)
  pre<-g-1; post<-T-g+1
  W_cs<-(n*post)/sum(n*post)                              # CS-NT simple (= imputation)
  w_s<-stacked_w(g,n,T)
  nc_imp<-sapply(seq_along(g),function(i) N_C + sum(n[g>g[i]]))   # initial control group
  mk<-function(est,W,n_c_vec,a_vec){ do.call(rbind,lapply(seq_along(g),function(i){
    data.frame(T=T,timing=timing,size=size,geom=sprintf("T%d|%s|%s",T,timing,size),
      estimator=est,cohort_time=g[i],pre_g=pre[i],post_g=post[i],n_g=n[i],n_c=n_c_vec[i],
      weight=W[i],MDE=mde_fn(pre[i],n[i],n_c_vec[i]),a_g=a_vec[i], row.names=NULL) })) }
  rbind(
    mk("CS",         W_cs, rep(N_C,K), (post+1)/2),
    mk("imputation", W_cs, nc_imp,     rep(T/2,K)),
    mk("stacked",    w_s,  rep(N_C,K), rep(T/2,K)) ) }

# ---- sweep the 12 geometries ------------------------------------------------
GRID<-expand.grid(T=c(10,25), timing=c("early","spread","late"), size=c("equal","concentrated"), stringsAsFactors=FALSE)
COH<-do.call(rbind,lapply(seq_len(nrow(GRID)),function(i) geom_rows(GRID$T[i],GRID$timing[i],GRID$size[i])))
for(th in THETA) COH[[sprintf("deltabar_%.2f",th)]]<-th/COH$a_g
write.csv(COH,file.path(OUT,"d7_sweep_cohort.csv"),row.names=FALSE)

# report the final adoption times per geometry
cat("=== adoption times per geometry (after dedup/shift) ===\n")
at<-unique(COH[,c("T","timing","size","geom")]); at<-at[at$size=="equal",]  # timing/T independent of size
for(i in seq_len(nrow(at))){ g<-unique(COH$cohort_time[COH$T==at$T[i]&COH$timing==at$timing[i]])
  cat(sprintf("  T=%2d %-6s : %s\n", at$T[i], at$timing[i], paste(sort(g),collapse=","))) }

# ---- summary per geometry x estimator ---------------------------------------
spear<-function(w,pg) if(length(w)>=3 && sd(w)>0) suppressWarnings(cor(w,pg,method="spearman")) else NA_real_
SUM<-do.call(rbind,lapply(split(COH,list(COH$geom,COH$estimator),drop=TRUE),function(d){
  r<-data.frame(T=d$T[1],timing=d$timing[1],size=d$size[1],geom=d$geom[1],estimator=d$estimator[1],
    U_inf=sum(d$weight[is.infinite(d$MDE)]), rho=spear(d$weight,d$pre_g), row.names=NULL)
  for(th in THETA){ db<-th/d$a_g; r[[sprintf("U_%.2f",th)]]<-sum(d$weight[d$MDE>db]) }; r }))
SUM<-SUM[order(SUM$timing,SUM$T,SUM$size,SUM$estimator),]
write.csv(SUM,file.path(OUT,"d7_sweep_summary.csv"),row.names=FALSE)
options(width=200); cat("\n=== summary (U(inf), U at break-even, rho) ===\n")
print(transform(SUM, U_inf=round(U_inf,3), U_0.10=round(U_0.10,3), U_0.25=round(U_0.25,3), U_0.50=round(U_0.50,3), rho=round(rho,2)),row.names=FALSE)

# ---- VALIDATION 2 & 3 -------------------------------------------------------
sw<-tapply(COH$weight,list(COH$geom,COH$estimator),sum)
cat(sprintf("\n[V2] weights sum to 1 within geometry x estimator: max|.-1| = %.2e\n", max(abs(sw-1),na.rm=TRUE)))
mono<-all(apply(SUM[,c("U_0.10","U_0.25","U_0.50")],1,function(u) all(diff(u)<=1e-12)))
ge<-all(SUM$U_0.50>=SUM$U_inf-1e-12 & SUM$U_0.10>=SUM$U_inf-1e-12)
cat(sprintf("[V3] U non-increasing in theta*: %s ; U(break-even) >= U(inf): %s\n",mono,ge))

# =============================================================================
# VALIDATION 1: P&P geometry through the same code
# =============================================================================
cat("\n=== [V1] P&P validation ===\n")
# CS: droppers cohorts/sizes, T=21, N_c=41
tab<-read.csv(file.path(OUT,"weight_checkability_table.csv")); cs<-tab[tab$estimator=="CS"&tab$sample=="droppers",]
cs<-cs[order(cs$cohort),]; W_cs_pp<-(cs$n_g*cs$post_g)/sum(cs$n_g*cs$post_g)
cat(sprintf("CS weights vs canonical CS-NT simple (droppers): max|diff| = %.2e\n", max(abs(W_cs_pp-cs$weight))))
# Stacked: fwl canonical w_s vs closed form with the stack's actual control count.
# Use the ESTIMATION window (observed outcome, 2000-2020, T=21) -- the sample the FWL w_s
# were computed on -- not the full 1987-2020 stack panel.
fw<-read.csv(file.path(OUT,"fwl_decomp_unweighted.csv")); st0<-read.csv(file.path(DATA,"stacked_fatal.csv"))
st<-st0[!is.na(st0$any.fatalities),]; Tpp<-length(unique(st$year))   # 21
sg<-sort(unique(st$cohort))
ng_st<-sapply(sg,function(g) length(unique(st$agency.id[st$cohort==g & st$treat==1])))
nc_st<-sapply(sg,function(g) length(unique(st$agency.id[st$cohort==g & st$treat==0])))
post_st<-sapply(sg,function(g) length(unique(st$year[st$year>=g])))   # observed post periods
Vs_cf<-(ng_st+nc_st)*Tpp*(ng_st/(ng_st+nc_st))*(1-ng_st/(ng_st+nc_st))*(post_st/Tpp)*(1-post_st/Tpp)
Vs_cf[nc_st==0]<-0                                        # no clean controls -> closed form gives 0
w_cf<-Vs_cf/sum(Vs_cf)
cmp<-merge(data.frame(stack=sg,w_cf=w_cf,nc=nc_st), fw[,c("stack","w_s")], by="stack")
cmp$diff<-cmp$w_cf-cmp$w_s
cat(sprintf("Stacked closed-form w_s vs canonical fwl w_s: correlation %.3f, L1 %.3f, max|diff| %.3f\n",
  cor(cmp$w_cf,cmp$w_s), sum(abs(cmp$diff)), max(abs(cmp$diff))))
cat("  largest gaps:\n"); print(transform(cmp[order(-abs(cmp$diff)),][1:4,c("stack","nc","w_cf","w_s","diff")], w_cf=round(w_cf,3),w_s=round(w_s,3),diff=round(diff,3)),row.names=FALSE)
cat("  attribution: (i) 2002: CONTROLS DEFINITION -- closed form uses nc=0 -> V_s=0, canonical gets V_s>0 from the within-treated adopter contrast;\n")
cat("               (ii) 2016 (and other single-treated reversals): the real stacked fit gives V_s=0 (collinear/no within-stack contrast), closed form gives nonzero;\n")
cat("               (iii) remaining gaps (2009, 2004): IMBALANCE -- closed form assumes a balanced common-pool stack, but the real stack mixes 41 never + ~700 always-treated controls.\n")

# =============================================================================
# FIGURE: U(inf) and U at theta*=0.25, 12 geometries grouped by timing
# =============================================================================
library(tidyr)
plt<-SUM[,c("T","timing","size","estimator","U_inf","U_0.25")]
plt<-pivot_longer(plt,cols=c("U_inf","U_0.25"),names_to="Utype",values_to="U")
plt$Utype<-factor(plt$Utype,levels=c("U_inf","U_0.25"),labels=c("U(infinity): testability only","U at break-even theta*=0.25"))
plt$timing<-factor(plt$timing,levels=c("early","spread","late"))
plt$xg<-factor(paste0("T",plt$T,"\n",ifelse(plt$size=="equal","eq","conc")),
  levels=c("T10\neq","T10\nconc","T25\neq","T25\nconc"))
plt$estimator<-factor(plt$estimator,levels=c("CS","imputation","stacked"))
g<-ggplot(plt,aes(xg,U,color=estimator,shape=estimator,group=estimator))+
  geom_point(size=3,position=position_dodge(width=0.5))+
  facet_grid(Utype~timing)+ ylim(0,1)+
  scale_color_manual(values=c(CS="#0072B2",imputation="#009E73",stacked="#D55E00"))+
  scale_shape_manual(values=c(CS=16,imputation=17,stacked=15))+
  labs(x="geometry (grouped by timing; T and size profile)",y="unchecked weight U",
    title="D7 analytic geometry sweep: unchecked weight by estimator",
    subtitle="rows: U(infinity) = untestable weight vs U at break-even theta*=0.25 ; cols: adoption timing")+
  theme_minimal(base_size=11)+theme(panel.grid.minor=element_blank(),legend.position="top",
    panel.spacing=unit(0.8,"lines"),axis.text.x=element_text(size=7))
ggsave(file.path(FIG,"fig_d7_geometry_sweep.pdf"),g,width=11,height=6,device=cairo_pdf)
ggsave(file.path(FIG,"fig_d7_geometry_sweep.png"),g,width=11,height=6,dpi=200)
file.copy(file.path(FIG,"fig_d7_geometry_sweep.png"),file.path(OUT,"fig_d7_geometry_sweep.png"),overwrite=TRUE)

# ---- CS-vs-stacked discrimination (for the writeup) -------------------------
cat("\n=== CS vs stacked gap in U (for RESULTS) ===\n")
w<-reshape(SUM[,c("geom","timing","T","size","estimator","U_inf","U_0.25")],idvar=c("geom","timing","T","size"),timevar="estimator",direction="wide")
w$gap_Uinf<-w$U_inf.CS-w$U_inf.stacked; w$gap_U25<-w$U_0.25.CS-w$U_0.25.stacked
print(transform(w[order(-w$gap_Uinf),c("geom","U_inf.CS","U_inf.stacked","gap_Uinf","U_0.25.CS","U_0.25.stacked","gap_U25")],
  U_inf.CS=round(U_inf.CS,2),U_inf.stacked=round(U_inf.stacked,2),gap_Uinf=round(gap_Uinf,2),
  U_0.25.CS=round(U_0.25.CS,2),U_0.25.stacked=round(U_0.25.stacked,2),gap_U25=round(gap_U25,2)),row.names=FALSE)

cat("\nWrote d7_sweep_cohort.csv, d7_sweep_summary.csv, fig_d7_geometry_sweep.*\n")
