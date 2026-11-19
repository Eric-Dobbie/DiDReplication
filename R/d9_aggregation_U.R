# =============================================================================
# d9_aggregation_U.R  --  unchecked-weight diagnostic U^k(theta*) under every
# aggregation choice (paper section 3.2 extension). Droppers sample.
#
# No new estimation. Per-cohort geometry (cohort, n_g, pre_g, post_g) and MDEs
# are the round-6 values in d6_round6_mde.csv; reference rows (stacked, Cengiz)
# are the round-6 single-weighting values in d6_unchecked_weight.csv. Aggregation
# weights are built from geometry to match did::aggte with default settings.
#
#   For estimator k, aggregation A, over post cells (g,t), t>=g:
#     W_g        = sum_{t>=g} w^A_{g,t}
#     a^{k,A}_g  = (1/W_g) sum_{t>=g} w^A_{g,t} (t - b^k_g)
#                  b_g = g-1 (CS-NT, CS-NYT, SA); mean(2000..g-1) (imputation)
#     deltabar   = theta* / a^{k,A}_g
#     U          = sum_g W_g 1{ MDE^k_g > deltabar }      (Inf MDE always counts)
#     U(inf)     = sum_g W_g 1{ MDE^k_g = Inf }
#
#   Aggregations (did::aggte defaults; cohorts weighted by P(G=g)):
#     simple      w_{g,t} ∝ n_g on every post cell          -> W_g ∝ n_g*post_g
#     cohort      P(G=g)/post_g per post cell               -> W_g  = P(G=g)
#     eventstudy  (1/E) P(G=g | G+e<=T) on cell (g,g+e)      overall dynamic
#     calendar    (1/C) P(G=g | G<=t)  on cell (g,t)         overall calendar
#
# SA uses never-treated controls => shares CS-NT geometry AND CS-NT MDEs.
# Outputs: tables/table_aggregation_U.tex, output/d9_aggregation_U.csv,
#          output/d9_eventstudy_by_horizon.csv
# =============================================================================
OUT<-Sys.getenv("DIDREP_OUT",unset="output"); TAB<-Sys.getenv("DIDREP_TAB",unset="tables")
dir.create(TAB,showWarnings=FALSE,recursive=TRUE)
T_CAL<-2020; THETA<-c(0.05,0.10)

md<-read.csv(file.path(OUT,"d6_round6_mde.csv"),stringsAsFactors=FALSE)
md<-md[md$sample=="droppers",c("estimator","cohort","n_g","pre_g","post_g","MDE")]
md$MDE<-suppressWarnings(as.numeric(md$MDE))        # "Inf" -> Inf, "NA" -> NA
# geometry list per estimator (CS-NT, CS-NYT, imputation present; SA := CS-NT)
ests_in<-unique(md$estimator)
geom<-setNames(lapply(ests_in,function(e){ d<-md[md$estimator==e,]; d[order(d$cohort),] }),ests_in)
geom[["SA"]]<-geom[["CS-NT"]]                        # never-treated controls

baseline<-function(est,g) if(est=="imputation") (2000+(g-1))/2 else g-1

# ---- per-cohort (W_g, a_g) for one estimator x aggregation ------------------
build<-function(est,agg){
  d<-geom[[est]]; g<-d$cohort; n<-d$n_g; MDE<-d$MDE; K<-length(g); b<-baseline(est,g)
  W<-numeric(K); a<-numeric(K)
  if(agg%in%c("simple","cohort")){
    tot<-if(agg=="simple") sum(n*d$post_g) else sum(n)
    for(i in seq_len(K)){ ts<-g[i]:T_CAL
      W[i]<-if(agg=="simple") n[i]*d$post_g[i]/tot else n[i]/tot
      a[i]<-mean(ts-b[i]) }                          # equal within cohort
  } else if(agg=="eventstudy"){
    gmin<-min(g); E<-T_CAL-gmin+1
    for(i in seq_len(K)){ es<-0:(T_CAL-g[i])
      wi<-sapply(es,function(e){ den<-sum(n[g+e<=T_CAL]); (1/E)*n[i]/den })
      W[i]<-sum(wi); a[i]<-sum(wi*((g[i]+es)-b[i]))/W[i] }
  } else {                                           # calendar
    gmin<-min(g); C<-T_CAL-gmin+1
    for(i in seq_len(K)){ ts<-g[i]:T_CAL
      wi<-sapply(ts,function(t){ den<-sum(n[g<=t]); (1/C)*n[i]/den })
      W[i]<-sum(wi); a[i]<-sum(wi*(ts-b[i]))/W[i] }
  }
  data.frame(cohort=g,n_g=n,W=W,a=a,MDE=MDE)
}
U_of<-function(p){
  Uinf<-sum(p$W[is.infinite(p$MDE)])
  out<-c(U_inf=Uinf,sumW=sum(p$W))
  for(th in THETA){ db<-th/p$a; checked<-is.finite(p$MDE)&(p$MDE<=db)
    out[[sprintf("U_%.2f",th)]]<-sum(p$W[!checked]) }
  out
}

ESTS<-c("CS-NT","CS-NYT","SA","imputation")
AGGS<-c("simple","cohort","eventstudy","calendar")
full<-do.call(rbind,lapply(ESTS,function(e) do.call(rbind,lapply(AGGS,function(a){
  u<-U_of(build(e,a)); data.frame(estimator=e,aggregation=a,U_inf=u[["U_inf"]],
    U_0.05=u[["U_0.05"]],U_0.10=u[["U_0.10"]],sumW=u[["sumW"]]) }))))

# reference rows from d6_unchecked_weight.csv (droppers)
uw<-read.csv(file.path(OUT,"d6_unchecked_weight.csv"),stringsAsFactors=FALSE)
refv<-function(kind,est,th=NA){ r<-uw[uw$kind==kind&uw$estimator==est&uw$sample=="droppers",]
  if(is.na(th)) as.numeric(r$value[r$theta_star=="NA"|is.na(r$theta_star)][1])
  else as.numeric(r$value[r$theta_star==th][1]) }
for(e in c("stacked","Cengiz")) full<-rbind(full,data.frame(estimator=e,aggregation="(single implied)",
  U_inf=refv("U_inf",e),U_0.05=refv("U_breakeven",e,"0.05"),U_0.10=refv("U_breakeven",e,"0.1"),sumW=NA))
write.csv(full,file.path(OUT,"d9_aggregation_U.csv"),row.names=FALSE)

# ---- per-horizon event study ------------------------------------------------
EH<-list()
for(est in ESTS){ d<-geom[[est]]; g<-d$cohort; n<-d$n_g; MDE<-d$MDE; gmin<-min(g); E<-T_CAL-gmin+1
  for(e in 0:(E-1)){ idx<-which(g+e<=T_CAL); den<-sum(n[idx]); b<-baseline(est,g[idx])
    a<-(g[idx]+e)-b; W<-n[idx]/den; m<-MDE[idx]
    Ue_inf<-sum(W[is.infinite(m)])
    Ue05<-sum(W[is.infinite(m)|m>THETA[1]/a]); Ue10<-sum(W[is.infinite(m)|m>THETA[2]/a])
    EH[[length(EH)+1L]]<-data.frame(estimator=est,horizon_e=e,cohort=g[idx],n_g=n[idx],
      W_ge=W,a_ge=a,MDE=m,untestable=as.integer(is.infinite(m)),
      delta_05=THETA[1]/a,delta_10=THETA[2]/a,
      checked_05=as.integer(is.finite(m)&m<=THETA[1]/a),checked_10=as.integer(is.finite(m)&m<=THETA[2]/a),
      Ue_inf=Ue_inf,Ue_05=Ue05,Ue_10=Ue10,sumW_e=sum(W)) } }
write.csv(do.call(rbind,EH),file.path(OUT,"d9_eventstudy_by_horizon.csv"),row.names=FALSE)

# ---- LaTeX table (panel sub-headers; no \multirow dependency) ---------------
c3<-function(x) sprintf("%.3f",x)
estlab<-c(`CS-NT`="Callaway--Sant'Anna, never-treated controls",`CS-NYT`="Callaway--Sant'Anna, not-yet-treated controls",
          SA="Sun--Abraham, never-treated controls",imputation="Imputation (BJS/Gardner)")
agglab<-c(simple="Simple",cohort="Cohort (group)",eventstudy="Event study, overall",calendar="Calendar, overall")
L<-c("\\begin{tabular}{lccc}","\\toprule",
     "Aggregation & $U(\\infty)$ & $U(\\theta^{*}{=}0.05)$ & $U(\\theta^{*}{=}0.10)$ \\\\","\\midrule")
for(i in seq_along(ESTS)){ est<-ESTS[i]
  L<-c(L,sprintf("\\multicolumn{4}{@{}l}{\\emph{%s}} \\\\",estlab[[est]]))
  for(a in AGGS){ r<-full[full$estimator==est&full$aggregation==a,]
    L<-c(L,sprintf("\\quad %s & %s & %s & %s \\\\",agglab[[a]],c3(r$U_inf),c3(r$U_0.05),c3(r$U_0.10))) }
  if(i<length(ESTS)) L<-c(L,"\\addlinespace") }
L<-c(L,"\\midrule","\\multicolumn{4}{@{}l}{\\emph{Reference rows (single implied weighting)}} \\\\")
for(pr in list(c("stacked","Stacked, full panel (agency $\\times$ stack FE)"),
               c("Cengiz","Cengiz $\\pm 4$, never-treated controls"))){
  r<-full[full$estimator==pr[1],]; L<-c(L,sprintf("\\quad %s & %s & %s & %s \\\\",pr[2],c3(r$U_inf),c3(r$U_0.05),c3(r$U_0.10))) }
L<-c(L,"\\bottomrule","\\end{tabular}","")
notes<-paste0("\\noindent Droppers sample (\\S4.3--4.5 analysis sample). $U$ is the aggregation weight on cohorts whose ",
"pre-trend slope test cannot detect a trend large enough to bias the estimand by $\\theta^{*}$ ",
"($\\mathrm{MDE}_g>\\theta^{*}/a_g$); $U(\\infty)$ is the weight on cohorts with no feasible pre-test ",
"($\\mathrm{MDE}_g=\\infty$). Per-cohort MDEs are the round-6 WCR values (\\texttt{d6\\_round6\\_mde.csv}) ",
"for each estimator's comparison group. Sun--Abraham uses never-treated controls and therefore shares ",
"Callaway--Sant'Anna (never-treated)'s geometry and MDEs. Baseline period $b_g$: $g-1$ for ",
"Callaway--Sant'Anna, Sun--Abraham; the mean of cohort $g$'s pre-periods for imputation. ",
"Weights reproduce \\texttt{did::aggte} with default settings (simple $\\propto n_g$ per post cell; ",
"group $= P(G{=}g)$; event-study and calendar equally weight their horizons/periods). Reference rows ",
"carry the single implied weighting of each design (round-6 values).")
L<-c(L,"\\smallskip",paste0("{\\footnotesize ",notes,"\\par}"))
writeLines(L,file.path(TAB,"table_aggregation_U.tex"))

# ---- validation -------------------------------------------------------------
num<-full[!is.na(full$sumW),]
cat(sprintf("[V1] simple CS-NT U_inf=%.6f imputation U_inf=%.6f (target 0.243243)\n",
  full$U_inf[full$estimator=="CS-NT"&full$aggregation=="simple"],
  full$U_inf[full$estimator=="imputation"&full$aggregation=="simple"]))
cat(sprintf("[V3] max|sumW-1| = %.2e\n",max(abs(num$sumW-1))))
cat(sprintf("[V4] U(theta*)>=U(inf) & non-increasing: %s\n",
  all(num$U_0.05>=num$U_inf-1e-9 & num$U_0.10>=num$U_inf-1e-9 & num$U_0.05>=num$U_0.10-1e-9)))
cat("Wrote tables/table_aggregation_U.tex, output/d9_aggregation_U.csv, output/d9_eventstudy_by_horizon.csv\n")
