# =============================================================================
# taskB3B6_validation.R
#   B6: what does n (used in the /n correction) count, and which agencies are in
#       the CS estimation sample?
#   B3: validate the hand-rolled Wald vs did's own pooled pre-test Wpval, and
#       confirm V_analytical/n is a UNIFORM rescaling (off-diagonals too) by
#       checking V_analytical == crossprod(inffunc)/n and comparing the analytical
#       correlation matrix to the multiplier-bootstrap correlation matrix.
# =============================================================================
suppressMessages({library(did); library(MASS)})
DATA_DIR <- Sys.getenv("DIDREP_DATA", unset="data"); OUT_DIR <- Sys.getenv("DIDREP_OUT", unset="output")
dta <- read.csv(file.path(DATA_DIR,"dta.csv")); dta$agency.num <- as.numeric(factor(dta$agency.id))
dta$year.changed <- ave(dta$year.changed, dta$agency.id, FUN=function(x){v<-x[!is.na(x)]; if(length(v)) min(v) else NA_real_})
dta$yc <- ifelse(is.na(dta$year.changed)&dta$no.req==0,0, ifelse(is.na(dta$year.changed)&dta$no.req==1,1987,dta$year.changed))

set.seed(0)
mp <- did::att_gt(yname="any.fatalities",tname="year",idname="agency.num",gname="yc",
                  xformla=~1,control_group="nevertreated",clustervars="agency.num",base_period="varying",data=dta)

cat("=========== B6: what is n ===========\n")
n <- mp$DIDparams$n
cat("mp$DIDparams$n =", n, "\n")
idata <- mp$DIDparams$data; uids <- unique(idata[[mp$DIDparams$idname]])
cat("unique units in att_gt estimation data =", length(uids), "\n")
# map back to gname composition
ug <- idata[[mp$DIDparams$gname]][!duplicated(idata[[mp$DIDparams$idname]])]
cat("estimation-sample units by group (0 = never-treated controls):\n")
print(table(ug))
cat(sprintf("=> n=%d = %d never-treated controls + %d treated (cohorts 2001-2020)\n",
            n, sum(ug==0), sum(ug>0)))
cat(sprintf("total agencies in panel = %d; dropped by att_gt = %d (716 already-treated: 700 always-no-req coded 1987 + pre-2000 changers)\n",
            length(unique(dta$agency.id)), length(unique(dta$agency.id)) - n))

cat("\n=========== B3: Wpval + machinery validation ===========\n")
cat("inference settings: bstrap =", mp$DIDparams$bstrap, " cband =", mp$DIDparams$cband,
    " biters =", mp$DIDparams$biters, "\n")
cat("did pooled pre-test Wpval =", ifelse(is.null(mp$Wpval), "NULL", as.character(mp$Wpval)), "\n")

# hand-rolled pooled Wald on ALL pre-treatment coefficients, V_analytical/n
att<-mp$att; G<-mp$group; Tt<-mp$t; V <- as.matrix(mp$V_analytical)/n
idx <- which(Tt < G & !is.na(att))
Vb <- V[idx,idx]; a <- att[idx]
sv <- svd(Vb)$d; rk <- sum(sv > sqrt(.Machine$double.eps)*max(sv))
cat(sprintf("\npooled pre-cells: %d coefficients, vcov rank = %d (rank-deficient: %s)\n",
            length(idx), rk, rk < length(idx)))
if (rk < length(idx)) {
  cat("=> pooled analytical pre-test vcov is SINGULAR (matches did's 'singular covariance' warning);\n")
  cat("   did returns Wpval = NULL for this reason. A full-rank pooled Wald is not computable either.\n")
}

# IF identity: V_analytical == crossprod(inffunc)/n  (uniform scaling, off-diagonals included)
IF <- mp$inffunc
Vif <- as.matrix(Matrix::crossprod(IF)) / n
Va  <- as.matrix(mp$V_analytical)
cat(sprintf("\nuniform-scaling check: max|V_analytical - crossprod(inffunc)/n| = %.3e (0 => V_analytical IS crossprod(IF)/n, so /n rescales the WHOLE matrix)\n",
            max(abs(Va - Vif))))

# analytical vs bootstrap correlation for a well-identified cohort block (2009)
blk <- which(G==2009 & Tt<2009 & !is.na(att))
Ca <- cov2cor(Va[blk,blk])
set.seed(1); Bm <- 2000; k <- length(blk)
IFb <- as.matrix(IF[,blk])
boot <- matrix(NA_real_, Bm, k)
for (b in 1:Bm) { w <- sample(c(-1,1), nrow(IFb), replace=TRUE); boot[b,] <- colSums(IFb*w)/n }
Cb <- cor(boot)
cat(sprintf("2009 pre-block: max|corr_analytical - corr_bootstrap| = %.3f (small => off-diagonal structure agrees)\n",
            max(abs(Ca - Cb))))
cat("analytical corr (upper 3x3):\n"); print(round(Ca[1:3,1:3],3))
cat("bootstrap  corr (upper 3x3):\n"); print(round(Cb[1:3,1:3],3))

sink(file.path(OUT_DIR,"taskB3B6_validation.txt"))
cat("n =", n, "=", sum(ug==0), "never-treated +", sum(ug>0), "treated\n")
cat("Wpval =", ifelse(is.null(mp$Wpval),"NULL",as.character(mp$Wpval)), "\n")
cat("pooled pre vcov rank", rk, "of", length(idx), "\n")
cat("max|V - crossprod(IF)/n| =", max(abs(Va-Vif)), "\n")
cat("2009 block max|corr_anal - corr_boot| =", max(abs(Ca-Cb)), "\n")
sink()
cat("\nWrote output/taskB3B6_validation.txt\n")
