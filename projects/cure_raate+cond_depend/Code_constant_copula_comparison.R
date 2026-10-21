library(QRM)
library(numDeriv)
library(parallel)

## ============================================================
## CONSTANT-COPULA COMPARISON
## n = 500, 20% susceptible censoring, 3 dependence settings
## Uses the same DGP and replication seeds as the main experiment.
## ============================================================
R <- 1000L
dc <- parallel::detectCores()
cores <- if(is.na(dc)) 2L else max(1L, dc - 1L)

n <- 1000L
target.censor <- .20
prob.trt <- c(.30,.30,.40)       # Trt 1, Trt 2, reference
p.true <- .70
a.true <- c(0,2,1)
b.true <- c(0,2,1)
sigma1.true <- sigma2.true <- 1

tau.scenarios <- list(
  null     = c(.40,.40,.40),
  moderate = c(.40,.45,.60),
  strong   = c(.40,.55,.70)
)

tau2gamma <- function(tau){
  phi <- 2*tau/(1-tau)
  c(gamma0=log(phi[1]),
    gamma1=log(phi[2]/phi[1]),
    gamma2=log(phi[3]/phi[1]))
}

## IDs 2, 6, 10 are exactly the n=500 / 20% rows
## from the original 12-scenario script, preserving its seeds.
scenarios <- data.frame(
  id=c(2L,6L,10L),
  dependence=c("null","moderate","strong"),
  stringsAsFactors=FALSE
)
G <- t(vapply(tau.scenarios, tau2gamma, numeric(3)))
scenarios <- cbind(scenarios, G[scenarios$dependence,,drop=FALSE])

## ============================================================
## DISTRIBUTION / COPULA FUNCTIONS
## ============================================================
margin <- function(t,mu,sigma){
  z <- (log(t)-mu)/sigma
  F <- pmin(pmax(QRM::pGumbel(z,0,1),1e-12),1-1e-12)
  logf <- QRM::dGumbel(z,0,1,log=TRUE)-log(sigma)-log(t)
  list(F=F,logf=logf)
}

log1pexp <- function(x){
  i <- x>35; out <- numeric(length(x))
  out[i] <- x[i]+log1p(exp(-x[i]))
  out[!i] <- log1p(exp(x[!i]))
  out
}
logexpm1 <- function(x){
  x <- pmax(x,1e-14); i <- x>35; out <- numeric(length(x))
  out[i] <- x[i]+log1p(-exp(-x[i]))
  out[!i] <- log(expm1(x[!i]))
  out
}

rclayton <- function(phi){
  u <- pmin(pmax(runif(length(phi)),1e-12),1-1e-12)
  w <- pmin(pmax(runif(length(phi)),1e-12),1-1e-12)
  lv <- -log1pexp(-phi*log(u)+
                    logexpm1(-phi/(1+phi)*log(w)))/phi
  cbind(u,pmin(pmax(exp(lv),1e-12),1-1e-12))
}

clayton <- function(u,v,phi){
  u <- pmin(pmax(u,1e-12),1-1e-12)
  v <- pmin(pmax(v,1e-12),1-1e-12)
  lu <- log(u); lv <- log(v)
  aa <- -phi*lu; bb <- -phi*lv; mm <- pmax(aa,bb)
  logA <- mm+log(exp(aa-mm)+exp(bb-mm)-exp(-mm))
  C1 <- exp(pmin((-phi-1)*lu+(-1/phi-1)*logA,0))
  C1 <- pmin(pmax(C1,1e-12),1-1e-12)
  logc <- log1p(phi)+(-1-phi)*(lu+lv)+(-2-1/phi)*logA
  list(C1=C1,logc=logc)
}

## ============================================================
## SAME EMPIRICAL U_C CALIBRATION AS MAIN EXPERIMENT
## ============================================================
set.seed(2026)
B.cal <- 300000L
trt.cal <- sample(1:3,B.cal,TRUE,prob.trt)
x1.cal <- as.integer(trt.cal==1)
x2.cal <- as.integer(trt.cal==2)
mu1.cal <- a.true[1]+a.true[2]*x1.cal+a.true[3]*x2.cal
T1.cal <- exp(mu1.cal+
                sigma1.true*QRM::qGumbel(runif(B.cal),0,1))

f.uc <- function(uc) mean(pmin(T1.cal/uc,1))-target.censor
hi <- 10
while(f.uc(hi)>0) hi <- 2*hi
uc <- uniroot(f.uc,c(1e-8,hi),tol=1e-8)$root
cat("Calibrated u_c =",uc,"\n")

## ============================================================
## DATA GENERATOR: IDENTICAL TO MAIN EXPERIMENT
## ============================================================
simulate.data <- function(gamma){
  trt <- sample(1:3,n,TRUE,prob.trt)
  x1 <- as.integer(trt==1)
  x2 <- as.integer(trt==2)

  mu1 <- a.true[1]+a.true[2]*x1+a.true[3]*x2
  mu2 <- b.true[1]+b.true[2]*x1+b.true[3]*x2
  phi <- exp(gamma[1]+gamma[2]*x1+gamma[3]*x2)

  UV <- rclayton(phi)
  T1 <- exp(mu1+sigma1.true*QRM::qGumbel(UV[,1],0,1))
  T2 <- exp(mu2+sigma2.true*QRM::qGumbel(UV[,2],0,1))

  susceptible <- rbinom(n,1,p.true)
  T1.obs <- ifelse(susceptible==1,T1,Inf)
  C <- runif(n,0,uc)

  d1 <- as.integer(T1.obs<=C)
  y1 <- pmin(T1.obs,C)
  d2 <- numeric(n); y2 <- rep(NA_real_,n)

  id1 <- which(d1==1)
  if(length(id1)){
    rem <- C[id1]-T1[id1]
    d2[id1] <- as.integer(T2[id1]<=rem)
    y2[id1] <- pmin(T2[id1],rem)
  }

  dat <- data.frame(x1,x2,y1,y2,d1,d2)
  attr(dat,"rates") <- c(
    censor1.sus=mean(C[susceptible==1]<T1[susceptible==1]),
    noevent1.all=mean(d1==0),
    censor2=if(length(id1)) mean(d2[id1]==0) else NA_real_
  )
  dat
}

## ============================================================
## RESTRICTED FIT: phi_i = exp(gamma0)
## ============================================================
par.names <- c("p","a0","a1","a2","sigma1",
               "b0","b1","b2","sigma2","gamma0")

natural.par <- function(z)
  setNames(c(plogis(z[1]),z[2:4],exp(z[5]),
             z[6:8],exp(z[9]),z[10]),par.names)

nll.i <- function(z,dat){
  p <- plogis(z[1])
  a <- z[2:4]; s1 <- exp(z[5])
  b <- z[6:8]; s2 <- exp(z[9])
  phi <- exp(z[10])

  M1 <- margin(dat$y1,a[1]+a[2]*dat$x1+a[3]*dat$x2,s1)
  ans <- numeric(nrow(dat))

  id0 <- which(dat$d1==0)
  ans[id0] <- -log1p(-p*M1$F[id0])

  id1 <- which(dat$d1==1)
  if(length(id1)){
    M2 <- margin(dat$y2[id1],
                 b[1]+b[2]*dat$x1[id1]+b[3]*dat$x2[id1],s2)
    CP <- clayton(M1$F[id1],M2$F,phi)

    j0 <- which(dat$d2[id1]==0)
    if(length(j0))
      ans[id1[j0]] <- -(log(p)+M1$logf[id1[j0]]+
                          log1p(-CP$C1[j0]))

    j1 <- which(dat$d2[id1]==1)
    if(length(j1))
      ans[id1[j1]] <- -(log(p)+M1$logf[id1[j1]]+
                          M2$logf[j1]+CP$logc[j1])
  }

  if(any(!is.finite(ans))) rep(1e10,nrow(dat)) else ans
}

fit.constant <- function(dat){
  start <- c(qlogis(.65),0,1.5,.8,0,0,1.5,.8,0,log(1))
  fit <- optim(start,function(z) sum(nll.i(z,dat)),method="BFGS",
               hessian=TRUE,control=list(maxit=1500,reltol=1e-9))
  if(fit$convergence!=0) stop("non-convergence")

  z <- fit$par
  H <- (fit$hessian+t(fit$hessian))/2
  J <- numDeriv::jacobian(function(x) nll.i(x,dat),z)
  Hi <- solve(H)
  Vz <- Hi%*%crossprod(J)%*%t(Hi)

  est <- natural.par(z)
  D <- diag(10)
  D[1,1] <- est["p"]*(1-est["p"])
  D[5,5] <- est["sigma1"]
  D[9,9] <- est["sigma2"]
  V <- D%*%Vz%*%t(D)
  se <- sqrt(pmax(diag(V),0)); names(se) <- par.names

  phi <- exp(z[10])
  tau <- phi/(phi+2)

  c(est,setNames(se,paste0("se.",par.names)),
    tau.constant=tau,loglik=-fit$value)
}

result.names <- c(par.names,paste0("se.",par.names),
                  "tau.constant","loglik")

## ============================================================
## PARALLEL MONTE CARLO
## ============================================================
one.rep <- function(r,sc){
  ## Same seed as matching row in the original 12-scenario run.
  set.seed(2025000L+10000L*sc$id+r)
  gamma <- unlist(sc[c("gamma0","gamma1","gamma2")])
  dat <- simulate.data(gamma)
  rates <- attr(dat,"rates")

  ans <- tryCatch(fit.constant(dat),error=function(e) NULL)
  if(is.null(ans))
    return(c(setNames(rep(NA_real_,length(result.names)),result.names),
             converged=0,rates))
  c(ans,converged=1,rates)
}

worker.objects <- c(
  "n","uc","prob.trt","p.true","a.true","b.true",
  "sigma1.true","sigma2.true","margin","log1pexp","logexpm1",
  "rclayton","clayton","simulate.data","par.names","result.names",
  "natural.par","nll.i","fit.constant","one.rep"
)

cl <- makeCluster(cores)
clusterEvalQ(cl,{library(QRM);library(numDeriv)})
clusterExport(cl,worker.objects,envir=.GlobalEnv)

raw <- vector("list",3)
for(s in 1:3){
  cat("Constant-copula fit:",scenarios$dependence[s],"\n")
  sc <- as.list(scenarios[s,])
  raw[[s]] <- do.call(rbind,
    parLapplyLB(cl,seq_len(R),one.rep,sc=sc))
}
stopCluster(cl)

## ============================================================
## SUMMARIES
## ============================================================
true.base <- setNames(c(p.true,a.true,sigma1.true,
                        b.true,sigma2.true),
                      par.names[1:9])

parameter.summary <- do.call(rbind,lapply(1:3,function(s){
  M <- raw[[s]]
  keep <- M[,"converged"]==1
  E <- M[keep,par.names[1:9],drop=FALSE]
  S <- M[keep,paste0("se.",par.names[1:9]),drop=FALSE]

  out <- do.call(rbind,lapply(1:9,function(j){
    eh <- E[,j]; sh <- S[,j]; th <- true.base[j]
    data.frame(
      dependence=scenarios$dependence[s],
      parameter=names(true.base)[j],
      true=th,mean=mean(eh),bias=mean(eh)-th,
      empirical.sd=sd(eh),
      rmse=sqrt(mean((eh-th)^2)),
      mean.se=mean(sh,na.rm=TRUE),
      coverage=mean(eh-1.96*sh<=th & eh+1.96*sh>=th,na.rm=TRUE)
    )
  }))
  out
}))

scenario.summary <- do.call(rbind,lapply(1:3,function(s){
  M <- raw[[s]]
  keep <- M[,"converged"]==1
  tau.true <- tau.scenarios[[scenarios$dependence[s]]]
  tau.hat <- M[keep,"tau.constant"]

  data.frame(
    dependence=scenarios$dependence[s],
    n=n,target.censor=target.censor,uc=uc,
    converged=sum(keep),convergence.rate=mean(keep),
    censor1.sus=mean(M[,"censor1.sus"],na.rm=TRUE),
    noevent1.all=mean(M[,"noevent1.all"],na.rm=TRUE),
    censor2=mean(M[,"censor2"],na.rm=TRUE),
    tau0.true=tau.true[1],tau1.true=tau.true[2],tau2.true=tau.true[3],
    mean.constant.tau=mean(tau.hat),
    sd.constant.tau=sd(tau.hat),
    mean.abs.tau.discrepancy=
      mean(rowMeans(abs(outer(tau.hat,tau.true,"-"))))
  )
}))

## gamma0 has a genuine truth only under the null setting.
M0 <- raw[[1]]
k0 <- M0[,"converged"]==1
g0.true <- tau2gamma(tau.scenarios$null)["gamma0"]
gamma0.null <- data.frame(
  true=g0.true,
  mean=mean(M0[k0,"gamma0"]),
  bias=mean(M0[k0,"gamma0"])-g0.true,
  empirical.sd=sd(M0[k0,"gamma0"]),
  mean.se=mean(M0[k0,"se.gamma0"]),
  coverage=mean(M0[k0,"gamma0"]-1.96*M0[k0,"se.gamma0"]<=g0.true &
                M0[k0,"gamma0"]+1.96*M0[k0,"se.gamma0"]>=g0.true)
)

write.csv(parameter.summary,"constant_parameter_summary.csv",row.names=FALSE)
write.csv(scenario.summary,"constant_scenario_summary.csv",row.names=FALSE)
write.csv(gamma0.null,"constant_gamma0_null_summary.csv",row.names=FALSE)
saveRDS(list(design=scenarios,raw=raw),
        "constant_copula_raw_results.rds")

cat("\nScenario summary:\n")
print(scenario.summary)
cat("\nNull gamma0 summary:\n")
print(gamma0.null)
