library(QRM);library(numDeriv);library(parallel)

## ============================================================
## 1. DESIGN
## ============================================================
R <- 1000L
dc <- parallel::detectCores()
cores <- if(is.na(dc)) 2L else max(1L, dc - 1L)
prob.trt <- c(.30, .30, .40)              # Trt 1, Trt 2, reference
p.true <- .70
a.true <- c(0, 2, 1); b.true <- c(0, 2, 1)
sigma1.true <- sigma2.true <- 1

tau.scenarios <- list(
  null     = c(.40, .40, .40),
  moderate = c(.40, .45, .60),            # current manuscript setting
  strong   = c(.40, .55, .70)
)

tau2gamma <- function(tau){
  phi <- 2 * tau / (1 - tau)
  c(gamma0 = log(phi[1]),
    gamma1 = log(phi[2] / phi[1]),
    gamma2 = log(phi[3] / phi[1]))
}

G <- t(vapply(tau.scenarios, tau2gamma, numeric(3)))
scenarios <- expand.grid(
  n = c(250L, 500L),
  censor = c(.20, .40),
  dependence = names(tau.scenarios),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
scenarios$id <- seq_len(nrow(scenarios))
scenarios <- scenarios[, c("id","n","censor","dependence")]
scenarios <- cbind(scenarios, G[scenarios$dependence, , drop = FALSE])

## ============================================================
## 2. DISTRIBUTION / COPULA FUNCTIONS
## ============================================================
margin <- function(t, mu, sigma){
  z <- (log(t) - mu) / sigma
  F <- pmin(pmax(QRM::pGumbel(z, 0, 1), 1e-12), 1 - 1e-12)
  logf <- QRM::dGumbel(z, 0, 1, log = TRUE) - log(sigma) - log(t)
  list(F = F, logf = logf)
}

log1pexp <- function(x){
  i <- x > 35
  out <- numeric(length(x))
  out[i] <- x[i] + log1p(exp(-x[i]))
  out[!i] <- log1p(exp(x[!i]))
  out
}

logexpm1 <- function(x){
  x <- pmax(x, 1e-14)
  i <- x > 35
  out <- numeric(length(x))
  out[i] <- x[i] + log1p(-exp(-x[i]))
  out[!i] <- log(expm1(x[!i]))
  out
}

rclayton <- function(phi){
  u <- pmin(pmax(runif(length(phi)), 1e-12), 1 - 1e-12)
  w <- pmin(pmax(runif(length(phi)), 1e-12), 1 - 1e-12)
  logv <- -log1pexp(-phi * log(u) +
                      logexpm1(-phi / (1 + phi) * log(w))) / phi
  cbind(u, pmin(pmax(exp(logv), 1e-12), 1 - 1e-12))
}

clayton <- function(u, v, phi){
  u <- pmin(pmax(u, 1e-12), 1 - 1e-12)
  v <- pmin(pmax(v, 1e-12), 1 - 1e-12)
  lu <- log(u); lv <- log(v)
  aa <- -phi * lu; bb <- -phi * lv; mm <- pmax(aa, bb)
  logA <- mm + log(exp(aa-mm) + exp(bb-mm) - exp(-mm))
  logC1 <- (-phi-1) * lu + (-1/phi-1) * logA
  C1 <- pmin(pmax(exp(pmin(logC1, 0)), 1e-12), 1 - 1e-12)
  logc <- log1p(phi) + (-1-phi)*(lu+lv) + (-2-1/phi)*logA
  list(C1 = C1, logc = logc)
}

## ============================================================
## 3. EMPIRICAL CALIBRATION OF U_C
##    Target = P(C < T1 | susceptible)
## ============================================================
set.seed(2026)
B.cal <- 300000L
trt.cal <- sample(1:3, B.cal, TRUE, prob.trt)
x1.cal <- as.integer(trt.cal == 1)
x2.cal <- as.integer(trt.cal == 2)
mu1.cal <- a.true[1] + a.true[2]*x1.cal + a.true[3]*x2.cal
T1.cal <- exp(mu1.cal +
                sigma1.true * QRM::qGumbel(runif(B.cal), 0, 1))

calibrate.uc <- function(target){
  f <- function(uc) mean(pmin(T1.cal / uc, 1)) - target
  hi <- 10
  while(f(hi) > 0) hi <- hi * 2
  uniroot(f, c(1e-8, hi), tol = 1e-8)$root
}

uc.map <- c(`0.2` = calibrate.uc(.20),
            `0.4` = calibrate.uc(.40))
scenarios$uc <- unname(uc.map[as.character(scenarios$censor)])

cat("\nCalibrated u_c values:\n")
print(data.frame(target.censor = c(.20,.40), uc = unname(uc.map)))

## ============================================================
## 4. DATA GENERATOR
## ============================================================
simulate.data <- function(n, gamma, uc){
  trt <- sample(1:3, n, TRUE, prob.trt)
  x1 <- as.integer(trt == 1)
  x2 <- as.integer(trt == 2)

  mu1 <- a.true[1] + a.true[2]*x1 + a.true[3]*x2
  mu2 <- b.true[1] + b.true[2]*x1 + b.true[3]*x2
  phi <- exp(gamma[1] + gamma[2]*x1 + gamma[3]*x2)

  UV <- rclayton(phi)
  T1 <- exp(mu1 + sigma1.true * QRM::qGumbel(UV[,1], 0, 1))
  T2 <- exp(mu2 + sigma2.true * QRM::qGumbel(UV[,2], 0, 1))

  susceptible <- rbinom(n, 1, p.true)
  T1.obs <- ifelse(susceptible == 1, T1, Inf)
  C <- runif(n, 0, uc)

  d1 <- as.integer(T1.obs <= C)
  y1 <- pmin(T1.obs, C)
  d2 <- numeric(n)
  y2 <- rep(NA_real_, n)

  id1 <- which(d1 == 1)
  if(length(id1)){
    rem <- C[id1] - T1[id1]
    d2[id1] <- as.integer(T2[id1] <= rem)
    y2[id1] <- pmin(T2[id1], rem)
  }

  dat <- data.frame(x1=x1, x2=x2, y1=y1, y2=y2, d1=d1, d2=d2)
  attr(dat, "rates") <- c(
    censor1.sus = mean(C[susceptible==1] < T1[susceptible==1]),
    noevent1.all = mean(d1 == 0),
    censor2 = if(length(id1)) mean(d2[id1] == 0) else NA_real_
  )
  dat
}

## ============================================================
## 5. ONE-STAGE MLE + SANDWICH
## ============================================================
par.names <- c("p","a0","a1","a2","sigma1",
               "b0","b1","b2","sigma2",
               "gamma0","gamma1","gamma2")
result.names <- c(par.names, paste0("se.",par.names),
                  "wald.p","tau0","tau1","tau2")

natural.par <- function(z)
  setNames(c(plogis(z[1]), z[2:4], exp(z[5]),
             z[6:8], exp(z[9]), z[10:12]), par.names)

nll.i <- function(z, dat){
  p <- plogis(z[1])
  a <- z[2:4]; s1 <- exp(z[5])
  b <- z[6:8]; s2 <- exp(z[9])
  g <- z[10:12]

  mu1 <- a[1] + a[2]*dat$x1 + a[3]*dat$x2
  M1 <- margin(dat$y1, mu1, s1)
  ans <- numeric(nrow(dat))

  id0 <- which(dat$d1 == 0)
  ans[id0] <- -log1p(-p * M1$F[id0])

  id1 <- which(dat$d1 == 1)
  if(length(id1)){
    mu2 <- b[1] + b[2]*dat$x1[id1] + b[3]*dat$x2[id1]
    M2 <- margin(dat$y2[id1], mu2, s2)
    phi <- exp(g[1] + g[2]*dat$x1[id1] + g[3]*dat$x2[id1])
    CP <- clayton(M1$F[id1], M2$F, phi)

    j0 <- which(dat$d2[id1] == 0)
    if(length(j0))
      ans[id1[j0]] <- -(log(p) + M1$logf[id1[j0]] +
                          log1p(-CP$C1[j0]))

    j1 <- which(dat$d2[id1] == 1)
    if(length(j1))
      ans[id1[j1]] <- -(log(p) + M1$logf[id1[j1]] +
                          M2$logf[j1] + CP$logc[j1])
  }

  if(any(!is.finite(ans))) rep(1e10, nrow(dat)) else ans
}

fit.model <- function(dat){
  start <- c(qlogis(.65), 0, 1.5, .8, 0,
             0, 1.5, .8, 0, log(1), 0, 0)

  fit <- optim(start, function(z) sum(nll.i(z, dat)), method="BFGS",
               hessian=TRUE, control=list(maxit=1500, reltol=1e-9))
  if(fit$convergence != 0) stop("non-convergence")

  z <- fit$par
  H <- (fit$hessian + t(fit$hessian))/2
  J <- numDeriv::jacobian(function(x) nll.i(x, dat), z)
  Hi <- solve(H)
  Vz <- Hi %*% crossprod(J) %*% t(Hi)

  est <- natural.par(z)
  D <- diag(12)
  D[1,1] <- est["p"]*(1-est["p"])
  D[5,5] <- est["sigma1"]
  D[9,9] <- est["sigma2"]
  V <- D %*% Vz %*% t(D)
  se <- sqrt(pmax(diag(V), 0))
  names(se) <- par.names

  W <- drop(t(z[11:12]) %*%
              solve(Vz[11:12,11:12,drop=FALSE]) %*% z[11:12])
  pW <- pchisq(W, 2, lower.tail=FALSE)

  c(est, setNames(se, paste0("se.",par.names)),
    wald.p=pW,
    tau0=exp(z[10])/(exp(z[10])+2),
    tau1=exp(z[10]+z[11])/(exp(z[10]+z[11])+2),
    tau2=exp(z[10]+z[12])/(exp(z[10]+z[12])+2))
}

## ============================================================
## 6. PARALLEL MONTE CARLO: ALL 12 SCENARIOS
## ============================================================
one.rep <- function(r, sc){
  set.seed(2025000L + 10000L*sc$id + r)
  gamma <- unlist(sc[c("gamma0","gamma1","gamma2")])
  dat <- simulate.data(sc$n, gamma, sc$uc)
  rates <- attr(dat, "rates")

  ans <- tryCatch(fit.model(dat), error=function(e) NULL)
  if(is.null(ans))
    return(c(setNames(rep(NA_real_, length(result.names)), result.names),
             converged=0, rates))

  c(ans, converged=1, rates)
}

worker.objects <- c(
  "prob.trt","p.true","a.true","b.true","sigma1.true","sigma2.true",
  "margin","log1pexp","logexpm1","rclayton","clayton",
  "simulate.data","par.names","result.names","natural.par",
  "nll.i","fit.model","one.rep"
)

cl <- makeCluster(cores)
clusterEvalQ(cl, {library(QRM); library(numDeriv)})
clusterExport(cl, worker.objects, envir=.GlobalEnv)

raw <- vector("list", nrow(scenarios))
for(s in seq_len(nrow(scenarios))){
  cat(sprintf("Scenario %02d/%02d: n=%d, censor=%.0f%%, %s\n",
              s, nrow(scenarios), scenarios$n[s],
              100*scenarios$censor[s], scenarios$dependence[s]))
  sc <- as.list(scenarios[s,])
  raw[[s]] <- do.call(rbind,
    parLapplyLB(cl, seq_len(R), one.rep, sc=sc))
}
stopCluster(cl)

## ============================================================
## 7. SUMMARIES
## ============================================================
true.base <- c(p.true, a.true, sigma1.true, b.true, sigma2.true)

summarize.one <- function(M, sc){
  keep <- M[,"converged"] == 1 & apply(M[,par.names,drop=FALSE],1,
                                       function(x) all(is.finite(x)))
  E <- M[keep, par.names, drop=FALSE]
  S <- M[keep, paste0("se.",par.names), drop=FALSE]
  gamma <- unlist(sc[c("gamma0","gamma1","gamma2")])
  truth <- setNames(c(true.base, gamma), par.names)

  parameter <- do.call(rbind, lapply(seq_along(par.names), function(j){
    eh <- E[,j]; sh <- S[,j]; th <- truth[j]
    data.frame(
      parameter=par.names[j], true=th,
      mean=mean(eh), bias=mean(eh)-th,
      empirical.sd=sd(eh),
      rmse=sqrt(mean((eh-th)^2)),
      mean.se=mean(sh,na.rm=TRUE),
      coverage=mean(eh-1.96*sh <= th & eh+1.96*sh >= th,na.rm=TRUE)
    )
  }))

  tau.true <- tau.scenarios[[sc$dependence]]
  tau <- do.call(rbind, lapply(0:2, function(j){
    x <- M[keep,paste0("tau",j)]
    data.frame(group=j, true=tau.true[j+1], mean=mean(x),
               bias=mean(x)-tau.true[j+1], empirical.sd=sd(x),
               rmse=sqrt(mean((x-tau.true[j+1])^2)))
  }))

  scenario <- data.frame(
    n=sc$n, target.censor=sc$censor, dependence=sc$dependence, uc=sc$uc,
    converged=sum(keep), convergence.rate=mean(keep),
    censor1.sus=mean(M[,"censor1.sus"],na.rm=TRUE),
    noevent1.all=mean(M[,"noevent1.all"],na.rm=TRUE),
    censor2=mean(M[,"censor2"],na.rm=TRUE),
    wald.rejection=mean(M[keep,"wald.p"] < .05,na.rm=TRUE)
  )

  list(parameter=parameter, tau=tau, scenario=scenario)
}

S <- lapply(seq_len(nrow(scenarios)),
            function(i) summarize.one(raw[[i]], as.list(scenarios[i,])))

parameter.summary <- do.call(rbind, lapply(seq_along(S), function(i)
  cbind(scenarios[i,c("id","n","censor","dependence")], S[[i]]$parameter)))

tau.summary <- do.call(rbind, lapply(seq_along(S), function(i)
  cbind(scenarios[i,c("id","n","censor","dependence")], S[[i]]$tau)))

scenario.summary <- do.call(rbind, lapply(S, `[[`, "scenario"))

write.csv(parameter.summary, "mc_parameter_summary.csv", row.names=FALSE)
write.csv(tau.summary, "mc_tau_summary.csv", row.names=FALSE)
write.csv(scenario.summary, "mc_scenario_summary.csv", row.names=FALSE)
saveRDS(list(design=scenarios, raw=raw), "mc_raw_results.rds")

cat("\nScenario summary:\n")
print(scenario.summary)
cat("\nFiles written: mc_parameter_summary.csv, mc_tau_summary.csv,",
    "mc_scenario_summary.csv, mc_raw_results.rds\n")
