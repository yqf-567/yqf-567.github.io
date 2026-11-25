# Three focused experiments for recurrent events with a terminal death event.
# Based on recurrent_design_simulations.R; requires only base R and parallel.
# Select all and run in R/RStudio, or use Rscript. See the companion README.

# =========================== MANUAL SETTINGS ===========================
RUN_MODE <- "full"       # CHOOSE "quick" (100) OR "full" (1000 replicates).
N_WORKERS <- 4L           # 1 = serial; >=2 = Windows-compatible parallel workers.
OUTPUT_DIR <- file.path(getwd(), "simulation_results_death", RUN_MODE)
SEED <- 20261117L
ALLOW_SERIAL_FALLBACK <- TRUE
# Optional: --mode=full --workers=4 --out=C:/data/death_results
# Developer execution check: --mode=smoke (10 replicates, ALL selected cases).
# =======================================================================

# Key changes:
# * Death is an OUTCOME, not ordinary censoring. W(t)=0 after death.
# * Only six scientifically chosen scenarios; two weighted-score sensitivities.
# * Exp 1: n=100,500, short/full observation; bias, 95% coverage, variance ratio.
# * Exp 2: B=300, short/H-optimal/full; duplicate durations are merged.
# * Exp 3: three dependent-death scenarios, m=50,200, B=300.
# * Costs: planned E(C)=L, or actual E(min(C,D)); identical cases are merged.
# * Full = 1000 replicates and 31 grid points (previously 5000 and 81).
# * Complete pilots observe each patient through death or tau, with W=0 later.
# * Nonparametric analysis: generation parameters are NOT given to the pilot.
# * Exactly three CSVs; reference diagnostics stay in memory and the console.
# * Simulation MCSEs exclude reference integration error, diagnosed separately.
# * Oracle comparisons are restricted to the specified finite duration grid.

options(stringsAsFactors=FALSE)
default_output_dir <- file.path(getwd(),"simulation_results_death",RUN_MODE)
for (arg in commandArgs(trailingOnly=TRUE)) {
  if (grepl("^--mode=",arg)) RUN_MODE <- sub("^--mode=","",arg)
  if (grepl("^--workers=",arg)) N_WORKERS <- as.integer(sub("^--workers=","",arg))
}
# Resolve the default folder AFTER a command-line mode override.
if (!any(grepl("^--out=",commandArgs(trailingOnly=TRUE)))) {
  # Keep a manually edited OUTPUT_DIR; update only the supplied default.
  if (identical(OUTPUT_DIR,default_output_dir))
    OUTPUT_DIR <- file.path(getwd(),"simulation_results_death",RUN_MODE)
}
for (arg in commandArgs(trailingOnly=TRUE))
  if (grepl("^--out=",arg)) OUTPUT_DIR <- sub("^--out=","",arg)
stopifnot(RUN_MODE %in% c("smoke","quick","full"),length(N_WORKERS)==1L,
          is.finite(N_WORKERS),N_WORKERS>=1L,N_WORKERS==as.integer(N_WORKERS))
settings <- switch(RUN_MODE,
  smoke=list(R=10L,grid=11L,ref_n=512L,ref_shifts=4L,cal_n=4096L,batch=5L),
  quick=list(R=100L,grid=21L,ref_n=4096L,ref_shifts=8L,cal_n=32768L,batch=25L),
  full=list(R=1000L,grid=31L,ref_n=16384L,ref_shifts=8L,cal_n=65536L,batch=50L))
cfg <- c(settings,list(mode=RUN_MODE,seed=SEED,tau=1,min_G=.1,
                      n=c(100L,500L),B=300,m=c(50L,200L)))
dir.create(OUTPUT_DIR,recursive=TRUE,showWarnings=FALSE)
if (!dir.exists(OUTPUT_DIR)) stop("Cannot create OUTPUT_DIR.")
OUTPUT_DIR <- normalizePath(OUTPUT_DIR,winslash="/",mustWork=TRUE)
cat(sprintf("\nMODE: %s | %d replicates/setting | requested workers: %d\nOUTPUT: %s\n",
            RUN_MODE,cfg$R,N_WORKERS,OUTPUT_DIR))
if (RUN_MODE!="full") cat("This is an execution/exploratory run, not the final study.\n")
cat(sprintf("Coverage MCSE near 95%% is approximately %.2f percentage points.\n",
            100*sqrt(.95*.05/cfg$R)))
RNGkind("L'Ecuyer-CMRG")
set.seed(cfg$seed)

x_grid <- seq(-log(cfg$min_G),0,length.out=cfg$grid)
duration <- function(x,tau) ifelse(x==0,tau,tau*(-expm1(-x))/ifelse(x==0,1,x))
design <- data.frame(j=seq_along(x_grid),x=x_grid,L=duration(x_grid,cfg$tau))
scores <- list(R2=c(0,1,0),W3=rep(1/3,3))
scenarios <- data.frame(
  id=c("no_death_early","independent_death20_early","dependent_death20_early",
       "dependent_death40_early","weibull_death20_early","dependent_death20_moderate"),
  mechanism=c(rep("exponential",4),"dependent_weibull","exponential"),
  mean_latent_T2=c(rep(.24,5),1),
  death_probability=c(0,.2,.2,.4,.2,.2),
  death_log_multiplier=c(0,0,rep(log(2),4)))
# Main score in all six; weighted sensitivity only for two distinct mechanisms.
cases <- data.frame(s=seq_len(nrow(scenarios)),score="R2")
cases <- rbind(cases,data.frame(s=c(3L,5L),score="W3"))
cases$key <- paste(scenarios$id[cases$s],cases$score,sep="__")
pilot_scenarios <- c("dependent_death20_early","weibull_death20_early",
                     "dependent_death20_moderate")

halton_points <- function(n) {
  radical <- function(base) {
    k <- seq_len(n); f <- 1/base; value <- numeric(n)
    while(any(k>0)) {value <- value+f*(k %% base); k <- k %/% base; f <- f/base}
    value
  }
  vapply(c(2L,3L,5L,7L,11L),radical,numeric(n))
}
shift_points <- function(H) pmin(pmax(sweep(H,2L,runif(ncol(H)),"+") %% 1,1e-12),1-1e-12)
latent_events <- function(U,sc) {
  if(sc$mechanism=="dependent_weibull") {
    gaps <- (-log(U[,1:3,drop=FALSE]))^(1/.7)/gamma(1+1/.7)
    gaps <- gaps*exp(.8*qnorm(U[,4L])-.8^2/2)
  } else gaps <- -log(U[,1:3,drop=FALSE])
  gaps <- gaps*sc$scale
  cbind(gaps[,1],gaps[,1]+gaps[,2],rowSums(gaps))
}
# Death hazard while alive is lambda_D*exp(alpha*min(N(t),3)). Only the first
# three potential recurrences are needed for both scores and this death hazard.
# Potential recurrences after death are NEVER counted as observed events.
hazard_exposure <- function(T,alpha,tau) {
  ends <- pmin(T,tau)
  widths <- cbind(ends[,1],ends[,2]-ends[,1],ends[,3]-ends[,2],tau-ends[,3])
  as.vector(widths %*% exp(alpha*(0:3)))
}
death_times <- function(T,U,lambda,alpha) {
  n <- nrow(T)
  if(lambda==0) return(rep(Inf,n))
  remaining <- -log(pmax(U,.Machine$double.xmin))
  D <- rep(Inf,n); active <- rep(TRUE,n); start <- numeric(n)
  for(k in 0:3) {
    end <- if(k<3) T[,k+1L] else rep(Inf,n)
    rate <- lambda*exp(alpha*k)
    capacity <- rate*(end-start)
    die <- active & remaining<=capacity
    D[die] <- start[die]+remaining[die]/rate
    active[die] <- FALSE
    remaining[active] <- remaining[active]-capacity[active]
    start <- end
  }
  D
}
histories_from_uniforms <- function(U,sc) {
  T <- latent_events(U,sc)
  D <- death_times(T,U[,5L],sc$lambda_D,sc$death_log_multiplier)
  # U_k=min(T_k,D): W(t)=sum_k w_k I(t<U_k). Equal endpoints at death are valid.
  list(ends=pmin(T,matrix(D,nrow(T),3L)),D=D)
}
histories <- function(n,sc) histories_from_uniforms(matrix(runif(n*5L),n,5L),sc)
area <- function(T,w,tau) as.vector(pmin(T,tau) %*% w)
primitive <- function(t,x,tau) if(x==0) t else (tau/x)*expm1(t*x/tau)
cost_integral <- function(D,x,tau) {
  d <- pmin(D,tau)
  if(x==0) d else (tau/x)*(-expm1(-d*x/tau))
}
observed_z <- function(T,w,x,tau,U) {
  n <- nrow(T)
  C <- if(x==0) rep(tau,n) else pmin(-log(pmax(U,.Machine$double.xmin))*tau/x,tau)
  as.vector(primitive(pmin(T,matrix(C,n,3L),tau),x,tau) %*% w)
}
# E_C[Z^2 | joint outcome history], integrated exactly between score changes.
# Using min(T_k,D) makes the original variance integration valid after death.
conditional_q <- function(T,w,x,tau) {
  n <- nrow(T); ends <- pmin(T,tau)
  starts <- cbind(rep(0,n),ends[,1],ends[,2]); widths <- ends-starts
  levels <- rev(cumsum(rev(w))); tail_area <- Q <- numeric(n)
  for(k in 3:1) {
    u <- starts[,k]; d <- widths[,k]; h <- levels[k]
    if(x==0) {i0 <- d; i1 <- d*d/2} else {
      a <- tau/x; z <- d/a; delta <- expm1(z)-z
      small <- abs(z)<1e-4
      delta[small] <- z[small]^2*(.5+z[small]*(1/6+z[small]*(1/24+z[small]/120)))
      i0 <- a*exp(u/a)*expm1(z); i1 <- a*a*exp(u/a)*delta
    }
    Q <- Q+2*h*(tail_area*i0+h*i1)
    tail_area <- tail_area+h*d
  }
  Q
}
quad <- function(f,tau) integrate(f,0,tau,rel.tol=1e-10,abs.tol=1e-12,subdivisions=300L)$value

self_check <- function() {
  latent <- rbind(c(.2,.5,.9),c(.7,1.2,1.5),c(1.1,1.5,2),c(.1,.1,.1))
  D <- c(.1,.9,Inf,.1)
  T <- pmin(latent,matrix(D,4,3))
  for(w in scores) {
    A <- area(T,w,1)
    stopifnot(max(abs(conditional_q(T,w,0,1)-A^2))<1e-10,
              max(abs(observed_z(T,w,0,1,c(.2,.5,.8,.4))-A))<1e-10)
    for(x in c(.3,2)) for(i in seq_len(nrow(T))) {
      cuts <- sort(unique(c(0,pmin(T[i,],1),1))); ans <- 0
      for(k in seq_len(length(cuts)-1L)) {
        f <- function(s) vapply(s,function(ss) 2*sum(w*(ss<T[i,]))*
          exp(ss*x)*sum(w*pmax(pmin(T[i,],1)-ss,0)),numeric(1))
        ans <- ans+integrate(f,cuts[k],cuts[k+1L],rel.tol=1e-10)$value
      }
      stopifnot(abs(ans-conditional_q(T,w,x,1)[i])<1e-8)
    }
  }
  uniform <- c(.2,.4,.6,.8)
  stopifnot(max(abs(death_times(latent,uniform,.3,0)+log(uniform)/.3))<1e-12)
  # Check the inverse death hazard both before and after recurrent events.
  dd <- death_times(latent,uniform,.3,log(2))
  for(i in seq_along(dd))
    stopifnot(abs(.3*hazard_exposure(latent[i,,drop=FALSE],log(2),dd[i])+log(uniform[i]))<1e-10)
  stopifnot(max(abs(vapply(design$x,function(x) cost_integral(Inf,x,1),numeric(1))-design$L))<1e-12,
            all(diff(design$L)>0))
  for(x in c(0,.3,2)) for(d in c(.1,.9,2))
    stopifnot(abs(cost_integral(d,x,1)-quad(function(t) exp(-x*t),min(d,1)))<1e-10)
}
self_check()
cat("Death, score integration, and cost checks passed. Preparing reference quantities...\n")

refs <- list(); reference_diagnostics <- list()
Hbase <- halton_points(cfg$ref_n)
for(s in seq_len(nrow(scenarios))) {
  sc <- as.list(scenarios[s,]); sc$scale <- sc$mean_latent_T2/2
  if(sc$death_probability==0) sc$lambda_D <- 0 else if(sc$death_log_multiplier==0) {
    sc$lambda_D <- -log1p(-sc$death_probability)/cfg$tau
  } else {
    set.seed(cfg$seed+100L+s)
    calT <- latent_events(shift_points(halton_points(cfg$cal_n)),sc)
    exposure <- hazard_exposure(calT,sc$death_log_multiplier,cfg$tau)
    sc$lambda_D <- uniroot(function(lambda) mean(-expm1(-lambda*exposure))-sc$death_probability,
      c(0,-log1p(-sc$death_probability)/cfg$tau),tol=1e-12)$root
  }
  score_names <- cases$score[cases$s==s]
  mu_blocks <- matrix(0,cfg$ref_shifts,length(score_names))
  q_blocks <- array(0,c(cfg$ref_shifts,cfg$grid,length(score_names)))
  cost_blocks <- matrix(0,cfg$ref_shifts,cfg$grid)
  death_blocks <- numeric(cfg$ref_shifts)
  for(b in seq_len(cfg$ref_shifts)) {
    set.seed(cfg$seed+1000L+100L*s+b)
    sample <- histories_from_uniforms(shift_points(Hbase),sc)
    death_blocks[b] <- mean(sample$D<=cfg$tau)
    for(j in seq_len(cfg$grid))
      cost_blocks[b,j] <- mean(cost_integral(sample$D,design$x[j],cfg$tau))
    for(k in seq_along(score_names)) {
      w <- scores[[score_names[k]]]
      mu_blocks[b,k] <- mean(area(sample$ends,w,cfg$tau))
      for(j in seq_len(cfg$grid))
        q_blocks[b,j,k] <- mean(conditional_q(sample$ends,w,design$x[j],cfg$tau))
    }
  }
  actual_cost <- colMeans(cost_blocks)
  cost_se <- apply(cost_blocks,2,sd)/sqrt(cfg$ref_shifts)
  for(k in seq_along(score_names)) {
    score <- score_names[k]; w <- scores[[score]]
    mu <- mean(mu_blocks[,k]); mu_se <- sd(mu_blocks[,k])/sqrt(cfg$ref_shifts)
    v <- colMeans(q_blocks[,,k])-mu^2
    vse <- apply(q_blocks[,,k]-2*mu*mu_blocks[,k],2,sd)/sqrt(cfg$ref_shifts)
    # Independent exponential recurrence/death gives an analytic survival anchor.
    exact <- sc$mechanism=="exponential" && sc$death_log_multiplier==0 && score=="R2"
    if(exact) {
      rate <- 1/sc$scale
      surv <- function(t) exp(-(rate+sc$lambda_D)*t)*(1+rate*t)
      mu_exact <- quad(surv,cfg$tau)
      v_exact <- vapply(design$x,function(x)
        2*quad(function(t) surv(t)*primitive(t,x,cfg$tau),cfg$tau)-mu_exact^2,numeric(1))
      if(abs(mu-mu_exact)>max(8*mu_se,.002))
        stop("Reference integration failed the independent-death analytic anchor.")
      mu <- mu_exact; mu_se <- 0; v <- v_exact; vse <- rep(0,cfg$grid)
      actual_cost <- vapply(design$x,function(x)
        cost_integral(cfg$tau,x+sc$lambda_D*cfg$tau,cfg$tau),numeric(1))
      cost_se <- rep(0,cfg$grid)
    }
    if(any(!is.finite(v)) || any(v<=0)) stop("Invalid reference variance.")
    if(any(actual_cost<=0) || any(actual_cost>design$L+1e-10)) stop("Invalid cost reference.")
    key <- paste(sc$id,score,sep="__")
    refs[[key]] <- list(sc=sc,score=score,w=w,mu=mu,mu_se=mu_se,v=v,vse=vse,
                       actual_cost=actual_cost,cost_se=cost_se)
    # Keep accuracy diagnostics separate from concise publication summaries.
    ratio <- mu_se/sqrt(min(v)/max(cfg$n)/cfg$R)
    reference_diagnostics[[key]] <- list(mu_mcse=mu_se,variance_mcse=vse,cost_mcse=cost_se,
      death_probability=mean(death_blocks),death_mcse=sd(death_blocks)/sqrt(cfg$ref_shifts),
      mean_reference_to_bias_mcse=ratio,lambda_D=sc$lambda_D,
      method=if(exact) "quadrature" else "randomized_Halton")
    cat(sprintf("  %-35s %-2s | death=%.3f | mu=%.5f | mean reference MCSE=%.2g\n",
                sc$id,score,mean(death_blocks),mu,mu_se))
    if(ratio>.25 && RUN_MODE!="smoke")
      warning("Reference mean uncertainty is material for ",key,"; increase ref_n for final bias assessment.",call.=FALSE)
    if(max(vse/v)>.02 && RUN_MODE!="smoke")
      warning("Reference variance uncertainty exceeds 2% for ",key,"; increase ref_n.",call.=FALSE)
  }
}
rm(Hbase,sample,mu_blocks,q_blocks,cost_blocks)

cost_for <- function(ref,basis) if(basis=="actual") ref$actual_cost else design$L
bases_for <- function(ref) if(ref$sc$death_probability==0) "both" else c("planned","actual")
sample_sizes <- function(B,cost) {
  # No upward rounding tolerance: floor preserves the expected-budget constraint.
  n <- floor(B/cost)
  if(any(n<2) || any(!is.finite(n))) stop("Budget gives invalid sample size.")
  as.integer(n)
}
metrics <- function(z,mu) {
  err <- mean(z)-mu; vhat <- var(z)/length(z)
  c(error=err,squared_error=err^2,covered=as.numeric(abs(err)<=qnorm(.975)*sqrt(vhat)))
}
pack <- function(M,meta) c(meta,list(R=nrow(M),sum=colSums(M),sumsq=colSums(M^2)))

jobs <- list()
set.seed(cfg$seed+50000L)
next_seed <- .Random.seed
add_jobs <- function(e,key,size,basis="none",m=NA_integer_) {
  for(first in seq.int(1L,cfg$R,by=cfg$batch)) {
    id <- length(jobs)+1L
    jobs[[id]] <<- list(id=id,e=e,key=key,size=size,basis=basis,m=m,
      count=min(cfg$batch,cfg$R-first+1L),seed=next_seed)
    next_seed <<- parallel::nextRNGStream(next_seed)
  }
}
for(key in names(refs)) for(n in cfg$n) add_jobs(1L,key,n)
for(key in names(refs)) for(basis in bases_for(refs[[key]])) add_jobs(2L,key,cfg$B,basis)
for(key in names(refs)) if(refs[[key]]$score=="R2" && refs[[key]]$sc$id %in% pilot_scenarios)
  for(basis in bases_for(refs[[key]])) for(m in cfg$m) add_jobs(3L,key,cfg$B,basis,m)

run_job <- function(job) {
  assign(".Random.seed",job$seed,envir=.GlobalEnv)
  ref <- refs[[job$key]]; sc <- ref$sc; w <- ref$w; tau <- cfg$tau; R <- job$count
  basic_names <- c("error","squared_error","covered")
  meta <- function(j,rule) list(group=paste(job$e,job$key,job$size,job$basis,job$m,j,sep="|"),
    e=job$e,key=job$key,size=job$size,basis=job$basis,m=job$m,j=j,rule=rule)
  if(job$e==1L) {
    js <- c(1L,cfg$grid)
    M <- lapply(js,function(j) matrix(0,R,3L,dimnames=list(NULL,basic_names)))
    for(r in seq_len(R)) {
      h <- histories(job$size,sc); U <- runif(job$size)
      for(k in seq_along(js))
        M[[k]][r,] <- metrics(observed_z(h$ends,w,design$x[js[k]],tau,U),ref$mu)
    }
    return(lapply(seq_along(js),function(k) pack(M[[k]],meta(js[k],"fixed"))))
  }
  B <- job$size; cost <- cost_for(ref,job$basis)
  ns <- sample_sizes(B,cost); H <- cost*ref$v; jH <- which.min(H)
  if(job$e==2L) {
    labels <- c(short=1L,H_optimal=jH,full=cfg$grid)
    js <- sort(unique(unname(labels)))
    rules <- vapply(js,function(j) paste(names(labels)[labels==j],collapse="+"),character(1))
    M <- lapply(js,function(j) matrix(0,R,3L,dimnames=list(NULL,basic_names)))
    for(r in seq_len(R)) {
      h <- histories(max(ns[js]),sc); U <- runif(nrow(h$ends))
      for(k in seq_along(js)) {
        idx <- seq_len(ns[js[k]])
        M[[k]][r,] <- metrics(observed_z(h$ends[idx,,drop=FALSE],w,design$x[js[k]],tau,U[idx]),ref$mu)
      }
    }
    return(lapply(seq_along(js),function(k) pack(M[[k]],meta(js[k],rules[k]))))
  }
  extra <- c("selected_L","H_loss_pct","budget_ratio")
  M <- matrix(0,R,length(basic_names)+length(extra),dimnames=list(NULL,c(basic_names,extra)))
  for(r in seq_len(R)) {
    pilot <- histories(job$m,sc); pmu <- mean(area(pilot$ends,w,tau))
    pv <- vapply(design$x,function(x) mean(conditional_q(pilot$ends,w,x,tau))-pmu^2,numeric(1))
    if(any(pv < -1e-10)) stop("Negative pilot variance beyond numerical tolerance.")
    pv <- pmax(pv,0)
    pc <- if(job$basis=="actual") vapply(design$x,function(x)
      mean(cost_integral(pilot$D,x,tau)),numeric(1)) else design$L
    j <- which.min(pc*pv)
    n <- sample_sizes(B,pc[j])  # Uses estimated cost, never the hidden true cost.
    h <- histories(n,sc)
    M[r,] <- c(metrics(observed_z(h$ends,w,design$x[j],tau,runif(n)),ref$mu),
                design$L[j],100*(H[j]/H[jH]-1),n*cost[j]/B)
  }
  list(pack(M,meta(NA_integer_,"pilot_min_H")))
}

run_all <- function() {
  cl <- NULL
  if(N_WORKERS>1L) {
    cat("Starting parallel workers...\n")
    cl <- tryCatch(parallel::makePSOCKcluster(N_WORKERS,master="127.0.0.1",setup_timeout=30),
      error=function(e) {
        if(!ALLOW_SERIAL_FALLBACK) stop(e)
        message("Parallel startup failed: ",conditionMessage(e),"\nContinuing SERIAL execution.")
        NULL
      })
    if(!is.null(cl)) {
      on.exit(parallel::stopCluster(cl),add=TRUE)
      parallel::clusterExport(cl,c("cfg","design","refs","latent_events","death_times",
        "histories_from_uniforms","histories","area","primitive","cost_integral","observed_z",
        "conditional_q","cost_for","sample_sizes","metrics","pack","run_job"),envir=environment(run_job))
    }
  }
  cat(if(is.null(cl)) "Execution: SERIAL\n" else paste0("Execution: PARALLEL (",N_WORKERS," workers)\n"))
  ans <- list(); labels <- c("VARIANCE VALIDATION","ALLOCATION CRITERION","PILOT SELECTION")
  wave <- max(4L,4L*N_WORKERS)
  for(e in 1:3) {
    these <- jobs[vapply(jobs,function(j) j$e==e,logical(1))]
    cat(sprintf("\nEXPERIMENT %d/3: %s | mode=%s | %d replicates/setting\n",e,labels[e],RUN_MODE,cfg$R))
    for(first in seq.int(1L,length(these),by=wave)) {
      idx <- first:min(length(these),first+wave-1L)
      batch <- if(is.null(cl)) lapply(these[idx],run_job) else parallel::parLapplyLB(cl,these[idx],run_job)
      ans <- c(ans,batch)
      cat(sprintf("  Experiment %d: %d/%d batches complete\n",e,max(idx),length(these)))
      flush.console()
    }
  }
  unlist(ans,recursive=FALSE)
}
batch_results <- run_all()

summarize_group <- function(parts) {
  z <- parts[[1L]]; ref <- refs[[z$key]]
  R <- sum(vapply(parts,function(p) p$R,numeric(1)))
  sums <- Reduce("+",lapply(parts,function(p) p$sum))
  squares <- Reduce("+",lapply(parts,function(p) p$sumsq))
  means <- sums/R
  variances <- pmax((squares-sums^2/R)/(R-1),0)
  ses <- sqrt(variances/R)
  stopifnot(identical(names(variances),names(sums)),all(is.finite(means)),all(is.finite(ses)))
  row <- list(mode=RUN_MODE,replicates=R,scenario=ref$sc$id,score=ref$score)
  if(z$e==1L) {
    row <- c(row,list(n=z$size,L=design$L[z$j],bias=means["error"],bias_mcse=ses["error"],
      coverage_95=means["covered"],variance_ratio=variances["error"]/(ref$v[z$j]/z$size)))
  } else if(z$e==2L) {
    cost <- cost_for(ref,z$basis); ns <- sample_sizes(z$size,cost); risk <- ref$v/ns
    row <- c(row,list(cost_basis=z$basis,rule=z$rule,L=design$L[z$j],n=ns[z$j],
      mse_ratio=means["squared_error"]/risk[z$j],
      efficiency_loss_pct=100*(risk[z$j]/min(risk)-1)))
  } else {
    cost <- cost_for(ref,z$basis); jH <- which.min(cost*ref$v)
    row <- c(row,list(cost_basis=z$basis,pilot_n=z$m,oracle_L=design$L[jH],
      mean_L=means["selected_L"],sd_L=sqrt(variances["selected_L"]),
      H_loss_pct=means["H_loss_pct"],H_loss_mcse=ses["H_loss_pct"],
      coverage_95=means["covered"],mean_budget_ratio=means["budget_ratio"]))
  }
  as.data.frame(lapply(row,unname),check.names=FALSE)
}

# Rerun write_summaries() in the same session to rebuild CSVs without simulations.
write_summaries <- function() {
  groups <- split(batch_results,vapply(batch_results,function(z) z$group,character(1)))
  filenames <- c("experiment_1_variance.csv","experiment_2_allocation.csv","experiment_3_pilot.csv")
  result <- vector("list",3L)
  for(e in 1:3) {
    selected <- groups[vapply(groups,function(g) g[[1L]]$e==e,logical(1))]
    tab <- do.call(rbind,lapply(selected,summarize_group)); rownames(tab) <- NULL
    numeric_columns <- vapply(tab,is.numeric,logical(1))
    stopifnot(all(tab$replicates==cfg$R),all(vapply(tab[numeric_columns],function(v) all(is.finite(v)),logical(1))))
    if("coverage_95" %in% names(tab)) stopifnot(all(tab$coverage_95>=0 & tab$coverage_95<=1))
    # Six significant digits, enough for small biases without long CSV decimals.
    tab[numeric_columns] <- lapply(tab[numeric_columns],signif,digits=6)
    dest <- file.path(OUTPUT_DIR,filenames[e]); write.csv(tab,dest,row.names=FALSE)
    cat(sprintf("WROTE experiment %d: %d rows x %d columns | %s\n",e,nrow(tab),ncol(tab),dest))
    result[[e]] <- tab
  }
  invisible(result)
}
summary_tables <- write_summaries()
cat("\nDONE:",RUN_MODE,"mode. Three concise CSVs written to",OUTPUT_DIR,"\n")
cat("Exp 2/3 use B=300; all times use tau=1. Full mode has 1000 replicates.\n")
cat("Exp 3 H_loss evaluates design quality at equal expected cost; mean_budget_ratio checks plug-in enrollment.\n")
cat("Reference diagnostics and batch_results remain in memory; no additional files are written.\n")
