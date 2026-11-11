library(mvtnorm)
library(numDeriv)
library(parallel)

m <- 3L
pos <- rbind(c(0,0),c(.31,0),c(.10,.35))
K <- exp(-as.matrix(dist(pos))/.45)
rtrue <- c(.55,.45,-.35)
labels <- c(paste0("intercept",1:3),paste0("beta",1:3),
            paste0("sigma",1:3),paste0("tau",1:3),
            "rho12","rho23","rho13|2")

Rvine <- function(z) {
  r <- tanh(z)
  a <- r[1]*r[2]+r[3]*sqrt((1-r[1]^2)*(1-r[2]^2))
  matrix(c(1,r[1],a,r[1],1,r[2],a,r[2],1),3,3)
}
# Gaussian likelihood with observed coordinates and lower-censored coordinates.
# Missing later gaps contribute neither densities nor survivor factors.
lcens <- function(mu,V,y,lo,exact,cens) {
  e <- which(exact); c <- which(cens)
  if (!length(e) && !length(c)) return(0)
  if (length(e)) {
    W <- V[e,e,drop=FALSE]
    d <- y[e]-mu[e]
    L <- chol(W)
    q <- backsolve(L, d, transpose=TRUE)
    ans <- -.5*(length(e)*log(2*pi)+2*sum(log(diag(L)))+sum(q^2))
  } else ans <- 0
  if (!length(c)) return(ans)
  if (length(e)) {
    M <- V[c,e,drop=FALSE] %*% solve(W)
    mc <- as.vector(mu[c]+M%*%(y[e]-mu[e]))
    vc <- V[c,c,drop=FALSE]-M%*%V[e,c,drop=FALSE]
  } else {
    mc <- mu[c]; vc <- V[c,c,drop=FALSE]
  }
  if (length(c)==1) {
    pr <- pnorm((lo[c]-mc)/sqrt(vc[1,1]),lower.tail=FALSE)
  } else {
    # Fixed seed ensures reproducible objective and numerical derivatives.
    set.seed(11323)
    pr <- as.numeric(pmvnorm(lower=lo[c],upper=rep(Inf,length(c)),
                            mean=mc,sigma=vc,
                            algorithm=GenzBretz(maxpts=1024,abseps=1e-5)))
  }
  ans+log(pmax(pr,1e-14))
}
covmat <- function(th,z) {
  P <- matrix(th,4,3)
  sg <- exp(P[3,]); tg <- exp(P[4,])
  # Order: site 1 gap 1:3, site 2 gap 1:3, ...
  kronecker(K,diag(tg^2,3))+kronecker(diag(m),diag(sg,3)%*%Rvine(z)%*%diag(sg,3))
}
block1 <- function(th,d) {
  P <- matrix(th,4,3)
  sum(vapply(1:3,function(j) {
    p <- P[,j]; mu <- p[1]+p[2]*d$x
    V <- exp(2*p[3])*diag(m)+exp(2*p[4])*K
    lcens(mu,V,d$y[,j],d$lo[,j],d$ex[,j],d$ce[,j])
  },0.0))
}
block2 <- function(th,z,d) {
  P <- matrix(th,4,3)
  mu <- as.vector(t(sweep(outer(d$x,P[2,]),2,P[1,],"+")))
  lcens(mu,covmat(th,z),as.vector(t(d$y)),as.vector(t(d$lo)),
        as.vector(t(d$ex)),as.vector(t(d$ce)))-block1(th,d)
}
# One subject-level Uniform(0,cu) censoring time; later gaps are missing.
generate <- function(n,tau,cu) {
  B <- n/m; x <- matrix(rnorm(n),B,m)
  R <- Rvine(atanh(rtrue))
  eps <- matrix(rmvnorm(n,sigma=R),n,3)*.85
  w <- array(0,c(B,m,3))
  for(j in 1:3) w[,,j] <- rmvnorm(B,sigma=tau^2*K)
  latent <- matrix(0,n,3)
  for(j in 1:3) latent[,j] <- c(x)*c(.45,-.35,.65)[j]+c(w[,,j])+eps[,j]
  g <- exp(latent)
  C <- runif(n,0,cu)
  y <- lo <- ex <- ce <- matrix(0,n,3)
  entered <- rep(TRUE,n)
  for(j in 1:3) {
    event <- entered & (g[,j]<=C)
    cens <- entered & !event
    y[event,j] <- latent[event,j]
    lo[cens,j] <- log(pmax(C[cens],1e-15))
    ex[event,j] <- 1; ce[cens,j] <- 1
    C <- pmax(C-g[,j],0)
    entered <- event
  }
  lapply(1:B,function(b) {
    ix <- ((b-1)*m+1):(b*m)
    list(x=c(x)[ix],y=y[ix,,drop=FALSE],lo=lo[ix,,drop=FALSE],
         ex=ex[ix,,drop=FALSE],ce=ce[ix,,drop=FALSE])
  })
}
# E[min(T/cu,1)] = target censoring fraction.
calibrate <- function(tau,target) {
  set.seed(991773+round(100*tau)+round(1000*target))
  N <- 100000
  x <- rnorm(N)
  eps <- rmvnorm(N,sigma=Rvine(atanh(rtrue)))*.85
  w <- matrix(rnorm(3*N,sd=tau),N,3)
  T <- rowSums(exp(outer(x,c(.45,-.35,.65))+eps+w))
  uniroot(function(cu) mean(pmin(T/cu,1))-target,c(.01,2000))$root
}
natural <- function(th,z) {
  P <- matrix(th,4,3)
  c(P[1,],P[2,],exp(P[3,]),exp(P[4,]),tanh(z))
}
one <- function(n,tau,cu) {
  dat <- generate(n,tau,cu)
  l1 <- function(th) -sum(vapply(dat,function(d) block1(th,d),0.0))
  th0 <- unlist(lapply(1:3,function(j)c(0,c(.45,-.35,.65)[j],log(.85),log(tau))))
  fit1 <- optim(th0,l1,method="L-BFGS-B",
                lower=rep(c(-5,-4,-2.3,-2.3),3),
                upper=rep(c(5,4,1.7,1.7),3),
                control=list(maxit=160))
  th <- fit1$par
  l2 <- function(z,th=th) -sum(vapply(dat,function(d)block2(th,z,d),0.0))
  fit2 <- optim(atanh(.7*rtrue),l2,method="BFGS",control=list(maxit=110))
  z <- fit2$par
  # Cluster stacked Murphy--Topel covariance; numerical differentiation.
  s1 <- t(vapply(dat,function(d) grad(function(a)block1(a,d),th),numeric(12)))
  s2 <- t(vapply(dat,function(d) grad(function(a)block2(th,a,d),z),numeric(3)))
  H11 <- hessian(l1,th)
  H22 <- hessian(function(a)l2(a,th),z)
  H21 <- jacobian(function(a)grad(function(v)l2(v,a),z),th)
  A <- rbind(cbind(H11,matrix(0,12,3)),cbind(H21,H22))
  S <- cbind(s1,s2)
  V <- solve(A,crossprod(S))%*%t(solve(A))
  J <- matrix(0,15,15)
  for(j in 1:3) {
    J[j,4*(j-1)+1] <- 1; J[j+3,4*(j-1)+2] <- 1
    J[j+6,4*(j-1)+3] <- exp(th[4*(j-1)+3])
    J[j+9,4*(j-1)+4] <- exp(th[4*(j-1)+4])
    J[j+12,j+12] <- 1-tanh(z[j])^2
  }
  se <- sqrt(pmax(diag(J%*%V%*%t(J)),0))
  est <- natural(th,z)
  valid <- all(is.finite(c(est,se))) && all(se>0) &&
           all(abs(z)<4.5) && kappa(V)<1e14
  list(est=est,se=se,valid=valid)
}
simulate <- function(n=1200,tau=.5,censor=.2,reps=500,cores=4) {
  cu <- calibrate(tau,censor)
  RNGkind("L'Ecuyer-CMRG"); set.seed(20261102)
  if (.Platform$OS.type=="windows") {
    cl <- makeCluster(cores); on.exit(stopCluster(cl))
    clusterExport(cl,c("m","K","rtrue","Rvine","lcens","covmat","block1",
       "block2","generate","natural","one"),envir=environment())
    clusterEvalQ(cl,{library(mvtnorm);library(numDeriv)})
    clusterSetRNGStream(cl,20261102)
    out <- parLapply(cl,seq_len(reps),function(i)
      tryCatch(one(n,tau,cu),error=function(e)NULL))
  } else {
    out <- mclapply(seq_len(reps),function(i)
      tryCatch(one(n,tau,cu),error=function(e)NULL),
      mc.cores=cores,mc.set.seed=TRUE)
  }
  ok <- Filter(function(a)!is.null(a)&&isTRUE(a$valid),out)
  truth <- c(rep(0,3),c(.45,-.35,.65),rep(.85,3),rep(tau,3),rtrue)
  if(length(ok)<2) stop("Fewer than two usable replications")
  E <- t(vapply(ok,`[[`,numeric(15),"est"))
  SE <- t(vapply(ok,`[[`,numeric(15),"se"))
  data.frame(n=n,tau=tau,censor=censor,attempted=reps,usable=length(ok),
             parameter=labels,true=truth,
             bias=colMeans(sweep(E,2,truth)),
             empirical_sd=apply(E,2,sd),
             mean_MT_se=colMeans(SE),
             coverage_95=colMeans(abs(sweep(E,2,truth))<=1.96*SE))
}
# Example: 4 scenarios x 500 replications. Adjust cores to your machine.
if (sys.nframe()==0L) {
  grid <- expand.grid(tau=c(.5,.9),censor=c(.2,.4))
  ans <- do.call(rbind,lapply(seq_len(nrow(grid)),function(k)
    simulate(n=1200,tau=grid$tau[k],censor=grid$censor[k],
             reps=500,cores=max(1L,detectCores()-1L))))
  write.csv(ans,"spatial_dvine_R_results.csv",row.names=FALSE)
}
