library(survival)
colonc=colon
colonc1=colonc[colonc[, "etype"]==1,]
colonc2=colonc[colonc[, "etype"]==2,]
rm_idx=which(colonc1[,"status"]==0&colonc2[,"status"]==1)
colonc1_clean=colonc1[-rm_idx, ,drop=F]
colonc2_clean=colonc2[-rm_idx, ,drop=F]
colonc2_clean$time=colonc2_clean$time-colonc1_clean$time
colonc1_clean$lev=ifelse(colonc1_clean$rx == "Lev", 1, 0)
colonc1_clean$lev5fu=ifelse(colonc1_clean$rx == "Lev+5FU", 1, 0)
colonc2_clean$lev=ifelse(colonc2_clean$rx == "Lev", 1, 0)
colonc2_clean$lev5fu=ifelse(colonc2_clean$rx == "Lev+5FU", 1, 0)
na_rows1=which(!complete.cases(colonc1_clean))
na_rows2=which(!complete.cases(colonc2_clean))
na_rows=union(na_rows1, na_rows2)
colonc1_clean=colonc1_clean[-na_rows,]
colonc2_clean=colonc2_clean[-na_rows,]
idx_bad=which(colonc2_clean$time==0&colonc2_clean$status==1)
colonc2_clean$time[idx_bad]=10
n=851
colonc2_clean_pos=subset(colonc2_clean,time>0)
idx=which(colonc2_clean$time>0)
colonc1_clean_pos=colonc1_clean[idx,]
m=444
rm(colonc,colonc1,colonc2)
rm(idx,idx_bad,na_rows,na_rows1,na_rows2,rm_idx)

#####################################################


library(survival)
library(numDeriv)

xvars=c("lev","lev5fu","sex","age","obstruct","perfor","adhere",
        "nodes","differ","extent","surg","node4")

form.x=as.formula(paste0("~",paste(xvars,collapse="+")))
form.surv=as.formula(paste0("Surv(time,status)~",paste(xvars,collapse="+")))

X=model.matrix(form.x,data=colonc1_clean)
q=ncol(X)
n=nrow(colonc1_clean)

y1=colonc1_clean$time
d1=colonc1_clean$status
y2=colonc2_clean$time
d2=colonc2_clean$status

id1=which(d1==1)

f1=survreg(form.surv,data=colonc1_clean,dist="loglogistic")
f2=survreg(form.surv,data=colonc2_clean_pos,dist="loglogistic")

########################################################################
# parameter locations
#
# 1                  = transformed p
# 2:(q+1)            = beta1
# q+2                = log sigma1
# (q+3):(2q+2)       = beta2
# 2q+3               = log sigma2
# (2q+4):(2q+6)      = gamma0,gamma1,gamma2
########################################################################

i.p=1
i.b1=2:(q+1)
i.s1=q+2
i.b2=(q+3):(2*q+2)
i.s2=2*q+3
i.g=(2*q+4):(2*q+6)

start=c(qlogis(.60),
        coef(f1),log(f1$scale),
        coef(f2),log(f2$scale),
        log(1),0,0)

########################################################################
# log-logistic margin
########################################################################

margin=function(t,mu,sigma){
  z=(log(t)-mu)/sigma
  F=plogis(z)
  logf=dlogis(z,log=T)-log(sigma)-log(t)
  list(F=F,logf=logf)
}

########################################################################
# Clayton quantities
########################################################################

clayton=function(u,v,phi){
  A=u^(-phi)+v^(-phi)-1
  C1=u^(-phi-1)*A^(-1/phi-1)
  logc=log(1+phi)+(-1-phi)*(log(u)+log(v))+
    (-2-1/phi)*log(A)
  list(C1=C1,logc=logc)
}

########################################################################
# individual negative log-likelihood
########################################################################

loglik.i=function(par){
  
  p=plogis(par[i.p])
  
  beta1=par[i.b1]
  sigma1=exp(par[i.s1])
  
  beta2=par[i.b2]
  sigma2=exp(par[i.s2])
  
  gamma0=par[i.g[1]]
  gamma1=par[i.g[2]]
  gamma2=par[i.g[3]]
  
  mu1=as.vector(X%*%beta1)
  M1=margin(y1,mu1,sigma1)
  
  u1=M1$F
  logf1=M1$logf
  
  answer=rep(0,n)
  
  ########################################################################
  # first event censored:
  # 1-pF0(t1)
  ########################################################################
  
  id0=which(d1==0)
  answer[id0]=-log(1-p*u1[id0])
  
  ########################################################################
  # first event observed
  ########################################################################
  
  X1=X[id1,,drop=F]
  mu2=as.vector(X1%*%beta2)
  
  M2=margin(y2[id1],mu2,sigma2)
  
  u2=M2$F
  logf2=M2$logf
  
  phi=exp(gamma0+
            gamma1*colonc1_clean$lev[id1]+
            gamma2*colonc1_clean$lev5fu[id1])
  
  CP=clayton(u1[id1],u2,phi)
  
  ########################################################################
  # first event observed, second censored:
  # p f0(t1){1-C1(u1,u2)}
  ########################################################################
  
  j0=which(d2[id1]==0)
  
  answer[id1[j0]]=-(log(p)+
                      logf1[id1[j0]]+
                      log(1-CP$C1[j0]))
  
  ########################################################################
  # both events observed:
  # p f0(t1)f2(t2)c_phi(u1,u2)
  ########################################################################
  
  j1=which(d2[id1]==1)
  
  answer[id1[j1]]=-(log(p)+
                      logf1[id1[j1]]+
                      logf2[j1]+
                      CP$logc[j1])
  
  answer
}

########################################################################
# full negative log-likelihood
########################################################################

loglik=function(par){
  sum(loglik.i(par))
}

########################################################################
# one-stage maximum likelihood
########################################################################

result=nlm(loglik,start,hessian=T,iterlim=1000)
result

est=result$estimate

########################################################################
# estimates
########################################################################

p.hat=plogis(est[i.p])

beta1.hat=est[i.b1]
sigma1.hat=exp(est[i.s1])

beta2.hat=est[i.b2]
sigma2.hat=exp(est[i.s2])

gamma.hat=est[i.g]

names(beta1.hat)=colnames(X)
names(beta2.hat)=colnames(X)
names(gamma.hat)=c("gamma0","gamma1","gamma2")

########################################################################
# sandwich estimator
#
# V=H^{-1}BH^{-1}
########################################################################

H=result$hessian

score=jacobian(loglik.i,est)

B=crossprod(score)

H.inv=solve(H)

V.eta=H.inv%*%B%*%t(H.inv)

########################################################################
# delta method to natural parameter scale
########################################################################

D=diag(length(est))

D[i.p,i.p]=p.hat*(1-p.hat)
D[i.s1,i.s1]=sigma1.hat
D[i.s2,i.s2]=sigma2.hat

V=D%*%V.eta%*%t(D)
se=sqrt(diag(V))

########################################################################
# parameter names
########################################################################

par.names=c("p",
            paste0("T1_",colnames(X)),
            "T1_scale",
            paste0("T2_",colnames(X)),
            "T2_scale",
            "gamma0","gamma1","gamma2")

names(se)=par.names

estimate=c(p.hat,
           beta1.hat,
           sigma1.hat,
           beta2.hat,
           sigma2.hat,
           gamma.hat)

names(estimate)=par.names

all.results=data.frame(
  parameter=par.names,
  estimate=estimate,
  sandwich.se=se,
  lower=estimate-1.96*se,
  upper=estimate+1.96*se
)

all.results

########################################################################
# 95% two-sided Wald tests for treatment-dependent association
########################################################################

g.id=(2*q+4):(2*q+6)

g.est=est[g.id]
g.se=sqrt(diag(V.eta[g.id,g.id]))

z=g.est/g.se
pvalue=2*pnorm(-abs(z))

dep.test=data.frame(
  parameter=c("gamma0","gamma1: Lev vs placebo","gamma2: Lev+5FU vs placebo"),
  estimate=g.est,
  sandwich.se=g.se,
  z=z,
  p.value=pvalue,
  lower=g.est-1.96*g.se,
  upper=g.est+1.96*g.se,
  significant=pvalue<.05
)

dep.test

########################################################################
# treatment-specific Clayton dependence
########################################################################

phi0=exp(g.est[1])
phi1=exp(g.est[1]+g.est[2])
phi2=exp(g.est[1]+g.est[3])

tau0=phi0/(phi0+2)
tau1=phi1/(phi1+2)
tau2=phi2/(phi2+2)

dependence=data.frame(
  group=c("Placebo","Lev","Lev+5FU"),
  phi=c(phi0,phi1,phi2),
  tau=c(tau0,tau1,tau2)
)

dependence
