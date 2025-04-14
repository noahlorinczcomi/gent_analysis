#rm(list=ls(all=TRUE))
library(data.table);library(magrittr);library(tidyr);library(dplyr);library(ggplot2)
library(mvnfast,lib='/home/lorincn/Rpkgs')
library(mvsusieR,lib='/home/lorincn/Rpkgs')
# library(snpsettest,lib='/home/lorincn/Rpkgs')
# source('/home/lorincn/Rpkgs/manual_snpsettestcode.R')
library(ACAT,lib='/home/lorincn/Rpkgs')
library(gent,lib='/home/lorincn/Rpkgs')
library(exset,lib='/home/lorincn/Rpkgs')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# gene-based association tests ####
# snpsettest (saddlepoint/numerical approximation to the null distribution of the VEGAS statistic)
#  - Method: https://github.com/HimesGroup/snpsettest/wiki/Statistical-test-in-snpsettest
#  - Code: https://github.com/HimesGroup/snpsettest
#  - Example usage: snpsettest(z^2,ld_eigenvalues)
# ACAT (SKAT with summary statistics using Cauchy combination test)
#  - Method: https://www.sciencedirect.com/science/article/pii/S0002929719300023?via%3Dihub 
#  - Code: https://github.com/yaowuliu/ACAT
#  - Example usage: ACAT(pvals)
# GATES
#  - Method: https://pmc.ncbi.nlm.nih.gov/articles/PMC3059433/
#  - Code: GATES() (see below)
#  - Example usage: GATES(pvals)
# GenT
#  - Method: https://dx.doi.org/10.2139/ssrn.5080346 
#  - Code: https://github.com/noahlorinczcomi/gent
#  - Example usage: gent(zs,ld)
# exset (exact set-based testing)
#  - Method: https://github.com/noahlorinczcomi/exset 
#  - Code: https://github.com/noahlorinczcomi/exset 
#  - Example usage: exset(sum(zs^2),ld_eigenvalues)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# functions ####
theme_set(
  theme_bw() +
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          strip.background=element_rect(fill="beige"),
          panel.background=element_rect(fill='#EBF9FF',colour='black'))
)
tr=function(x) sum(diag(x))
ar1=function(n,rho) rho^toeplitz(0:(n-1))
CS=function(n,rho) matrix(rho,n,n)+(1-rho)*diag(n)
posDefifyCorrMat=function(ndmatrix,epsilon=1e-4) {
  # ndmatrix: pxp square negative-definite matrix of measurement error correlations/covariances
  # since it is a correlation matrix, abs max off-diagonal value is 1
  deco=eigen(ndmatrix)
  eig0=min(deco$values)
  eigP=max(deco$values)
  if(eig0<epsilon) {
    mu=max(c(epsilon,(eig0+eigP)/2))
    alpha=(mu-epsilon)/(mu-eig0) # optimal choice given by Choi et al: https://doi.org/10.1016/j.jmva.2018.12.002
    Sigmahat=alpha*ndmatrix+(1-alpha)*mu*diag(nrow(ndmatrix)) # the solution
    deco=eigen(Sigmahat)
  } else {
    alpha=1
    Sigmahat=ndmatrix
  }
  return(list(Sigmahat=Sigmahat,alpha=alpha,deco=deco))
}
pchisqsum <- function (x, df, a, lower.tail = TRUE) {
  sat <- satterthwaite(a, df)
  guess <- pchisq(x / sat$scale, sat$df, lower.tail = lower.tail)
  for (i in seq(length = length(x))) {
    lambda <- rep(a, df)
    sad <- sapply(x, saddle, lambda = lambda)
    if (lower.tail) sad <- 1 - sad
    guess <- ifelse(is.na(sad), guess, sad)
  }
  return(guess)
}
satterthwaite <- function(a, df) {
  if (any(df > 1)) {
    a <- rep(a, df)
  }
  tr <- mean(a)
  tr2 <- mean(a^2) / (tr^2)
  list(scale = tr * tr2, df = length(a) / tr2)
}
saddle <- function(x, lambda) {
    d <- max(lambda)
    lambda <- lambda / d
    x <- x / d
    k0 <- function(zeta) {
      -sum(log(1 - 2 * zeta * lambda)) / 2
    }
    kprime0 <- function(zeta) {
      sapply(zeta, function(zz) sum(lambda / (1 - 2 * zz * lambda)))
    }
    kpprime0 <- function(zeta) {
      2 * sum(lambda^2 / (1 - 2 * zeta * lambda)^2)
    }
    if (any(lambda < 0)) {
        lmin <- max(1 / (2 * lambda[lambda < 0])) * 0.99999
    } else if (x > sum(lambda)) {
        lmin <- -0.01
    } else {
        lmin <- -length(lambda) / (2 * x)
    }
    lmax <- min(1 / (2 * lambda[lambda > 0])) * 0.99999
    hatzeta <- uniroot(function(zeta) kprime0(zeta) - x, lower = lmin,
                       upper = lmax, tol = 1e-08)$root
    w <- sign(hatzeta) * sqrt(2 * (hatzeta * x - k0(hatzeta)))
    v <- hatzeta * sqrt(kpprime0(hatzeta))
    if (abs(hatzeta) < 1e-04) {
      NA
    }  else {
      pnorm(w + log(v / w) / w, lower.tail = FALSE)
    }
}
snpsettest=function(chisquares,ld_eigenvalues) pchisqsum(sum(chisquares),df=1,ld_eigenvalues)
meff=function(lambda,M=length(lambda)) sum(+(lambda>=1)+lambda*(lambda<1))
GATES=function(pvalues,ld) {
  M=length(pvalues)
  ix=order(pvalues)
  p=pvalues[ix]
  R=ld[ix,ix]
  mej=c()
  for(j in 1:M) {
    Rj=R[1:j,1:j]
    mej[j]=meff(eigen(Rj)$values)
  }
  me=tail(mej,1)
  Pg=min(me*p/mej)
  Pg
}
simdata=function(m,rho,z0=rep(0,m),nref=505,niter=1000) {
  R0=ar1(m,rho)
  Z=mvnfast::rmvn(niter,z0,R0)
  R=rWishart(niter,max(c(m+1,nref)),R0)
  for(i in 1:niter) R[,,i]=cov2cor(R[,,i])
  list(Q=rowSums(Z),Z=Z,R=R,R0=R0)
}
simQs=function(es,vs,R,ngenes=ncol(R),niter=1000) {
  ## R: correlation between GenT statistics
  # normal copula
  rates=es/vs
  shapes=es*rates
  X=mvnfast::rmvn(niter,rep(0,ngenes),R)
  p=pnorm(X)
  Qs=matrix(nr=niter,nc=ngenes)
  for(i in 1:ngenes) Qs[,i]=qgamma(p[,i],shape=shapes[i],rate=rates[i])
  Qs
}
atransform=function(Q,null_mean,null_variance) {
  # function to transform gene-based test statistics to an asymptotically normal distribution
  m=null_mean
  z=Q/m-1
  z*m/sqrt(null_variance)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# simulation comparing GenT, exset, ACACT, and VEGAS/snpsettest (numerical approximation) ####
## Type I error
ms=c(10,100,500)
rhos=c(0.1,0.5,0.9)
niter=1000
RES1=RES2=RES3=RES4=RES5=matrix(nr=length(ms),nc=length(rhos))
T1=T2=T3=T4=T5=RES1
for(i in 1:length(ms)) {
  m=ms[i]
  cat('m=',m,'\n',sep='')
  for(j in 1:length(rhos)) {
    rho=rhos[j]
    cat(' rho=',rho,'\n',sep='')
    LD0=ar1(m,rho)
    LD=rWishart(niter,505,LD0) # 505 is the size of the European 1000 Genomes Phase 3 cohort
    for(o in 1:niter) LD[,,o]=cov2cor(LD[,,o])
    Z=mvnfast::rmvn(niter,rep(0,m),LD0)
    Z2=Z^2
    P=pchisq(Z2,1,lower.tail=FALSE)
    res1=res2=res3=res4=res5=c() # H0 rejection indicator
    t1=t2=t3=t4=t5=c() # simulation time
    for(iter in 1:niter) {
      if(iter%%(niter*0.1)==0) cat('  ',round(iter/niter*100),'% complete\n',sep='')
      LDiter=LD[,,iter]
      # ld_eigenvalues=eigen(LDiter)$values # using estimated LD
      ld_eigenvalues=eigen(LD0)$values # using true LD
      # GenT
      t0=Sys.time()
      res1[iter]=gent(Z[iter,],LDiter)$pval<0.05
      t1[iter]=as.numeric(difftime(Sys.time(),t0,units="secs"))
      # exset
      t0=Sys.time()
      res2[iter]=exset(sum(Z2[iter,]),ld_eigenvalues,iterations=1e6,type1_error=0.05)=='reject H0'
      t2[iter]=as.numeric(difftime(Sys.time(),t0,units="secs"))
      # ACAT
      t0=Sys.time()
      res3[iter]=ACAT(P[iter,])<0.05
      t3[iter]=as.numeric(difftime(Sys.time(),t0,units="secs"))
      # snpsettest (VEGAS)
      t0=Sys.time()
      res4[iter]=snpsettest(Z2[iter,],ld_eigenvalues)<0.05
      t4[iter]=as.numeric(difftime(Sys.time(),t0,units="secs"))
      # GATES
      t0=Sys.time()
      res5[iter]=GATES(P[iter,],LDiter)<0.05
      t5[iter]=as.numeric(difftime(Sys.time(),t0,units="secs"))
    }
    RES1[i,j]=mean(res1)
    RES2[i,j]=mean(res2)
    RES3[i,j]=mean(res3)
    RES4[i,j]=mean(res4)
    RES5[i,j]=mean(res5)
    T1[i,j]=mean(t1)
    T2[i,j]=mean(t2)
    T3[i,j]=mean(t3)
    T4[i,j]=mean(t4)
    T5[i,j]=mean(t5)
  }
}
# organize
todf=function(RES,TX) data.frame(x=c(RES),time=c(TX),m=rep(ms,length(rhos)),rho=rep(rhos,each=length(ms)))
df=bind_rows(
  todf(RES1,T1) %>% mutate(method='GenT'), # Gamma approximation
  todf(RES2,T2) %>% mutate(method='exset'), # exact test
  todf(RES3,T3) %>% mutate(method='ACAT'), # Cauchy combination test
  todf(RES4,T4) %>% mutate(method='VEGAS'), # saddlepoint/Satterthwaite approximation
  todf(RES5,T5) %>% mutate(method='GATES') # adjusted SIMES method
)
# save df









