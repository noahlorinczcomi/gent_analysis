rm(list=ls(all=TRUE))
library(mvnfast)
source('simulations/gent/functions.R')
##################################################################################
# changing ms, ngwas, LD density
ms=c(50,100,150)
ngwas=c(10000,50000,100000)
rhos=c(0.1,0.5,0.9)
h2=0
nref=500
niter=10000
#m.=n.=r.=1
for(m. in 1:length(ms)) {
  for(n. in 1:length(ngwas)) {
    for(r. in 1:length(rhos)) {
      # generate data under null
      R0=ar1(ms[m.],rhos[r.]) # true LD matrix
      Th0=solve(R0) # inverse of true LD matrix
      b=rep(0,ms[m.]) # under H0: b_j=0 for all j=1,...,ms[m.]
      bhat=lapply(1:niter,function(h) c(rmvn(1,b,Th0/ngwas[n.])))
      betahat=lapply(1:niter,function(h) c(R0%*%bhat[[h]]))
      betahat=do.call(rbind,betahat)
      s2hat=rchisq(prod(dim(betahat)),ngwas[n.]-1)/(ngwas[n.]-1)/ngwas[n.]
      zhat=betahat/sqrt(s2hat)
      stats=rowSums(zhat^2)
      # estimate null distribution using approximation
      R=rWishart(niter,nref,R0) # LD matrix estimates from reference panel
      R=lapply(1:niter,function(h) cov2cor(R[,,h]))
      trRRs=lapply(1:niter,function(h) 2*tr(R[[h]]%*%R[[h]]))
      trRRs=unlist(trRRs)
      xi=ms[m.]/mean(trRRs)
      alpha=ms[m.]*xi
      # plot
      main=paste0('# SNPs=',ms[m.],'; # ppl.=',ngwas[n.]/1e3,'K; LD rho=',rhos[r.])
      linehist(stats,linecol='black',breaks=30,main=main,xlab='simulated statistics')
      dd=density(stats)
      q0=qgamma(0.95,shape=alpha,rate=xi)
      lab.=paste0('Type I error: ',round(sum(stats>q0)/niter*100,2),'%')
      text(x=median(stats),median(dd$y)*0.2,label=lab.,cex=4/5)
      curve(dgamma(x,shape=alpha,rate=xi),min(stats),max(stats),add=T,col='blue')
      legend('topright',legend=c('empirical density','theoretical density'),lty=c(1,1),col=c('black','blue'),cex=2/3.5)
    }
  }
}
