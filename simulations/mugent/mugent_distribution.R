rm(list=ls(all=TRUE))
library(mvnfast);library(mvnfast)
source('simulations/mugent/functions.R')
#################################################################################
niter=1000
ms=c(50,100,150)
ps=c(2,5,7)
rhos=c(0.9,0.1)
RES=array(dim=c(length(ms),length(ps),length(rhos)))
par(mfrow=c(6,3),mar=c(1.5,1.5,1.5,1.5))
for(i in 1:length(ms)) {
  z0=rep(0,ms[i])
  for(j in 1:length(ps)) {
    for(k in 1:length(rhos)) {
      # different LD matrices for each population
      ldlist=list();for(o in 1:ps[j]) ldlist[[o]]=ar1(ms[i],rhos[k]^o)
      stats=pvals=c()
      for(iter in 1:niter) {
        z=matrix(nr=ps[j],nc=ms[i]); 
        for(o in 1:ps[j]) z[o,]=rmvn(1,z0,ldlist[[o]])
        out=mugent(t(z),ldlist)
        stats[iter]=out$stat
        pvals[iter]=out$p
      }
      RES[i,j,k]=sum(pvals<0.05)/niter
      main.=paste0('# SNPs=',ms[i],'; LD rho=',rhos[k],'; ',ps[j],' pops.')
      linehist(stats,breaks=50,xlab='null test statistics',main=main.)
      curve(dgamma(x,shape=out$alpha,rate=out$beta),min(stats),max(stats),add=T,col='blue')
      legend('topright',cex=1/2.5,legend=c('empiricial density','theoretical density'),lty=c(1,1),col=c('black','blue'))
      x=density(stats)
      text(x=mean(stats),y=max(x$y)/10,labels=paste0('Type I error: ',sum(pvals<0.05)/niter*100,'%'),cex=2/3)
    }
  }
}












