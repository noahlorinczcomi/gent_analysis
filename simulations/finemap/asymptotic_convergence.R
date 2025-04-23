rm(list=ls(all=TRUE))
library(mvnfast)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
ar1=function(n,rho) rho^toeplitz(0:(n-1))
atransform=function(Q,null_mean,null_variance) {
  # function to transform gene-based test statistics to an asymptotically normal distribution
  m=null_mean
  z=Q/m-1
  z*m/sqrt(null_variance)
}
tr=function(x) sum(diag(x))
linehist=function(x,col,lwd=1,linecol='black',fill=NA,...) {
  h=hist(x,border=NA,col=fill,pr=T,...)
  w=h$breaks[2]-h$breaks[1]
  n=length(h$breaks)
  # Draw top horizontal lines
  for(i in 1:(n-1)) lines(x=c(h$mids[i]-w/2,h$mids[i]+w/2),y=rep(h$density[i],2),lwd=lwd,col=linecol)
  # Draw vertical lines at the start of each bar
  for(i in 2:n) lines(x=rep(h$breaks[i],2),y=c(h$density[i-1],h$density[i]),lwd=lwd,col=linecol)
  # Draw the left-most and right-most vertical lines
  lines(x=rep(h$breaks[1],2),y=c(0,h$density[1]),lwd=lwd,col=linecol)
  lines(x=rep(tail(h$breaks,1),2),y=c(0,tail(h$density,1)),lwd=lwd,col=linecol)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
ms=c(50,100,500,1000)
niter=10000
par(mfrow=c(1,length(ms)),mar=c(4,2,2,2))
for(i in 1:length(ms)) {
  R=ar1(ms[i],0.5)
  e=ms[i]
  v=2*tr(R%*%R)
  Z=mvnfast::rmvn(niter,rep(0,ms[i]),R)
  Q=rowSums(Z^2)
  Zt=atransform(Q,e,v)
  linehist(Zt,breaks=50,xlab='transformed statistic',main=paste0(ms[i],' SNPs'),xlim=c(-5,5))
  curve(dnorm(x),col='red',add=T)
}
par(mfrow=c(1,1))
