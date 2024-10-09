# AR1
ar1=function(n,rho=0.9) rho^toeplitz(0:(n-1))
# CS
cs=function(n,rho=0.9) {r=matrix(rho,n,n);diag(r)=1;r}
# trace function
tr=function(x) sum(diag(x))
# function to return P-values for sum(zs^2) using method of moments to find null distribution
gent=function(zs=NULL,LD,A=NULL,chisquares=NULL) {
  # find null distribution
  if(is.null(A)) {
    mu=nrow(LD)
    trASAS=tr(LD%*%LD)
  } else {
    A=as.matrix(A)
    mu=sum(diag(A%*%LD))
    trASAS=tr(A%*%LD%*%A%*%LD)
  }
  sigma2=2*trASAS
  beta=(mu/trASAS)/2
  alpha=beta*mu
  # P-values for zs (should be in LD)
  if(!is.null(chisquares)) y=sum(chisquares) else y=sum(zs^2)
  mu_h1=mu+y
  if(!is.null(zs)) sigma2_h1=sigma2+4*t(zs)%*%LD%*%zs else sigma2_h1=NULL
  pval=pgamma(y,shape=alpha,rate=beta,lower.tail=FALSE)
  list(pval=pval,shape=alpha,rate=beta,mu_h0=mu,sigma2_h0=sigma2,mu_h1=mu_h1,sigma2_h1=sigma2_h1)
}
linehist=function(x,lwd=1,linecol='black',...) {
  h=hist(x,border=NA,col=NA,pr=T,...)
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
# transformd statistic F

f=function(z,R) {
  n=length(z)
  d=eigen(R)
  U=d$vectors
  D=diag(d$values)
  ei=d$values
  vi=2*ei^2
  c0i=ei/vi
  y=t(U)%*%z
  fi=c0i*y^2
  fis=sum(fi)
  pval=pgamma(fis,shape=n/2,rate=1,lower.tail=FALSE)
  list(pval=pval,stat=fis)
}

