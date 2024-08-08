ar1=function(n,rho=0.9) rho^toeplitz(0:(n-1))
tr=function(x) sum(diag(x))
anova_mugent=function(chisquares,LDlist,p) {
  # does NOT all LD matrices are the same across populations
  chisquares=c(chisquares)
  m=length(chisquares)
  R1=matrix(0,nr=m,nc=m)
  for(ll in 1:p) R1=R1+LDlist[[ll]]^2; R1=R1/p
  D=diag(sqrt(2*p),m)
  j=rep(1,m)
  mu=m*p
  variance=c(t(j)%*%D%*%R1%*%D%*%j) # actually same as sum(D%*%R1%*%D)
  beta=mu/variance
  alpha=beta*mu
  p=pgamma(sum(chisquares),shape=alpha,rate=beta,lower.tail=FALSE)
  out=list(null_mean=mu,null_variance=variance,gamma_alpha=alpha,gamma_beta=beta,pvalue=p)
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