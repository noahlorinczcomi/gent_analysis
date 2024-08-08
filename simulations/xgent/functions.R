tr=function(x)sum(diag(x))
ar1=function(n,rho) rho^toeplitz(0:(n-1))
xvegas=function(z,R,Z_xqtls) {
  z=as.matrix(z);Z_xqtls=as.matrix(Z_xqtls)
  m=nrow(R);p=ncol(Z_xqtls)
  L=matrix(0,nrow=nrow(Z_xqtls),ncol=nrow(Z_xqtls))
  for(i in 1:p) L=L+Z_xqtls[,i]%*%t(Z_xqtls[,i])
  L=L/sqrt(m*p)
  mu=tr(R%*%R)
  variance=2*tr(L%*%R%*%L%*%R)
  beta=mu/variance
  alpha=mu*beta
  stat=t(z)%*%L%*%z
  pval=pgamma(c(stat),shape=alpha,rate=beta,lower.tail=FALSE)
  list(pval=pval,stat=stat,alpha=alpha,beta=beta,mu=mu,variance=variance)
}
vegas=function(z,R) {
  z=as.matrix(z);
  m=nrow(R)
  beta=m/tr(R%*%R)/2
  alpha=m*beta
  stat=t(z)%*%z
  pgamma(c(stat),shape=alpha,rate=beta,lower.tail=FALSE)
}
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