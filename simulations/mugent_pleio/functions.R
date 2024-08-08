# gaussian (parametric) mixture of N(0,1) and N(0,sigma2)
pcausalsnp=function(stats,nullvar=1,mu1=0,sd1=1.0001,niter=30,verbose=T) {
  if(verbose) cat('Assuming `stats` are Z-statistics')
  m=length(stats)
  nullsd=sqrt(nullvar)
  d0=dnorm(stats,0,nullsd,log=T)
  d1=dnorm(stats,mu1,sd1,log=T)
  g=rep(0,m)
  g[d1>d0]=1
  k=0;eps=floor(m*0.01);error=eps+1
  while(k<niter & error>eps) {
    k=k+1
    xbar0=0 # fixed
    xbar1=0 # fixed
    sbar0=nullsd # fixed
    sbar1=mad(stats[g==1]) # not fixed
    d0=dnorm(stats,0,nullsd,log=T)
    d1=dnorm(stats,xbar0,sbar1,log=T)
    g0=+(d1>d0)
    error=sum(g!=g0)
    g=g0
  }
  p=dnorm(stats,xbar1,sbar1)/(dnorm(stats,xbar1,sbar1)+dnorm(stats,0,nullsd))
  if(abs(stats[which.max(p)])<abs(stats[which.min(p)])) p=max(p)-p
  if(min(p)>0.1) p=1-p
  p
}
# function to estimate P(SNP is associated with 3 traits)
multicausalsnp=function(Z,Rho=diag(ncol(Z)),nullvar=1) {
  # Z: a mx3 matrix of marginal Z-statistics for SNP->X (1st col), SNP->M (2nd col), SNP->Y (3rd col)
  # Rho: 3x3 matrix of correlations between GWAS estimation errors (identity matrix if all GWAS cohorts are independent)
  G=eigen(Rho)
  G=G$vectors%*%diag(1/sqrt(G$values))%*%t(G$vectors)
  Z=Z%*%G
  ## marginal causal probabilities (semi-parametric)
  D1=apply(Z,2,function(h) npdf(h,ks::hpi(h),kernel=KnGaussian))
  ## marginal causal probabilities (parametric)
  # D1=apply(Z,2,function(h) pcausalsnp(h,verbose=FALSE,nullvar=nullvar))
  D0=dnorm(Z,0,sqrt(nullvar))
  P=D1/(D1+D0)
  # joint probabilities
  pjoint=apply(P,1,prod)
  if(abs(Z[which.max(pjoint)])<abs(Z[which.min(pjoint)])) pjoint=max(pjoint)-pjoint
  if(min(pjoint)>0.1) pjoint=1-pjoint
  list(PH1=pjoint,marginalPs=P)
}
### new functions
# infers SNP is associated with all three traits adjusting for #independent SNPs
f=function(Z,LD,type1=0.05,reg=FALSE) {
  # Z: mxp matrix of Z-statistics (rows are SNPs, columns are traits)
  # LD: estimated matrix of LD correlations between m SNPs
  # type1: nominal Type I error rate which will be corrected
  Z=as.matrix(Z); LD=as.matrix(LD)
  m=nrow(Z);p=ncol(Z)
  G=eigen(LD);d=G$values
  meff=sum(d>=1+d*(d<1))
  # meff=m
  q_ac=sqrt(qchisq((1-type1)^(1/(p*meff)),1))
  checkf=apply(Z,1,function(h) all(abs(h)>q_ac))
  fvres=fv(q_ac,p,meff,type1)
  list(result=checkf,q=q_ac,nullvar=fvres)
  # result: TRUE/FALSE vector (TRUE means SNP associated with all traits; FALSE means the opposite)
  # q: critical value to infer a SNP is associated with all traits
  # nullvar: variance of multivariate null distribution which achieves level typeI at q_ac
}
# (internal) calculates null variance to achieve desired corrected Type I error rate
fv=function(q_ac,p,m,type1=0.05) {
  require(mvtnorm)
  # ac=1-(1-type1)^(1/(p*m))
  ac=type1
  I=diag(p)
  s=seq(1,10,length.out=100)
  d=c()
  for(i in 1:length(s)) {
    d[i]=abs(q_ac-qnorm(1-ac,0,sqrt(s[i])))
    # d[i]=abs(q_ac-qmvnorm(1-ac,mean=rep(0,p),sigma=s[i]*I)$quantile)
    if(i>2) {if(tail(d,1)>tail(d,2)[1] & tail(d,3)[1]>tail(d,2)[1]) break}
  }
  if(i==length(s)) warning('could not find null variance; try changing function')
  s[i-1] # variance of null distribution which satisfies F_0(q_ac)=1-ac with mean 0
}
# regularized LD (minimizes a penalty while maximizing a likelihood)
tr=function(x)sum(diag(x))
regld=function(Z,LD,nlambdas=20,doplot=FALSE) {
  m=nrow(LD)
  lams=seq(0.01,0.5,length.out=nlambdas)
  I=diag(m);Z0=Z*0;Ip=diag(ncol(as.matrix(Z)))
  pen=c()
  maxf=norm(LD-diag(m),'f');maxf=1
  maxd=2*sqrt(m)/log(sum(abs(LD)<lams[1])+2);maxd=1
  for(i in 1:length(lams)) {
    r=soft(LD,lams[i]); diag(r)=1
    s=sum(r==0)
    pen[i]=norm(LD-r,'f')+2*log(m^2-1+1)-matrixNormal::dmatnorm(Z,Z0,r,Ip,log=TRUE)
  }
  if(doplot) {plot(lams,pen,type='b',pch=4,cex=2/3,xlab=bquote(lambda),ylab='loss',col='red');abline(v=lams[which.min(pen)])}
  r=soft(LD,lams[which.min(pen)]); diag(r)=1
  if(min(eigen(r)$values)<0) r=posadj(r)
  r
}
mugent_pleio=function(Z,ldlist,alpha=0.05) {
  p=ncol(Z);meff=f_meff(ldlist)
  q0=qchisq(1-(1-(1-alpha)^(1/(2*meff)))^(1/p),1)
  boo=apply(Z^2,1,function(h) all(h>q0))
  list(result=any(boo),q0=q0)
}
meff_fun=function(A) {
  d=eigen(A)$values
  sum((d>=1)+d*(d<1))
}
simfun=function(tausqrt,mu,R,niter=1000) {
  z=rmvn(niter,mu,R)
  res=apply(z,1,function(h) all(abs(h)>tausqrt))
  sum(res)/niter
}
ar1=function(n,rho=0.9) rho^toeplitz(0:(n-1))