ar1=function(n,rho=0.9) rho^toeplitz(0:(n-1))
pos=function(x) ifelse(x<0,0,x)
# penalized estimation of number of nonzero means
penfun=function(Z,LD,ngwas) {
  m=nrow(Z);p=ncol(Z)
  counter=0
  Ip=diag(p)
  pens=c()
  ds=list()
  for(k in 0:p) {
    cs=combn(1:p,k)
    for(j in 1:ncol(cs)) {
      counter=counter+1
      d=rep(0,p);d[cs[,j]]=1;D=diag(c(d))
      like=matrixNormal::dmatnorm(Z,Z%*%D,LD,Ip,log=TRUE)/m
      reg=sum(d)
      pen=reg-like
      pens=c(pens,pen)
      ds[[counter]]=d
    }
  }
  # plot(pens)
  dix=ds[[which.min(pens)]]
  list(pens=pens,d_minpen=dix)
}
# soft thresholding
S=function(x,lambda) sapply(x,function(h) sign(h)*max(c(0,abs(h)-lambda)))
# penalized mean estimation
penmu=function(z,LD,upper=qnorm(0.975),reg_constant=1) {
  lams=seq(0,upper,0.05)
  pens=likes=regs=c()
  for(i in 1:length(lams)) {
    regz=S(z,lams[i])
    like=dmvn(z,regz,LD,log=TRUE)
    reg=reg_constant*sum(regz!=0)
    pens[i]=reg-like
    likes[i]=like
    regs[i]=reg
  }
  S(z,lams[which.min(pens)])
}
# matrix trace
tr=function(x)sum(diag(x))
# find parameters of gamma distribution
gampars=function(z,LD,null=TRUE) {
  if(null) {
    mu=length(c(z))
    sigma=2*tr(LD%*%LD)
  } else {
    muv=penmu(z,LD)
    mu=m+sum(muv^2)
    sigma=2*tr(LD%*%LD)+4*t(muv)%*%LD%*%muv
  }
  rate=mu/sigma
  shape=mu^2/sigma
  list(shape=shape,rate=rate)  
}
# deconvolute gamma mixture
sf=function(x) pos(2*(x-0.5))
mixp=function(z,LD) {
  alt_pars=apply(z,2,function(h) gampars(h,LD,null=FALSE));alt_pars=do.call(rbind,alt_pars)
  null_pars=apply(z,2,function(h) gampars(h,LD,null=TRUE));null_pars=do.call(rbind,null_pars)
  alt_pars=matrix(c(unlist(alt_pars)),nr=p,nc=2);colnames(alt_pars)=c('shape','rate')
  null_pars=matrix(c(unlist(null_pars)),nr=p,nc=2);colnames(null_pars)=c('shape','rate')
  stats=colSums(z^2)
  p_alts=c()
  for(i in 1:nrow(alt_pars)) {
    p1=dgamma(stats[i],shape=alt_pars[i,1],rate=alt_pars[i,2])
    p0=dgamma(stats[i],shape=null_pars[i,1],rate=null_pars[i,2])
    p_alts[i]=p1/(p1+p0)
  }
  p_alts=sf(p_alts) # shrinkage
  pboth=prod(p_alts) # P(mu1!=0,mu2!=0)
  peither=sum(p_alts)-prod(p_alts) # P(mu1!=0 or mu2!=0)
  pone=peither-pboth # P(mu1!=0 or mu2!=0 and !(mu1!=0 and mu2!=0))
  pone # P(only one is associated)
}
# function to generate simulation data
datgen=function(m,p,pcausal,LD,h2,ngwas_min,ngwas_max,niter) {
  ngwas=seq(ngwas_min,ngwas_max,length.out=p); ngwas=floor(ngwas)
  Th=solve(LD)
  zlist=b0list=row_ixs=col_ixs=list()
  for(iter in 1:niter) {
    row_ix=sample(1:m,1); row_ixs[[iter]]=row_ix
    col_ix=sample(1:p,pcausal); col_ixs[[iter]]=col_ix
    b0=matrix(0,m,p); b0list[[iter]]=b0
    b0[row_ix,col_ix]=sqrt(h2)
    z=matrix(nr=m,nc=p)
    for(i in 1:p) {
      bhat=b0[,i]+c(rmvn(1,rep(0,m),Th/ngwas[i]))
      betahat=LD%*%bhat
      shat=sqrt(rchisq(m,ngwas[i]-1)/(ngwas[i]-1)/ngwas[i])
      z[,i]=betahat/shat
    }
    zlist[[iter]]=z
  }
  list(zlist=zlist,b0list=b0list,row_ixs=row_ixs,col_ixs=col_ixs)
}