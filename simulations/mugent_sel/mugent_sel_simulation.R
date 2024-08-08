# rm(list=ls(all=TRUE))
# library(mvnfast);library(matrixNormal)
# 
# # Gamma likelihood
# S=function(x,lambda) sapply(x,function(h) sign(h)*max(c(0,abs(h)-lambda)))
# penmu=function(z,LD) {
#   lams=seq(0,qnorm(0.975),0.05)
#   pens=likes=regs=c()
#   for(i in 1:length(lams)) {
#     regz=S(z,lams[i])
#     like=dmvn(z,regz,LD,log=TRUE)
#     reg=sum(regz!=0)
#     pens[i]=reg-like
#     likes[i]=like
#     regs[i]=reg
#   }
#   # plot(likes)
#   # plot(regs)
#   # plot(pens)
#   S(z,lams[which.min(pens)])
# }
# gamlike=function(Z,LD,mu=NULL,D=NULL,log=FALSE) {
#   dens=c()
#   if(is.null(D)) D=matrix(1,nrow(Z),ncol(Z))
#   for(i in 1:ncol(Z)) {
#     if(!is.null(mu)) pmu=mu[,i] else pmu=penmu(Z[,i],LD)*D[,i]
#     mui=sum(diag(LD))+t(pmu)%*%pmu # hypothesized true mean
#     trSS=tr(LD%*%LD)
#     sigma2=2*trSS+4*t(pmu)%*%LD%*%pmu
#     beta=mui/sigma2
#     alpha=beta*mui
#     dens[i]=dgamma(sum(Z[,i]^2),shape=alpha,rate=beta,log=log) # likelihood 
#   }
#   dens
# }
# # Matrix normal likelihood
# matnormlike=function(Z,LD,mu,log=FALSE) matrixNormal::dmatnorm(Z,mu,LD,diag(ncol(Z)),log=log)
# # P(gene is associated in some populations but not all)
# ppartial_gene=function(Z,LD,prior_density=dbinom(0:ncol(Z),ncol(Z),0.1),alpha=1) {
#   # Z: mxp matrix of marginal SNP associations with gene in p populations
#   # LD: LD matrix between m SNPs in Z
#   # alpha: shrinkage parameter for likelihood estimation
#   Z=as.matrix(Z)
#   m=nrow(Z);p=ncol(Z)
#   ## same approach as for SNP but use Gamma likelihood
#   pd=prior_density
#   null=prod(gamlike(Z,LD,mu=Z*0))*pd[1]
#   full=prod(gamlike(Z,LD,mu=NULL,D=NULL))*tail(pd,1)
#   # null=matnormlike(Z,LD,Z*0)*pd[1]
#   # full=matnormlike(Z,LD,Z*alpha)*tail(pd,1)
#   partial=c()
#   k=0
#   for(s in 2:(p-1)) {
#     cs=combn(1:p,s)
#     for(j in 1:ncol(cs)) {
#       k=k+1
#       ixi=cs[,j];d=rep(0,p);d[ixi]=1;D=matrix(d,m,p,byrow=T)
#       partial[k]=prod(gamlike(Z,LD,mu=NULL,D=D))*pd[1+sum(d)]
#       # partial[k]=matnormlike(Z,LD,Z*alpha*D,log=FALSE)*pd[1+sum(d)]
#       names(partial)[k]=sum(d)
#     }
#   }
#   partial=sum(partial)
#   sumlike=null+full+partial
#   partial/sumlike
# }
# ar1=function(n,rho=0.9)rho^toeplitz(0:(n-1))
# m=100;p=5;pcausal=3
# LD=ar1(m)
# h2=0.0001
# ngwas=floor(seq(3e4,1e5,length.out=p))
# row_ix=sample(1:m,1)
# col_ix=sample(1:p,pcausal)
# z0=matrix(0,m,p)
# z0[row_ix,col_ix]=sqrt(ngwas[col_ix]*h2)
# niter=100
# r=c()
# for(iter in 1:niter) {
#   U=t(rmvn(p,rep(0,m),LD))
#   U=apply(U,2,function(h) (h-mean(h))/sd(h))
#   Z=z0+U
#   # res=penfun(Z,LD)
#   # r[iter]=sum(res$d_minpen)>0 & sum(res$d_minpen)<p
#   r[iter]=ppartial_gene(Z,LD,prior_density=dbinom(0:p,p,0.1),alpha=1)
#   if(iter%%round(floor(niter*0.1))==0) cat(round(iter/niter*100),'% complete\n',sep='')
# }
# plot(r,cex=1/3,ylim=c(0,1))
# table(r>0.9)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # function to generate simulation data
# datgen=function(m,p,pcausal,h2,ngwas_min,ngwas_max,niter) {
#   ngwas=seq(ngwas_min,ngwas_max,length.out=p); ngwas=floor(ngwas)
#   LD=ar1(m)
#   Th=solve(LD)
#   zlist=b0list=row_ixs=col_ixs=list()
#   for(iter in 1:niter) {
#     row_ix=sample(1:m,1); row_ixs[[iter]]=row_ix
#     col_ix=sample(1:p,pcausal); col_ixs[[iter]]=col_ix
#     b0=matrix(0,m,p); b0list[[iter]]=b0
#     b0[row_ix,col_ix]=sqrt(h2)
#     z=matrix(nr=m,nc=p)
#     for(i in 1:p) {
#       bhat=b0[,i]+c(rmvn(1,rep(0,m),Th/ngwas[i]))
#       betahat=LD%*%bhat
#       shat=sqrt(rchisq(m,ngwas[i]-1)/(ngwas[i]-1)/ngwas[i])
#       z[,i]=betahat/shat
#     }
#     zlist[[iter]]=z
#   }
#   list(zlist=zlist,b0list=b0list,row_ixs=row_ixs,col_ixs=col_ixs)
# }
# dg=datgen(m=100,p=5,pcausal=3,h2=0.0001,ngwas_min=3e4,ngwas_max=1e5,niter=1000)
# z=dg$zlist[[1]]
# 
# # prior_prob=0.2
# # nchain=5000;burnin=100
# # chain=list()
# # chain[[1]]=rep(0,p)
# # for(i in 2:nchain) {
# #   gamma0=chain[[i-1]]
# #   gamma1=rbinom(p,1,prior_prob)
# #   like0=c(); for(j in 1:p) like0[j]=dmvn(z[,j],penmu(z[,j],LD)*gamma0[j],LD)
# #   like1=c(); for(j in 1:p) like1[j]=dmvn(z[,j],penmu(z[,j],LD)*gamma1[j],LD)
# #   dprior0=dbinom(gamma0,1,prior_prob)
# #   dprior1=dbinom(gamma1,1,prior_prob)
# #   like0=prod(like0);like1=prod(like1)
# #   dprior0=prod(dprior0);dprior1=prod(dprior1)
# #   a=dprior1*like1/(dprior1*like1+dprior0*like0)
# #   if(a>runif(1)) {
# #     chain[[i]]=gamma1
# #   } else {
# #     chain[[i]]=gamma0
# #   }
# #   if(i%%round(floor(nchain*0.01))==0) cat(round(i/nchain*100),'% complete\n',sep='')
# # }
# chainmat=do.call(rbind,chain)
# chainmat=chainmat[-c(1:burnin),]
# barplot(colSums(chainmat)/(nchain-burnin),ylim=c(0,1));abline(h=0.9,lty=2)
# 
# 
# 
# 
# f=function(lamvec,z,LD) {
#   # z=apply(z,2,function(h) h/mad(h))
#   m=nrow(z);p=ncol(z);Ip=diag(p)
#   D=diag(lamvec)
#   regz=z%*%D
#   # regz=apply(z,2,function(h) penmu(h,LD))
#   like=dmatnorm(z,regz,LD,Ip,log=TRUE)/m
#   reg=sum(lamvec)
#   pen=reg-like-priord # penalty to minimize
#   pen
# }
# of=function(f,z,LD,method='L-BFGS',lower=rep(0,p),upper=rep(1,p)) {
#   optim(par=rep(0.1,p),f=f,z=z,LD=LD,method='L-BFGS',lower=rep(0,p),upper=rep(1,p))
# }
# 
# p=5;pcausal=4
# niter=100
# dg=datgen(m=100,p=p,pcausal=pcausal,h2=0.0005,ngwas_min=3e4,ngwas_max=1e5,niter=niter)
# truth_gene=guess_gene=matrix(0,nr=niter,nc=p)
# truth_inference=guess_inference=c()
# for(iter in 1:niter) {
#   zi=dg$zlist[[iter]]
#   res=of(f,zi,LD)
#   col_ix=dg$col_ixs[[iter]]
#   if(length(col_ix)>0) truth_gene[iter,dg$col_ixs[[iter]]]=1
#   par_ix=which(res$par!=0)
#   guess_gene[iter,par_ix]=1
#   guess_inference[iter]=ifelse(length(par_ix) %in% c(2:(p-1)),'distinct','non-distinct')
#   truth_inference[iter]=ifelse(pcausal>0 & pcausal<p,'distinct','non-distinct')
#   if(iter%%round(floor(niter*0.01))==0) cat(round(iter/niter*100),'% complete\n',sep='')
# }
# # barplot(colSums(truth_gene==guess_gene)/niter,ylim=c(0,1))
# rf=function() {
#   if(pcausal %in% c(0,p)) { # if truly non-distinct
#     cat('False positive rate:',sum(guess_inference=='distinct')/niter,'\n')
#   }
#   if(pcausal %in% c(1:(p-1))) { # if truly distinct
#     cat('True positive rate:',sum(guess_inference=='distinct')/niter,'\n')
#   }
# }
# rf()
# 
# ### how about assessing P(any population is associated) and P(all populations are associated) and doing 1-(their sum) ... this is P(some populations are associated), which implies population specificity
# p=5;pcausal=4
# niter=100
# dg=datgen(m=100,p=p,pcausal=pcausal,h2=0.0005,ngwas_min=3e4,ngwas_max=1e5,niter=niter)
# pnone=function(z,LD) {
#   s
# }
# pall=function(z,LD) {
#   s
# }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ### pairwise assessments, comparing
# # P(d=0)
# # P(d=1)
# # P(d=2)
# # function to generate simulation data
# p=2;pcausal=1
# h2=0.0005
# ngwas=50000
# datgen=function(m,p,pcausal,h2,ngwas_min,ngwas_max,niter) {
#   ngwas=seq(ngwas_min,ngwas_max,length.out=p); ngwas=floor(ngwas)
#   LD=ar1(m)
#   Th=solve(LD)
#   zlist=b0list=row_ixs=col_ixs=list()
#   for(iter in 1:niter) {
#     row_ix=sample(1:m,1); row_ixs[[iter]]=row_ix
#     col_ix=sample(1:p,pcausal); col_ixs[[iter]]=col_ix
#     b0=matrix(0,m,p); b0list[[iter]]=b0
#     b0[row_ix,col_ix]=sqrt(h2)
#     z=matrix(nr=m,nc=p)
#     for(i in 1:p) {
#       bhat=b0[,i]+c(rmvn(1,rep(0,m),Th/ngwas[i]))
#       betahat=LD%*%bhat
#       shat=sqrt(rchisq(m,ngwas[i]-1)/(ngwas[i]-1)/ngwas[i])
#       z[,i]=betahat/shat
#     }
#     zlist[[iter]]=z
#   }
#   list(zlist=zlist,b0list=b0list,row_ixs=row_ixs,col_ixs=col_ixs)
# }
# dg=datgen(m=100,p=p,pcausal=pcausal,h2=h2,ngwas_min=ngwas,ngwas_max=ngwas,niter=1000)
# 
# z=dg$zlist[[1]]
# prior_prob=0.2
# nchain=1000;burnin=100
# chain=list()
# chain[[1]]=rep(0,p)
# for(i in 2:nchain) {
#   gamma0=chain[[i-1]]
#   gamma1=rbinom(1,p,prior_prob)
#   like0=c(); for(j in 1:p) like0[j]=dmvn(z[,j],penmu(z[,j],LD)*gamma0[j],LD)
#   like1=c(); for(j in 1:p) like1[j]=dmvn(z[,j],penmu(z[,j],LD)*gamma1[j],LD)
#   dprior0=dbinom(gamma0,p,prior_prob)
#   dprior1=dbinom(gamma1,p,prior_prob)
#   like0=prod(like0);like1=prod(like1)
#   # dprior0=prod(dprior0);dprior1=prod(dprior1)
#   a=dprior1*like1/(dprior1*like1+dprior0*like0)
#   if(a>runif(1)) {
#     chain[[i]]=gamma1
#   } else {
#     chain[[i]]=gamma0
#   }
#   if(i%%round(floor(nchain*0.01))==0) cat(round(i/nchain*100),'% complete\n',sep='')
# }
# chainmat=do.call(rbind,chain)[-c(1:burnin),]
# dec=rowSums(chainmat)
# (ps=colSums(chainmat)/(niter-burnin)) # P(d_i==1 | data)
# rollprop=function(x,k=100) zoo::rollsum(x,k,align='right')/k
# matplot(apply(chainmat,2,rollprop),cex=2/3)
# cat(dg$col_ixs[[1]])
# 
# sum(rowSums(chainmat)==0)/(niter-burnin) # P(d=0)
# sum(rowSums(chainmat)==1)/(niter-burnin) # P(d=1)
# sum(rowSums(chainmat)==2)/(niter-burnin) # P(d=2)
# 
# 
# 
# 
# p=2;pcausal=1
# h2=0.0001
# ngwas=10000
# dg=datgen(m=100,p=p,pcausal=pcausal,LD=LD,h2=h2,ngwas_min=ngwas,ngwas_max=ngwas,niter=1000)
# z=dg$zlist[[1]]
# library(matrixNormal)
# f=function(lamvec_s,z,LD) {
#   # z=apply(z,2,function(h) h/mad(h))
#   m=nrow(z);p=ncol(z);Ip=diag(p)
#   D=diag(lamvec_s[1:2])
#   regz=z%*%D
#   # regz=apply(z,2,function(h) penmu(h,LD))
#   like=dmatnorm(z,regz,LD,Ip,log=TRUE)/m
#   reg=sum(D)*lamvec_s[3]
#   pen=reg-like # penalty to minimize
#   pen
# }
# optim(par=c(rep(0.1,p),2),f=f,z=z,LD=LD,method='L-BFGS',lower=c(rep(0,p),1),upper=c(rep(1,p),4))
# L=matrix(nr=niter,nc=p)
# for(iter in 1:niter) {
#   # always make first gene the nonzero one
#   ix=dg$col_ixs[[iter]]
#   if(length(ix)>0) {
#     if(pcausal==1 & ix==2) ziter=dg$zlist[[iter]][,c(2,1)] else ziter=dg$zlist[[iter]]
#   } else {
#     ziter=dg$zlist[[iter]]
#   }
#   L[iter,]=optim(par=rep(0.1,p),f=f,z=ziter,LD=LD,method='L-BFGS',lower=rep(0,p),upper=rep(1,p))$par
#   if(iter%%floor(niter*0.1)==0) cat(round(iter/niter*100),'% complete\n',sep='')
# }
# matplot(L,ylim=c(0,1))
# colMeans(L!=0)
# 
# 
# 
# 
# rm(list=ls(all=TRUE))
# ##
# datgen=function(m,p,pcausal,LD,h2,ngwas_min,ngwas_max,niter) {
#   ngwas=seq(ngwas_min,ngwas_max,length.out=p); ngwas=floor(ngwas)
#   Th=solve(LD)
#   zlist=b0list=row_ixs=col_ixs=list()
#   for(iter in 1:niter) {
#     row_ix=sample(1:m,1); row_ixs[[iter]]=row_ix
#     col_ix=sample(1:p,pcausal); col_ixs[[iter]]=col_ix
#     b0=matrix(0,m,p); b0list[[iter]]=b0
#     b0[row_ix,col_ix]=sqrt(h2)
#     z=matrix(nr=m,nc=p)
#     for(i in 1:p) {
#       bhat=b0[,i]+c(rmvn(1,rep(0,m),Th/ngwas[i]))
#       betahat=LD%*%bhat
#       shat=sqrt(rchisq(m,ngwas[i]-1)/(ngwas[i]-1)/ngwas[i])
#       z[,i]=betahat/shat
#     }
#     zlist[[iter]]=z
#   }
#   list(zlist=zlist,b0list=b0list,row_ixs=row_ixs,col_ixs=col_ixs)
# }
# f=function(lamvec,z,LD) {
#   # z=apply(z,2,function(h) h/mad(h))
#   m=nrow(z);p=ncol(z);Ip=diag(p)
#   D=diag(lamvec)
#   regz=z%*%D
#   # regz=apply(z,2,function(h) penmu(h,LD))
#   like=dmatnorm(z,regz,LD,Ip,log=TRUE)/m
#   reg=sum(lamvec)*1
#   pen=reg-like # penalty to minimize
#   pen
# }
# ar1=function(n,rho)rho^toeplitz(0:(n-1))
# ##
# library(mvnfast);library(matrixNormal)
# m=100
# niter=100
# LD=ar1(m,0.9)
# ngwas=c(10000,50000,100000)
# p=2
# pcausal=0:p
# h2=seq(0.001,0.005,length.out=3)
# RES=array(dim=c(length(ngwas),length(pcausal),length(h2)))
# # first pop is associated, second is not
# for(i in 1:length(ngwas)) {
#   for(j in 1:length(pcausal)) {
#     for(k in 1:length(h2)) {
#       dg=datgen(m=100,p=p,pcausal=pcausal[j],LD=LD,h2=h2[k],ngwas_min=ngwas[i],ngwas_max=ngwas[i],niter=niter)
#       res=c()
#       for(iter in 1:niter) {
#         ix=dg$col_ixs[[iter]]
#         if(length(ix)==1) {
#           if(pcausal[j]==1 & ix==2) ziter=dg$zlist[[iter]][,c(2,1)] else ziter=dg$zlist[[iter]]
#         } else if(length(ix) %in% c(0,2)) {
#           ziter=dg$zlist[[iter]]
#         }
#         # L=optim(par=rep(0.1,p),f=f,z=ziter,LD=LD,method='L-BFGS',lower=rep(0,p),upper=rep(1,p))$par
#         L=optim(par=c(rep(0.1,p),2),f=f,z=ziter,LD=LD,method='L-BFGS',lower=c(rep(0,p),2),upper=c(rep(1,p),2.001))$par
#         res[iter]=sum(L[1:2]>0)==pcausal[j]
#       }
#       RES[i,j,k]=sum(res)/niter
#       cat(i,':',j,':',k,'\n',sep='')
#     }
#   }
# }
# pdf=data.frame(
#   res=c(RES),
#   ngwas=rep(rep(ngwas,length(pcausal)),length(h2)),
#   pcausal=rep(rep(pcausal,each=length(ngwas)),length(h2)),
#   h2=rep(h2,each=length(ngwas)*length(pcausal))
# )
# pdf$ngwas=paste0('GWAS n=',pdf$ngwas/1e3,'K')
# pdf$ngwas=factor(pdf$ngwas,levels=paste0('GWAS n=',ngwas/1e3,'K'))
# pdf$h2=paste0(pdf$h2*100,'%')
# pdf$h2=factor(pdf$h2,levels=paste0(h2*100,'%'))
# ggplot(pdf,aes(x=factor(pcausal),y=res,color=factor(h2))) +
#   facet_wrap(~ngwas) +
#   geom_point() +
#   geom_path(aes(group=factor(h2)),linetype='dotted') +
#   lims(y=c(0,1)) +
#   theme_classic() +
#   scale_color_brewer('SNP h2 (%)',palette='Set1') +
#   labs(x='number of causal populations',y='P(correctly inferred pcausal)') +
#   theme(legend.position='bottom') +
#   theme(panel.border=element_rect(colour="black",fill=NA,linewidth=1),
#         strip.background=element_rect(fill="beige"))
# 






## mixture distribution

S=function(x,lambda) sapply(x,function(h) sign(h)*max(c(0,abs(h)-lambda)))
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
  # plot(likes)
  # plot(regs)
  # plot(pens)
  S(z,lams[which.min(pens)])
}
tr=function(x)sum(diag(x))
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
stats=colSums(z^2)
xlyl=800
plot(1,2,xlim=c(-1,1)*xlyl,ylim=c(-1,1)*xlyl,type='n')
mu0=m
sigma20=2*tr(LD%*%LD)
rate0=mu0/sigma20
shape0=mu0^2/sigma20
q0=qgamma(1-0.05,shape=shape0,rate=rate0)
polygon(x=c(-1,1,1,-1)*q0,y=c(-1,-1,1,1)*q0,col='gray90')
abline(h=0,v=0,lty=1)
points(stats[1],stats[2],pch=4)
# P(all not associated)
alt_pars=apply(z,2,function(h) gampars(h,LD,null=FALSE));alt_pars=do.call(rbind,alt_pars)
null_pars=apply(z,2,function(h) gampars(h,LD,null=TRUE));null_pars=do.call(rbind,null_pars)
alt_pars=matrix(c(unlist(alt_pars)),nr=p,nc=2);colnames(alt_pars)=c('shape','rate')
null_pars=matrix(c(unlist(null_pars)),nr=p,nc=2);colnames(null_pars)=c('shape','rate')
p_alts=c()
for(i in 1:nrow(alt_pars)) {
  p1=dgamma(stats[i],shape=alt_pars[i,1],rate=alt_pars[i,2])
  p0=dgamma(stats[i],shape=null_pars[i,1],rate=null_pars[i,2])
  p_alts[i]=p1/(p1+p0)
}
p_alts=abs(p_alts-0.5)*2 # a kind of normalization
pboth=prod(p_alts) # P(mu1!=0,mu2!=0)
peither=sum(p_alts)-prod(p_alts) # P(mu1!=0 or mu2!=0)
pone=peither-pboth # P(mu1!=0 or mu2!=0 and !(mu1!=0 and mu2!=0))
cat('P(only one is associated)=',round(pone,3),'\n',sep='')
# plot the normalization function
s=seq(0,1,0.01)
pos=function(x) ifelse(x<0,0,x)
sf1=function(x) pos(2*(x-0.5))
sf2=function(x,lambda=1) pos(1-exp(-(x-0.5)*lambda))
plot(1,1,type='n',xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='P',ylab='f(P)')
bks=seq(0,1,0.25)
axis(side=1,at=bks,labels=bks)
axis(side=2,at=bks,labels=bks)
abline(v=bks,h=bks,lty=2,col='gray70')
lines(s,sf1(s))
lines(s,sf2(s,lambda=10),col='blue')


marginal_pchisqmix=function(z2,ncp_start=1e-4,max_iter=30,eps=1/length(z2)) {
  d0=dchisq(z2,1)
  d1=dchisq(z2,1,ncp=ncp_start)
  g0=+(d1>d0)
  k=0;error=eps+1
  while(error>eps & k<max_iter) {
    k=k+1
    ncp=mean(z2[g0==1])-1
    d1=dchisq(z2,1,ncp=ncp)
    g=+(d1>d0)
    error=sum(g!=g0)
    g0=g 
  }
  d1/(d1+d0)
}

#############################
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

m=100
LD=ar1(m,0.9)
p=2;pcausal=1
h2=0.0005
ngwas=10000
niter=1000
dg=datgen(m=m,p=p,pcausal=pcausal,LD=LD,h2=h2,ngwas_min=ngwas,ngwas_max=ngwas,niter=niter)
res=c()
xlyl=700
plot(1,2,xlim=c(0,1)*xlyl,ylim=c(0,1)*xlyl,xlab='Population 1 statistic',ylab='Population 2 statistic',type='n',main=expression('h'^2*'=0.05%; GWAS n=50K'))
mu0=m
sigma20=2*tr(LD%*%LD)
rate0=mu0/sigma20
shape0=mu0^2/sigma20
q0=qgamma(1-0.05,shape=shape0,rate=rate0)
polygon(x=c(-1,1,1,-1)*q0,y=c(-1,-1,1,1)*q0,col=NULL,lty=2,border='black')
abline(h=0,v=0,lty=1)
cols=paletteer::paletteer_c("ggthemes::Temperature Diverging", 100)
colbins=seq(0,1,length.out=length(cols))
for(iter in 1:niter) {
  zi=dg$zlist[[iter]]
  stats=colSums(zi^2)
  mp=mixp(zi,LD)
  res[iter]=mp
  # coli=cols[which.min(abs(mp-colbins))]
  # points(stats[1],stats[2],pch=19,cex=1,col=alpha(coli,0.75))
}
# legend('topright',title='P(one population)',cex=2/3,legend=c(1.0,0.5,0.0),pch=c(19,19,19),col=c(tail(cols,1),cols[floor(length(cols)/2)],cols[1]))
# polygon(x=c(-1,1,1,-1)*q0,y=c(-1,-1,1,1)*q0,col=NULL,lty=2,border='black')
linehist=function(x,col,lwd=1,linecol='black',...) {
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
linehist(res,breaks=10,lwd=1,linecol='black',xlab='probabilities of H1',xaxt='n',yaxt='n',ylab=NULL)
axis(side=1,at=c(0,0.5,1),labels=c(0,0.5,1))




rm(list=ls(all=T))
source('simulations/mugent_sel/functions.R')
m=100
LD=ar1(m,0.9)
p=2;pcausal=0:p
h2=0.0001*c(1:5)
ngwas=c(1e4,5e4,1e5)
niter=100
RES=array(dim=c(length(ngwas),length(pcausal),length(h2)))
# first pop is associated, second is not
for(i in 1:length(ngwas)) {
  for(j in 1:length(pcausal)) {
    for(k in 1:length(h2)) {
      dg=datgen(m=100,p=p,pcausal=pcausal[j],LD=LD,h2=h2[k],ngwas_min=ngwas[i],ngwas_max=ngwas[i],niter=niter)
      res=c()
      for(iter in 1:niter) {
        res[iter]=mixp(dg$zlist[[iter]],LD)
      }
      RES[i,j,k]=median(res)
      cat(i,':',j,':',k,'\n',sep='')
    }
  }
}
pdf=data.frame(
  res=c(RES),
  ngwas=rep(rep(ngwas,length(pcausal)),length(h2)),
  pcausal=rep(rep(pcausal,each=length(ngwas)),length(h2)),
  h2=rep(h2,each=length(ngwas)*length(pcausal))
)
pdf$ngwas=paste0('GWAS n=',pdf$ngwas/1e3,'K')
pdf$ngwas=factor(pdf$ngwas,levels=paste0('GWAS n=',ngwas/1e3,'K'))
pdf$h2=paste0(pdf$h2*100,'%')
pdf$h2=factor(pdf$h2,levels=paste0(h2*100,'%'))
ggplot(pdf,aes(x=factor(pcausal),y=res,color=factor(h2))) +
  facet_wrap(~ngwas) +
  geom_point() +
  geom_path(aes(group=factor(h2)),linetype='dotted') +
  lims(y=c(0,1)) +
  theme_classic() +
  scale_color_brewer('SNP h2 (%)',palette='Set1') +
  labs(x='number of causal populations',y='P(distinct pop. inferred)') +
  theme(legend.position='bottom') +
  theme(panel.border=element_rect(colour="black",fill=NA,linewidth=1),
        strip.background=element_rect(fill="beige"))




















