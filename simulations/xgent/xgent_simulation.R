rm(list=ls(all=TRUE))
library(mvnfast);library(ggplot2);library(dplyr);library(RColorBrewer)
source('simulations/xgent/functions.R')
###########################################################################################
# Type I error and power with changing xQTL and disease h2
m=100
p=3 # number of xQTL types
ngwas=10000 # sample size of disease GWAS
nxqtl=rep(1000,p) # sample size of xQTL GWASs
LD=ar1(m,0.5) # correlation between GWAS SNPs
S=diag(p+1) # sample overlap correlation
K=kronecker(S,LD) # total correlation (row-wise and column-wise) between SNPs
niter=10000 # number of iterations
h2_gwas=round(seq(0,0.0025,length.out=6),4) # h2 of disease GWAS explained by cis gene SNPs
h2_xqtls=round(seq(0,0.025,length.out=6),3) # h2 of xQTL type explained by cis gene SNPs
ncausal=1 # number of causal SNPs
SIG1=SIG2=matrix(nr=length(h2_xqtls),nc=length(h2_gwas))
EG1=EG2=SIG1 # EG1 and EG2 simply track the h2s used
k=0
eg=matrix(nr=prod(dim(SIG1)),nc=2)
for(hx. in 1:length(h2_xqtls)) {
  for(hg. in 1:length(h2_gwas)) {
    k=k+1
    h2s=c(h2_gwas[hg.],rep(h2_xqtls[hx.],p))
    B=matrix(0,m,p+1)
    for(i in 1:(p+1)) B[seq(1,m,length.out=ncausal),i]=sqrt(h2s[i]/ncausal)
    Z0=B*sqrt(matrix(c(ngwas,nxqtl),m,p+1,byrow=T))
    Zall=rmvn(niter,c(Z0),K)
    sig1=sig2=c()
    for(iter in 1:niter) {
      Z=Zall[iter,]
      Z=matrix(c(Z),m,p+1)
      vegas_fit=vegas(Z[,1],LD)
      xvegas_fit=xvegas(Z[,1],LD,Z[,-1])
      sig1[iter]=vegas_fit<0.05
      sig2[iter]=xvegas_fit$pval<0.05
    }
    SIG1[hx.,hg.]=mean(sig1)
    SIG2[hx.,hg.]=mean(sig2)  
    EG1[hx.,hg.]=h2s[1] # gwas
    EG2[hx.,hg.]=h2s[2] # xqtl
    cat(k,'of',prod(dim(SIG1)),'done\n')
  }
}
eg=cbind(c(EG1),c(EG2))
pdf=data.frame(vegas=c(SIG1),xvegas=c(SIG2),h2gwas=eg[,1],h2xqtl=eg[,2])
pdf=data.frame(pow=c(pdf$vegas,pdf$xvegas),type=rep(c('GenT','xGenT'),each=nrow(pdf)),h2gwas=rep(pdf$h2gwas,2),h2xqtl=rep(pdf$h2xqtl,2))
pdf %>% filter(type=='xGenT') %>%
  ggplot(aes(h2xqtl,h2gwas,fill=pow)) +
  geom_tile(color='black') +
  # facet_wrap(~type) +
  theme_classic() +
  theme(legend.position='bottom',panel.background=element_rect(colour="black",linewidth=1)) +
  scale_fill_gradient('Power of xGenT',low='white',high='darkgreen') +
  labs(x='heritability of xQTL',y='heritability of disease',title=paste0('GWAS n=',as.integer(ngwas/1e3),'K')) +
  scale_x_continuous(breaks=h2_xqtls,labels=paste0(h2_xqtls*100,'%')) +
  scale_y_continuous(breaks=h2_gwas,labels=paste0(h2_gwas*100,'%'))
##########################################################################################
# power and number of xQTL phenotypes
m=100
p=10
p_causal=seq(1,10) # number of causal xQTL types
p_noncausal=p-p_causal # number of non-causal xQTL types
ngwas=50000 # sample size of disease GWAS
nxqtl=rep(1000,p) # sample size of xQTL GWASs
LD=ar1(m,0.75) # correlation between GWAS SNPs
S=diag(p+1) # sample overlap correlation
K=kronecker(S,LD) # total correlation (row-wise and column-wise) between SNPs
niter=100 # number of iterations
h2_gwas=0.0005 # h2 of disease GWAS explained by cis gene SNPs
h2_xqtls=0.005 # h2 of xQTL type explained by cis gene SNPs
ncausal_snps=c(1,2,3) # number of causal SNPs
SIG1=SIG2=matrix(nr=length(p_causal),nc=length(ncausal_snps))
EG1=EG2=SIG1 # EG1 and EG2 simply track the h2s used
k=0
eg=matrix(nr=prod(dim(SIG1)),nc=2)
for(i in 1:length(p_causal)) {
  for(j in 1:length(ncausal_snps)) {
    k=k+1
    h2x=rep(0,p)
    h2x[1:p_causal[i]]=sqrt(h2_xqtls/ncausal_snps[j])
    h2s=c(h2_gwas,h2x)
    B=matrix(0,m,p+1)
    for(s. in 1:(p+1)) B[seq(1,m,length.out=ncausal_snps[j]),s.]=sqrt(h2s[s.]/ncausal_snps[j])
    Z0=B*sqrt(matrix(c(ngwas,nxqtl),m,p+1,byrow=T))
    Zall=rmvn(niter,c(Z0),K)
    sig1=sig2=c()
    for(iter in 1:niter) {
      Z=Zall[iter,]
      Z=matrix(c(Z),m,p+1)
      vegas_fit=vegas(Z[,1],LD)
      xvegas_fit=xvegas(Z[,1],LD,Z[,-1])
      sig1[iter]=vegas_fit<0.05
      sig2[iter]=xvegas_fit$pval<0.05
    }
    SIG1[i,j]=mean(sig1)
    SIG2[i,j]=mean(sig2)  
    # EG1[i,j]= # number causal xQTL types
    # EG2[i,j]= # number causal SNPs
    cat(k,'of',prod(dim(SIG1)),'done\n')
  }
}
pdf=data.frame(vegas=c(SIG1),xvegas=c(SIG2),mcausal=rep(ncausal_snps,each=length(p_causal)),pcause=rep(p_causal,length(ncausal_snps)))
pdf %>%
  ggplot(aes(factor(pcause),xvegas,color=factor(mcausal))) +
  geom_path(aes(group=factor(mcausal))) +
  geom_point() +
  theme_classic() +
  theme(legend.position='bottom',
        panel.background=element_rect(colour="black",linewidth=1)) +
  scale_color_manual(NULL,labels=paste0(ncausal_snps,' causal SNPs'),values=brewer.pal('Dark2',n=length(ncausal_snps))) + 
  labs(x='number of causal xQTL phenotypes',y='power of xGenT',title=paste0('GWAS n=',as.integer(ngwas/1e3),'K')) +
  lims(y=c(0,1))
# comparing power of GenT to xGenT
pdflong=data.frame(x=c(pdf$vegas,pdf$xvegas),m=rep(pdf$mcausal,2),p=rep(pdf$pcause,2)) %>% mutate(method=rep(c('GenT','xGenT'),each=nrow(pdf)))
ggplot(pdflong %>% filter(m==1),aes(factor(p),y=x,color=method,fill=method)) +
  geom_ribbon(aes(ymin=x-2*sqrt(x*(1-x)/niter),ymax=x+2*sqrt(x*(1-x)/niter),group=method),alpha=0.15) +
  geom_path(aes(group=method)) +
  geom_point(show.legend=FALSE) +
  scale_color_manual('',values=c('black','#47AD2B')) +
  scale_fill_manual('',values=c('black','#47AD2B')) +
  theme_classic() +
  lims(y=c(-0.01,1.01)) +
  theme(legend.position='none',
        legend.position.inside=c(0.5,0.15),
        panel.background=element_rect(colour="black",linewidth=1)) +
  guides(color=guide_legend(nrow=1)) +
  labs(x='number of xQTL phenotypes',y='power')
################################################################################
# comparing power of xGenT to GenT
ngwas=10000
nxqtl=1000
h2=0.0001
pxqtls=5
h2_xqtl=seq(h2/2,0.01,length.out=10)
ms=c(50,100,250)
niter=1000
RES1=RES2=matrix(nr=length(h2_xqtl),nc=length(ms))
for(i in 1:length(h2_xqtl)) {
  for(j in 1:length(ms)) {
    causalix=floor(ms[j]/2)
    LD=ar1(ms[j],0.9);Th=solve(LD)
    power1=power2=c()
    for(iter in 1:niter) {
      b0=rep(0,ms[j]);b0[causalix]=sqrt(h2)
      l0=matrix(0,ms[j],pxqtls);l0[causalix,]=sqrt(h2_xqtl[i])
      bhat=b0+c(rmvn(1,rep(0,ms[j]),Th/ngwas));betahat=LD%*%bhat
      l0=apply(l0,2,function(h) h+c(rmvn(1,rep(0,ms[j]),Th/nxqtl)))
      zx=betahat*sqrt(ngwas)
      zy=l0*sqrt(nxqtl)
      Zy=matrix(0,ms[j],ms[j]);for(i. in 1:pxqtls) Zy=Zy+zy[,i.]%*%t(zy[,i.])
      L=Zy
      gentres=gent(z=zx,LD=LD)
      xgentres=gent(z=zx,LD=LD,A=L)
      power1[iter]=gentres$pval<0.05
      power2[iter]=xgentres$pval<0.05
    }
    RES1[i,j]=mean(power1)
    RES2[i,j]=mean(power2)
    cat(i,':',j,'\n',sep='')
  }
}

