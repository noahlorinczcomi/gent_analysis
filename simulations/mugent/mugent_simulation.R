rm(list=ls(all=TRUE))
library(RColorBrewer);library(ggplot2);library(mvnfast);library(ggplot2);library(corrplot)
source('simulations/mugent/functions.R')
#########################################################################################
# Power
## changing genetic correlation and heritability exlpained
# source('/home/lorincn/isilon/Cheng-Noah/software/corefunctions/functions.R')
ar1=function(n,rho=0.9) rho^toeplitz(0:(n-1))
niter=10000
ngwas=10000
M=100
mcausal=2
LD=ar1(M,0.9); # assume same LD structure for now (doesnt matter for MuGenT anyway)
LDlist=list(); for(i in 1:2) LDlist[[i]]=LD
Th=solve(LD)
p=2
rgs=c(-0.9,-0.5,0,0.5,0.9)
mshared_causal_snps=c(0,1,2)
h2=0.0005
D=diag(sqrt(h2/mcausal),2)
RES=matrix(nr=length(rgs),nc=length(mshared_causal_snps))
for(i in 1:length(rgs)) {
  Rg=matrix(rgs[i],2,2);diag(Rg)=1
  Gc=D%*%Rg%*%D
  for(j in 1:length(mshared_causal_snps)) {
    if(mshared_causal_snps[j]==0) { # they dont share any of 2 causal SNPs
      ix1=sample(1:M,mcausal,replace=FALSE)
      ix2=sample(c(1:M)[-ix1],mcausal,replace=FALSE)
    } else if(mshared_causal_snps[j]==1) { # they share 1 of 2 causal SNPs
      ix1=sample(1:M,mcausal,replace=FALSE)
      ix2=c(ix1[1],sample(c(1:M)[-ix1[1]],1,replace=FALSE))
    } else { # they share 2 of 2 causal SNPs
      ix1=sample(1:M,mcausal,replace=FALSE)
      ix2=ix1
    }
    ps=c()
    for(iter in 1:niter) {
      # needs to be genetic correlation between shared causal SNPs only - not non-shared
      # rows are causal SNPs, columns are populations
      if(mshared_causal_snps[j]==0) { # they dont share any of 2 causal SNPs
        b0=rmvn(mcausal,c(0,0),D^2)
      } else if(mshared_causal_snps[j]==1) { # they share 1 of 2 causal SNPs
        b0=rbind(
          rmvn(1,c(0,0),Gc), # must be first because the shared one is the first index
          rmvn(1,c(0,0),D^2)
          )
      } else { # they share 2 of 2 causal SNPs
        b0=rmvn(mcausal,c(0,0),Gc)
      }
      b01=b02=rep(0,M)
      b01[ix1]=b0[,1] # true joint effects in population 1
      b02[ix2]=b0[,2] # true joint effects in population 2
      adj1=h2/sum(b01^2);b01=sqrt(adj1)*b01
      adj2=h2/sum(b02^2);b02=sqrt(adj2)*b02
      bhat1=c(rmvn(1,b01,Th/ngwas)) # estimated joint effects in population 1
      bhat2=c(rmvn(1,b02,Th/ngwas)) # estimated joint effects in population 2
      betahat1=LD%*%bhat1 # estimated marginal effects in population 1
      betahat2=LD%*%bhat2 # estimated marginal effects in population 2
      s21=rchisq(M,ngwas-1)/(ngwas-1)/ngwas
      s22=rchisq(M,ngwas-1)/(ngwas-1)/ngwas
      z1=betahat1/sqrt(s21)
      z2=betahat2/sqrt(s22)
      ps[iter]=mugent(cbind(z1,z2),LDlist)$p
    }
    RES[i,j]=sum(ps<0.05)
    cat(i,':',j,'\n',sep='')  
  }
}
RES=RES/niter
pdf=data.frame(x=c(RES),rg=rep(rgs,length(mshared_causal_snps)),mshared_causal_snps=rep(mshared_causal_snps,each=length(rgs)))
# plot
ggplot(pdf,aes(rg,x,color=factor(mshared_causal_snps))) +
  geom_point() +
  geom_line(lwd=1) +
  scale_color_manual(NULL,labels=paste0(mshared_causal_snps,' shared causal SNP',c('s','','s')),values=brewer.pal('Dark2',n=length(mshared_causal_snps))) +
  theme_classic() +
  labs(x='genetic correlation between populations',y='power',title=paste0('GWAS n=',as.integer(ngwas/1e3),'K')) +
  lims(y=c(0,1)) +
  theme(
    # legend.position='inside',
    # legend.position.inside=c(0.6,0.4),
    # legend.text=element_text(size=10),
    legend.position='none',
    panel.border=element_rect(color='black',fill=NA,linewidth=1)
  ) +
  scale_x_continuous(breaks=rgs,labels=rgs,limits=c(-1,1))
#########################################################################################
# Type I error
## changing number of populations and SNPs tested
niter=100
ngwas=50000
Ms=c(50,100,250)
ps=2:6
RES=matrix(nr=length(Ms),nc=length(ps))
for(i in 1:length(Ms)) {
  LD=ar1(Ms[i],0.9);Th=solve(LD)
  for(j in 1:length(ps)) {
    LDlist=lapply(1:ps[j],function(h) LD)
    b=rep(0,Ms[i]*ps[j])
    iK=kronecker(diag(ps[j]),Th)
    bhat=rmvn(niter,b,iK/ngwas)
    pvals=c()
    for(iter in 1:niter) {
      bi=matrix(c(bhat[iter,]),Ms[i],ps[j],byrow=FALSE)
      betahat=LD%*%bi
      s2=matrix(rchisq(Ms[i]*ps[j],ngwas-1)/(ngwas-1)/ngwas,Ms[i],ps[j])
      z=betahat/sqrt(s2)
      pvals[iter]=mugent(z,LDlist)$p
    }
    RES[i,j]=sum(pvals<0.05)
    cat(i,':',j,'\n',sep='')
  }
}
RES=RES/niter
pdf=data.frame(x=c(RES),M=rep(Ms,length(ps)),p=rep(ps,each=length(Ms)))
ggplot(pdf,aes(factor(p),x,color=factor(M))) +
  geom_hline(yintercept=0.05) +
  geom_hline(yintercept=0.05+c(1,-1)*2*sqrt(0.05*0.95/niter),linetype='dashed') +
  geom_line(lwd=1) +
  geom_path(aes(group=factor(M)),lwd=1/3,linetype='dotted') +
  geom_point() +
  scale_color_manual(NULL,labels=paste0(Ms,' SNPs'),values=brewer.pal('Dark2',n=length(Ms))) +
  theme_classic() +
  labs(x='number of populations',y='false positive rate',title=paste0('GWAS n=',as.integer(ngwas/1e3),'K')) +
  theme(
    legend.position='inside',
    legend.position.inside=c(0.5,0.2),
    legend.text=element_text(size=10),
    # legend.position='none',
    panel.border=element_rect(color='black',fill=NA,linewidth=1)
  ) +
  lims(y=c(0,0.2))
#########################################################################################
# Type I error
## when LD structures between populations are different
### 3 populations
p=3
bd=function(A,B) { # block diagonal
  n1=nrow(A);p1=ncol(A)
  n2=nrow(B);p2=ncol(B)
  zero1=matrix(0,nr=n2,nc=p1)
  zero2=matrix(0,nr=n1,nc=p2)
  cbind(rbind(A,zero1),rbind(zero2,B))
}
# changing size of LD reference panels
niter=1000
ngwas=50000
nrefs=c(500,1000,5000,10000)
Ms=c(50,100,250)
RES=matrix(nr=length(Ms),nc=length(nrefs))
for(i in 1:length(Ms)) {
  # define LD
  LD1=ar1(Ms[i],0.9)
  LD2=diag(Ms[i])
  LD3=matrix(0.2,Ms[i],Ms[i]);diag(LD3)=1
  # construct meta LD matrix
  metaLD=bd(LD1,LD2)
  metaLD=bd(metaLD,LD3)
  metaLDi=bd(solve(LD1),solve(LD2))
  metaLDi=bd(metaLDi,solve(LD3))
  ldlist=list(LD1=LD1,LD2=LD2,LD3=LD3)
  # iterate over varying sizes of LD reference panels
  for(j in 1:length(nrefs)) {
    bhat=rmvn(niter,rep(0,Ms[i]*p),metaLDi/ngwas)
    s2=rchisq(niter*Ms[i]*p,ngwas-1)/(ngwas-1)/ngwas
    s2=array(s2,dim=c(Ms[i],p,niter))
    betahat=t(metaLD%*%t(bhat))
    LDhat1=rWishart(niter,nrefs[j],LD1)
    LDhat2=rWishart(niter,nrefs[j],LD2)
    LDhat3=rWishart(niter,nrefs[j],LD3)
    pvals=c()
    for(iter in 1:niter) {
      LDhatlist=list(LD1=LDhat1[,,iter],LD2=LDhat2[,,iter],LD3=LDhat3[,,iter])
      LDhatlist=lapply(LDhatlist,cov2cor)
      betahati=c(betahat[iter,])
      betahati=matrix(betahati,nr=Ms[i],nc=p,byrow=FALSE)
      s2i=s2[,,iter]
      zi=betahati/sqrt(s2i)
      pvals[iter]=mugent(zi,ldlist)$p
    }
    RES[i,j]=sum(pvals<0.05)
    cat(i,':',j,'\n',sep='')
  }
}
RES=RES/niter
pdf=data.frame(m=rep(Ms,length(nrefs)),nref=rep(nrefs,each=length(Ms)),x=c(RES))
ggplot(pdf,aes(x=factor(nref),y=x,color=factor(m))) +
  geom_hline(yintercept=0.05) +
  geom_hline(yintercept=0.05+c(1,-1)*2*sqrt(0.05*0.95/niter),linetype='dashed') +
  geom_line(lwd=1) +
  geom_path(aes(group=factor(m)),lwd=1) +
  theme_classic() +
  scale_color_manual(NULL,labels=paste0(Ms,' SNPs'),values=brewer.pal('Dark2',n=length(Ms))) +
  labs(x='size of LD reference panel',y='false positive rate',title=paste0('GWAS n=',as.integer(ngwas/1e3),'K')) +
  theme(
    legend.position='inside',
    legend.position.inside=c(0.5,0.2),
    legend.text=element_text(size=10),
    panel.border=element_rect(color='black',fill=NA,linewidth=1)
  ) +
  lims(y=c(0,0.2)) +
  scale_x_discrete(breaks=nrefs,labels=paste0(nrefs/1e3,'K'))



