rm(list=ls(all=TRUE))
source('simulations/gent/functions.R')
library(RColorBrewer);library(mvnfast);library(ggplot2)
#########################################################################################
# Type I error
## changing LD density and changing number of SNPs tested
niter=1000
ngwas=50000
Ms=c(50,250,500)
ld_rhos=seq(0,0.9,0.1)
RES=matrix(nr=length(Ms),nc=length(ld_rhos))
for(i in 1:length(Ms)) {
  for(j in 1:length(ld_rhos)) {
    # LD=ar1(Ms[i],ld_rhos[j]);Th=solve(LD) # AR1
    LD=cs(Ms[i],ld_rhos[j]);Th=solve(LD) # compound symmetry
    b=rep(0,Ms[i]) # true joint effect sizes
    bhat=rmvn(niter,b,Th/ngwas) # estimated joint effect sizes
    betahat=t(LD%*%t(bhat))
    s2=rchisq(Ms[i]*niter,ngwas-1)/(ngwas-1)/ngwas
    s2=matrix(s2,nr=niter,nc=Ms[i])
    z=betahat/sqrt(s2)
    ps=apply(z,1,function(h) gent(h,LD)$pval)
    RES[i,j]=sum(ps<0.05)
    cat(i,':',j,'\n',sep='')
  }
}
RES=RES/niter
pdf=data.frame(x=c(RES),m=rep(Ms,length(ld_rhos)),r=rep(ld_rhos,each=length(Ms)))
# plot
ggplot(pdf,aes(r,x,color=factor(m))) +
  geom_hline(yintercept=0.05) +
  geom_hline(yintercept=0.05+c(1,-1)*2*sqrt(0.05*0.95/niter),linetype='dashed') +
  geom_line(lwd=1) +
  scale_color_manual(NULL,labels=paste0(Ms,' SNPs'),values=brewer.pal('Dark2',n=length(Ms))) +
  lims(y=c(0,0.2)) +
  theme_classic() +
  labs(x='density of LD matrix',y='false positive rate',title=paste0('GWAS n=',as.integer(ngwas/1e3),'K')) +
  theme(panel.border=element_rect(color='black',fill=NA,linewidth=1),
        legend.position='inside',
        legend.position.inside=c(0.5,0.7)) +
  scale_x_continuous(breaks=c(0,0.1,0.25,0.5,0.75,0.9),labels=c(0,0.1,0.25,0.5,0.75,0.9),limits=c(0,0.9))
########################################################################################
# Power 
## changing GWAS sample size, polygenicity, and heritability
# parameters
niter=1000
ngwas=10000 # change this manually for 10k and 50k
M=100
mcausals=1:10
h2s=c(0.00005,0.00025,0.0005)
# LD=ar1(M,0.9) # AR1
LD=cs(M,0.9) # compound symmetry
Th=solve(LD)
RES=matrix(nr=max(mcausals),nc=length(h2s))
for(i in 1:max(mcausals)) {
  for(j in 1:length(h2s)) {
    # true joint effect sizes
    b=rep(0,M)
    ix=sample(1:M,mcausals[i]) # put nonzero causal effects in random locations
    b[ix]=sqrt(h2s[j]/mcausals[i])
    bhat=rmvn(niter,b,Th/ngwas) # estimated joint effect sizes
    betahat=t(LD%*%t(bhat))
    s2=rchisq(M*niter,ngwas-1)/(ngwas-1)/ngwas
    s2=matrix(s2,nr=niter,nc=M)
    z=betahat/sqrt(s2)
    ps=apply(z,1,function(h) gent(h,LD)$pval)
    RES[i,j]=sum(ps<0.05)
    cat(i,':',j,'\n',sep='')
  }
}
RES=RES/niter
pdf=data.frame(x=c(RES),mcausal=rep(mcausals,length(h2s)),h2=rep(h2s,each=max(mcausals)))
# plot
ggplot(pdf,aes(factor(mcausal),x,color=factor(h2))) +
  geom_path(aes(group=factor(h2)),lwd=1/2,linetype='dotted') +
  geom_point() +
  scale_color_manual(NULL,labels=paste0('h2=',h2s*100,'%'),values=brewer.pal('Dark2',n=length(h2s))) +
  lims(y=c(0,1)) +
  theme_classic() +
  labs(x='number of causal gene SNPs',y='power of GenT',title=paste0('GWAS n=',as.integer(ngwas/1e3),'K')) +
  theme(panel.border=element_rect(color='black',fill=NA,linewidth=1),
        legend.position='inside',
        legend.position.inside=c(0.5,0.5),
        legend.text=element_text(size=12)
  )
########################################################################################
# Small LD ref size
## changing LD reference and GWAS sample sizes
# parameters
niter=1000
ngwas=10000
nrefs=c(500,1000,5000,10000) # sample size of LD reference panel
Ms=c(50,100,250)
RES=matrix(nr=length(Ms),nc=length(nrefs))
for(i in 1:length(Ms)) {
  # LD=ar1(Ms[i],0.9) # AR1
  LD=cs(Ms[i],0.9) # CS
  Th=solve(LD)
  for(j in 1:length(nrefs)) {
    LDhat=rWishart(niter,nrefs[j],LD)
    b=rep(0,Ms[i]) # no causal effects
    bhat=rmvn(niter,b,Th/ngwas) # estimated joint effect sizes
    betahat=t(LD%*%t(bhat))
    s2=rchisq(Ms[i]*niter,ngwas-1)/(ngwas-1)/ngwas
    s2=matrix(s2,nr=niter,nc=Ms[i])
    z=betahat/sqrt(s2)
    ps=c();for(kk in 1:niter) ps[kk]=gent(z[kk,],cov2cor(LDhat[,,kk]/nrefs[j]))$pval
    RES[i,j]=sum(ps<0.05)
    cat(i,':',j,'\n',sep='')
  }
}
RES=RES/niter
pdf=data.frame(x=c(RES),M=rep(Ms,length(nrefs)),nref=rep(nrefs,each=length(Ms)))
# plot
library(ggplot2)
ggplot(pdf,aes(x=factor(nref),y=x,color=factor(M))) +
  geom_hline(yintercept=0.05+c(1,-1)*2*sqrt(0.05*0.95/niter),linetype='dashed') +
  geom_hline(yintercept=0.05) +
  geom_path(aes(group=factor(M)),lwd=1) +
  scale_color_manual(NULL,labels=paste0(Ms,' SNPs'),values=brewer.pal('Dark2',n=length(Ms))) +
  theme_classic() +
  theme(panel.border=element_rect(color='black',fill=NA,linewidth=1),
        legend.position='inside',
        legend.position.inside=c(0.5,0.7),
        legend.text=element_text(size=10)
  ) +
  guides(colour=guide_legend(nrow=3)) +
  labs(x='size of LD reference panel',y='false positive rate',title=paste0('GWAS n=',as.integer(ngwas/1e3),'K')) +
  scale_x_discrete(breaks=nrefs,labels=paste0(nrefs/1e3,'K')) +
  scale_y_continuous(breaks=c(0,0.05,0.1,0.25),
                     labels=c('0.00','0.05','0.10','0.25'),
                     limits=c(0,0.25))
