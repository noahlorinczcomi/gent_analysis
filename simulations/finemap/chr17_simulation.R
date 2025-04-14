#rm(list=ls(all=TRUE))
library(data.table);library(magrittr);library(tidyr);library(dplyr);library(ggplot2)
library(mvnfast,lib='/home/lorincn/Rpkgs')
library(mvsusieR,lib='/home/lorincn/Rpkgs')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# functions ####
theme_set(
  theme_bw() +
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          strip.background=element_rect(fill="beige"),
          panel.background=element_rect(fill='#EBF9FF',colour='black'))
)
tr=function(x) sum(diag(x))
ar1=function(n,rho) rho^toeplitz(0:(n-1))
CS=function(n,rho) matrix(rho,n,n)+(1-rho)*diag(n)
posDefifyCorrMat=function(ndmatrix,epsilon=1e-4) {
  # ndmatrix: pxp square negative-definite matrix of measurement error correlations/covariances
  # since it is a correlation matrix, abs max off-diagonal value is 1
  deco=eigen(ndmatrix)
  eig0=min(deco$values)
  eigP=max(deco$values)
  if(eig0<epsilon) {
    mu=max(c(epsilon,(eig0+eigP)/2))
    alpha=(mu-epsilon)/(mu-eig0) # optimal choice given by Choi et al: https://doi.org/10.1016/j.jmva.2018.12.002
    Sigmahat=alpha*ndmatrix+(1-alpha)*mu*diag(nrow(ndmatrix)) # the solution
    deco=eigen(Sigmahat)
  } else {
    alpha=1
    Sigmahat=ndmatrix
  }
  return(list(Sigmahat=Sigmahat,alpha=alpha,deco=deco))
}
simdata=function(m,rho,z0=rep(0,m),nref=505,niter=1000) {
  R0=ar1(m,rho)
  Z=mvnfast::rmvn(niter,z0,R0)
  R=rWishart(niter,max(c(m+1,nref)),R0)
  for(i in 1:niter) R[,,i]=cov2cor(R[,,i])
  list(Q=rowSums(Z),Z=Z,R=R,R0=R0)
}
simQs=function(es,vs,R,ngenes=ncol(R),niter=1000) {
  ## R: correlation between GenT statistics
  # normal copula
  rates=es/vs
  shapes=es*rates
  X=mvnfast::rmvn(niter,rep(0,ngenes),R)
  p=pnorm(X)
  Qs=matrix(nr=niter,nc=ngenes)
  for(i in 1:ngenes) Qs[,i]=qgamma(p[,i],shape=shapes[i],rate=rates[i])
  Qs
}
atransform=function(Q,null_mean,null_variance) {
  # function to transform gene-based test statistics to an asymptotically normal distribution
  m=null_mean
  z=Q/m-1
  z*m/sqrt(null_variance)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# simulation for fine-mapping with GenT ####
## need to simulate multiple gene-based association test statistics
## use an actual matrix of gene correlations around AD lead genes
gent.Rho=readRDS('/mnt/isilon/w_gmi/chengflab/Cheng-Noah/reference_data/gent_stat_ld/full_matrices/full_matricesEUR.Rds')
blocks=readRDS('/mnt/isilon/w_gmi/chengflab/Cheng-Noah/manuscripts/aging_and_brain_ROIs/plot_data/GenT_EUR_correlation_blocks.Rds')
## use chromosome 17
chr=17
gent.Rho=gent.Rho[[paste0('chr',chr)]]
gent.Rho=posDefifyCorrMat(gent.Rho)$Sigmahat
gent.Rho=cov2cor(gent.Rho)
blocks=blocks[[paste0('chr',chr)]]
gwasn_vector=c(5e4,1e5,5e5)
rdf=data.frame()
for(i in 1:length(blocks)) {
  spp=paste(blocks[[i]],collapse=',')
  cat('block ',i,' (',spp,')\n',sep='')
  for(j in 1:length(gwasn_vector)) {
    cat(' ',gwasn_vector[j],'\n',sep='')
    ngenes=length(blocks[[i]])
    if(ngenes<2) next
    ncausalgenes=1
    causalix=round(seq(1,ngenes,length.out=ncausalgenes+2))[-c(1,ncausalgenes+2)]
    causalgenes=blocks[[i]][causalix]
    gwasn=gwasn_vector[j]
    h2=0.01 # h2 explained by each causal gene in this locus
    mcausalsnps=2 # number of causal SNPs in each [causal] gene set in this locus
    taus=rep(0,ngenes) # explained heritability divided by number of causal SNPs
    if(length(causalix)>0) taus[causalix]=h2/mcausalsnps
    m=370 # same m for all genes (mean of number tested SNPs in AD analysis)
    LD=ar1(m,0.75) # LD between SNPs for all genes
    R=gent.Rho[blocks[[i]],blocks[[i]]]
    null_mean=m
    null_variance=2*tr(LD%*%LD)
    es=m+taus*gwasn
    vs=2*tr(LD%*%LD)+2*(gwasn*taus)^2+4*gwasn*taus
    niter=1000
    Qs=simQs(es,vs,R,ngenes=ngenes,niter=niter)
    Ys=atransform(Qs,null_mean,null_variance)
    FINEMAPS=matrix(nr=niter,nc=ngenes)
    colnames(Qs)=colnames(Ys)=colnames(FINEMAPS)=rownames(R)
    # pb=txtProgressBar(min=0,max=niter,style=3)
    for(iter in 1:niter) {
      # setTxtProgressBar(pb,iter)
      if(iter %in% seq(1,niter,100)) cat('  ',round(100*iter/niter),'%\n',sep='')
      fit=susie_rss(z=Ys[iter,],R=R,n=gwasn,L=3,max_iter=500)
      FINEMAPS[iter,]=fit$pip
    }
    # close(pb)
    pdf=data.frame(x=c(FINEMAPS),gene=rep(colnames(FINEMAPS),each=niter))
    pdfcut=pdf %>% 
      group_by(gene) %>% 
      summarise(meanpip=mean(x),
                meanpip09=mean(x>0.9)) %>% 
      mutate(iscausal=gene %in% causalgenes)
    rownames(pdfcut)=pdfcut$gene
    pdfcut=pdfcut[blocks[[i]],] %>% mutate(blockpos=1:nrow(pdfcut))
    rdf=rbind(rdf,pdfcut %>% mutate(block=i,gwasn=gwasn))
  }
}
# save rdf








