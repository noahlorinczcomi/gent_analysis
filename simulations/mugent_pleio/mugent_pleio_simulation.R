rm(list=ls(all=TRUE))
library(dplyr);library(ggplot2);library(mvnfast)
source('simulations/mugent_pleio/functions.R')
##################################################################################
##################################################################################
## example of usage (generating data)
m=100 # number of SNPs
p=5 # number of traits
ngwass=c(1e4,5e4,1e5)
LD=0.9^toeplitz(0:(m-1)) # simulated LD matrix
ncausal=1 # number of simulated causal SNPs (you can change)
niter=1000
h2s=c(0,seq(0.0001,0.005,length.out=11))
RES=matrix(nr=length(h2s),nc=length(ngwass))
for(i in 1:length(h2s)) {
  for(j in 1:length(ngwass)) {
    res=c()
    for(iter in 1:niter) {
      B0=matrix(0,m,p) # matrix of true SNP causal effect sizes to be filled in
      if(ncausal>0) {
        causalix=seq(1,m,length.out=ncausal)
        if(ncausal==1) causalix=floor(m/2)
        B0[causalix,]=sqrt(h2s[i]/ncausal) # SNP is associated with all traits
      }
      Z0=B0*sqrt(ngwass[j]) # 'true' Z-statistics
      Z=t(rmvn(p,c(Z0[,1]),LD)) # simulated matrix of observed Z-statistics
      res[iter]=any(f(Z,LD)$result)
    }
    RES[i,j]=sum(res)/niter
    cat(i,':',j,'\n',sep='')
  }
}
pdf=data.frame(h2=rep(h2s,length(ngwass)),gwasn=rep(ngwass,each=length(h2s)),power=c(RES))
pdf$h2=paste0(round(pdf$h2*100,2),'%')
pdf$h2=factor(pdf$h2,levels=paste0(round(h2s*100,2),'%'))
pdf$gwasn=paste0(pdf$gwasn/1e3,'K')
pdf$gwasn=factor(pdf$gwasn,levels=paste0(ngwass/1e3,'K'))
ggplot(pdf,aes(factor(h2),power,color=factor(gwasn))) +
  annotate('rect',xmin=2/3,xmax=4/3,ymin=0,ymax=1,fill='gray90') +
  annotate('text',x='0%',y=1/2,label='type I error',angle=90) +
  geom_path(aes(group=factor(gwasn)),linetype='dotted',position=position_dodge(0.5)) +
  geom_point(size=1.5,position=position_dodge(0.5)) +
  theme_classic() +
  theme_classic() +
  lims(y=c(0,1)) +
  scale_color_brewer('GWAS sample size',palette='Set1') +
  labs(x='heritability explained by pleiotropic SNP (%)',y='power to detect pleiotropy') +
  theme(panel.border=element_rect(colour="black",fill=NA,linewidth=1)) +
  theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1)) +
  # theme(legend.position=c(0.75,0.25))
  theme(legend.position='none')









