# rm(list=ls(all=TRUE))
library(data.table);library(dplyr);library(susieR,lib='/home/lorincn/Rpkgs')
###############
pop='EUR'
gene='SYK'
# gwasfp='/home/lorincn/beegfs/lorincn/data/phenotype/T2D/1kgSNPs_Suzuki_EUR_Metal_LDSC-CORR_Neff.v2.txt.gz'
gwasfp='/home/lorincn/beegfs/lorincn/data/phenotype/AD/AD_Bellenguez_2022N_b37.tsv.gz'
# gwasfp='/home/lorincn/beegfs/lorincn/data/phenotype/MDD/MDD_PGC_UKB.txt.gz'
# gwasfp='/home/lorincn/beegfs/lorincn/data/phenotype/PD/Kim_PD.txt.gz'
# gwasfp='/home/lorincn/beegfs/lorincn/data/phenotype/ALS/ALS_RheenenW.EUR_2021.tsv.gz'
# gwasfp='/home/lorincn/beegfs/lorincn/data/phenotype/SCZ/SCZ_PGC_EUR.tsv.gz' # SCZ
kb_window=500 # window around the gene center in which SuSiE will be performed
n=450000 # GWAS sample size (450K for AD; 807K for MDD; 611K for PD; 139K for ALS)
## T2D N's
# EUR: 752000
# AFR: 129000
# EAS: 247000
# SAS: 42000
# HIS: 73000
###############
gene.=gene
genepositions=fread('/home/lorincn/isilon/Cheng-Noah/manuscripts/druggable_genes/data/genepositions_hg19.txt')
genei=genepositions %>% filter(symbol==gene.) %>% head(.,1)
# starti=genei$start;endi=genei$end
gwas=fread(gwasfp) %>% as.data.frame() %>% filter(chr==genei$chr)
gwas=gwas %>% filter(abs(position-genei$mid)<(kb_window*1e3))
if(pop=='HIS') ldpop='AMR' else ldpop=pop
ldref=paste0('/home/lorincn/beegfs/lorincn/data/reference_panels/1kg.v3/',ldpop)
bim=fread(paste0(ldref,'.bim')) %>% filter(V1==genei$chr)
colnames(bim)[1:2]=c('chr','rsid')
gwas=gwas %>% filter(rsid %in% bim$rsid) %>% arrange(position)
dupes=gwas %>% group_by(rsid) %>% summarise(n=n()) %>% filter(n>1) %>% pull(rsid)
gwas=gwas %>% filter(!(rsid %in% dupes))
m=nrow(gwas)
gwas=left_join(gwas,bim %>% select(rsid,ea=V5),by='rsid') %>% mutate(beta=ifelse(effect_allele==ea,beta,-beta)) %>% arrange(position)
randir=paste(sample(c(letters,1:10),replace=T),collapse='')
setwd('/home/lorincn/beegfs/lorincn/working_space')
dir.create(randir)
setwd(randir)
write.table(gwas$rsid,'temp.txt',row.names=FALSE,quote=FALSE,col.names=FALSE)
system(paste0('plink --bfile ',ldref,' --r square --extract temp.txt --out temp'),intern=T)
ld=fread('temp.ld') %>% as.matrix()
ld=0.99*ld+0.01*diag(m)
res=susie_rss(z=gwas$beta/gwas$se,R=ld,n=n)
max(res$pip)
cat(gwas$rsid[res$pip>0.9],'\n')
cat('Pop: ',pop,'\nMax PIP: ',max(res$pip),'\nMax PIP SNPs: ',gwas$rsid[which.max(res$pip)],'\nPIP lead GWAS SNP: ',res$pip[which.min(gwas$pval)],'\n',sep='')
######################################
# If you want to save these data to plot manually with GGplot
setwd('~/isilon/Cheng-Noah')
out=list(susie_res=res,gwas=gwas,ld=ld)
saveRDS(out,'temp.Rds')
######################################
# LOCALLY
setwd('/Volumes/chengflab/Cheng-Noah')
data=readRDS('temp.Rds')
res=data$susie_res
gwas=data$gwas
ld=data$ld
library(dplyr);library(ggplot2);library(ggrepel)
ix=which.min(gwas$pval)
lead=gwas$rsid[ix]
lddf=data.frame(r2=c(ld[ix,])^2,rsid=gwas$rsid)
gwas=left_join(gwas,lddf,by='rsid')
gwas$pip=res$pip
ggplot(gwas,aes(position/1e6,-log10(pval),color=r2)) +
  annotate('rect',xmin=93564069/1e6,xmax=93660831/1e6,ymin=-Inf,ymax=Inf,fill='gold') +
  # annotate('text',x=(93564069+93660831)/2/1e6-0.07,y=7.5,label='SYK gene',angle=90) +
  # annotate('text',x=93.2,y=7.5,label='P<5E-8') +
  geom_point() +
  scale_color_steps2(expression(r^2),
                     low='blue',mid='green',high='red',
                     breaks=(0:4)/4,midpoint=0.5) +
  geom_point(aes(position/1e6,-log10(pval),color=r2),
             shape=17,size=4,
             data=gwas %>% 
               filter(pip>0.9,pval<0.01)) +
  theme_bw() +
  lims(y=c(0,-log10(5e-9))) +
  geom_hline(yintercept=-log10(5e-8),linetype='dashed') +
  theme(panel.border=element_rect(color='black',fill=NA,linewidth=1)) +
  labs(x='Chromosome 9 (Mb)',y=expression('-log'[10]*'(AD GWAS P-value)')) +
  geom_text_repel(aes(position/1e6,-log10(pval),label=rsid),
                  data=gwas %>% 
                    filter(pip==1,pval<0.01),
                  color='black',
                  # nudge_x=.15,
                  # xlim=c(93.3,93.9),
                  # ylim=c(3,8),
                  box.padding=0.5,
                  nudge_y=1,
                  force=100,
                  segment.curvature=0.2,
                  segment.ncp=1,
                  segment.angle=-90) +
  theme(legend.position='none',
        legend.position.inside=c(0.85,0.7),
        legend.box.background=element_rect(color="black"))








################################################################################
# RIPK2 plot
################################################################################
## Ran part of figure4plots.R for this
library(locuscomparer);library(ggpubr)
merged=readRDS('/Volumes/chengflab/Cheng-Noah/manuscripts/druggable_genes/working_space/merged_xqtl_ad_gwas_data_for_novel_gene.Rds')
ld=readRDS('/Volumes/chengflab/Cheng-Noah/manuscripts/druggable_genes/working_space/ld_xqtl_ad_gwas_data_for_novel_gene.Rds')
library(magrittr)
longdf=merged %$% data.frame(
  rsid=rep(rsid,4),
  position=rep(position,4),
  z=c(ad_beta/ad_se,
      brain_cerebellum_beta/brain_cerebellum_se,
      brain_cortex_beta/brain_cortex_se,
      pqtl_beta/as.numeric(pqtl_se)),
  trait=rep(c('AD','Cerebellum eQTL','Cortex eQTL','Cortex pQTL'),each=nrow(merged))
) %>%
  mutate(pval=pchisq(z^2,1,lower.tail=FALSE))
# assign LD with lead SNPs
ts=unique(longdf$trait)
rdf=data.frame()
for(i in 1:length(ts)) {
  dati=longdf %>% filter(trait==ts[i])
  lead_ix=which.max(abs(dati$z))
  lddf=data.frame(rsid=merged$rsid,r2=c(ld[,lead_ix])^2)
  dati=left_join(dati,lddf,by='rsid')
  rdf=rbind(rdf,dati)
}
longdf=rdf;rm(rdf)
ripk2_finemaps=c('rs7815039','rs12156013','rs13254219','rs218935','rs218938','rs425721')
longdf %>%
  filter((position/1e6)>90.25,(position/1e6)<91.25) %>%
  ggplot(aes(position/1e6,-log10(pval),color=r2)) +
  annotate('rect',xmin=90769975/1e6,xmax=90803291/1e6,ymin=-Inf,ymax=Inf,fill='gold') +
  # annotate('text',x=93.2,y=7.5,label='P<5E-8') +
  geom_point() +
  facet_wrap(~trait,nrow=2) +
  scale_color_steps2(expression(r^2),
                     low='blue',mid='green',high='red',
                     breaks=(0:4)/4,midpoint=0.5) +
  theme_bw() +
  # lims(y=c(0,-log10(5e-9))) +
  geom_hline(yintercept=-log10(5e-8),linetype='dashed') +
  theme(panel.border=element_rect(color='black',fill=NA,linewidth=1)) +
  labs(x='Chromosome 8 (Mb)',y=expression('-log'[10]*'(P-value)')) +
  geom_text_repel(aes(position/1e6,-log10(pval),label=rsid),
                  data=longdf %>%
                    filter(rsid %in% ripk2_finemaps,trait=='AD'),
                  color='black',
                  # nudge_x=.15,
                  # xlim=c(93.3,93.9),
                  # ylim=c(3,8),
                  box.padding=0.5,
                  nudge_y=1,
                  force=50,
                  segment.curvature=0.2,
                  segment.ncp=1,
                  segment.angle=-90) +
  theme(legend.position='none',
        legend.position.inside=c(0.15,0.7),
        legend.box.background = element_rect(colour = "black"),
        strip.background =element_rect(fill="beige")) +
  scale_x_continuous(breaks=c(90.25,90.75,91.20),
                     labels=c(90.25,90.75,91.20))


################################################################################
# NTRK1 plot
## Ran part of figure4plots.R for this
merged=readRDS('/Volumes/chengflab/Cheng-Noah/manuscripts/druggable_genes/working_space/merged_xqtl_ad_gwas_data_for_novel_gene.Rds')
ld=readRDS('/Volumes/chengflab/Cheng-Noah/manuscripts/druggable_genes/working_space/ld_xqtl_ad_gwas_data_for_novel_gene.Rds')
library(magrittr)
ix=which((merged$position/1e6)>156.5 & (merged$position/1e6)<157.5)
merged=merged[ix,]
ld=ld[ix,ix]
longdf=merged %$%
  data.frame(
  rsid=rep(rsid,4),
  position=rep(position,4),
  z=c(ad_beta/ad_se,
      brain_spinal_cord_beta/brain_spinal_cord_se,
      brain_cortex_beta/brain_cortex_se,
      pqtl_beta/as.numeric(pqtl_se)),
  trait=rep(c('AD','Spinal cord eQTL','Cortex eQTL','Cortex pQTL'),each=nrow(merged))
) %>%
  mutate(pval=pchisq(z^2,1,lower.tail=FALSE))
# assign LD with lead SNPs
ts=unique(longdf$trait)
rdf=data.frame()
for(i in 1:length(ts)) {
  dati=longdf %>% filter(trait==ts[i])
  lead_ix=which.max(abs(dati$z))
  lddf=data.frame(rsid=merged$rsid,r2=c(ld[,lead_ix])^2)
  dati=left_join(dati,lddf,by='rsid')
  rdf=rbind(rdf,dati)
}
longdf=rdf;rm(rdf)
ntrk1_finemaps=c('rs61816272','rs3765787','rs1758186','rs928392','rs10157555')
longdf %>%
  na.omit() %>%
  ggplot(aes(position/1e6,-log10(pval),color=r2)) +
  annotate('rect',xmin=156785432/1e6,xmax=156851642/1e6,ymin=-Inf,ymax=Inf,fill='gold') +
  geom_point() +
  facet_wrap(~trait,nrow=4,scales='free_x') +
  scale_color_steps2(expression(r^2),
                     low='blue',mid='green',high='red',
                     breaks=(0:4)/4,midpoint=0.5) +
  theme_bw() +
  geom_hline(yintercept=-log10(5e-8),linetype='dashed') +
  theme(panel.border=element_rect(color='black',fill=NA,linewidth=1)) +
  labs(x='Chromosome 8 (Mb)',y=expression('-log'[10]*'(P-value)')) +
  geom_text_repel(aes(position/1e6,-log10(pval),label=rsid),
                  data=longdf %>%
                    filter(rsid %in% ntrk1_finemaps,trait=='AD'),
                  color='black',
                  # nudge_x=.15,
                  # xlim=c(93.3,93.9),
                  # ylim=c(3,8),
                  box.padding=0.5,
                  nudge_y=1,
                  force=50,
                  segment.curvature=0.2,
                  segment.ncp=1,
                  segment.angle=-90) +
  theme(legend.position='none',
        legend.position.inside=c(0.15,0.7),
        legend.box.background = element_rect(colour = "black"),
        strip.background =element_rect(fill="beige"))


################################################################################
### PPI with NTRK1
library(dplyr)
setwd('/Volumes/chengflab/Cheng-Noah/zhendong_network_project')
ppi=readRDS('PPI_matrix.Rds')
diag(ppi)=1
ix=which(rownames(ppi)=='NTRK1')
ppic=ppi[,ix]
ppic=ppic[ppic>0]
ixn=names(ppic)
ix=which(rownames(ppi) %in% ixn)
ppi=ppi[ix,ix]
rm(ppic)
diag(ppi)=0
ix=which(rownames(ppi)=='NTRK1')
ppi=ppi[,ix]
ppi=as.matrix(ppi)
ppi=ppi[-which(rownames(ppi)=='NTRK1'),,drop=FALSE]
colnames(ppi)='NTRK1'
library(GGally);library(network)
adgent=readRDS('/Volumes/chengflab/Cheng-Noah/GenT_results/AD_EUR.Rds')
adgent=adgent %>% 
  mutate(fdr=p.adjust(vegas,'BH')) %>%
  filter(gene %in% rownames(ppi))
net=network(ppi,directed=FALSE)
network.vertex.names(net)=rownames(ppi)
sig_adgenes=adgent %>% filter(fdr<0.05) %>% pull(gene)
# gene_class=ifelse(rownames(ppi) %in% sig_adgenes,'a','b')
# gene_class[which(rownames(ppi)=='NTRK1')]='c'
net %v% "class"=network.vertex.names(net) %in% sig_adgenes
p=ggnet2(net,
         mode='circrand',
         color="class",
         size="degree",
         label=NULL,  # We'll handle labels separately with geom_label_repel
         edge.size=1/8,
         edge.alpha=1/2,
         edge.color='#FAC42A')
p +
  geom_point(aes(x=x,y=y),
             color='red',
             data=p$data %>% 
               filter(label %in% sig_adgenes)) +
  scale_color_manual('',
                     labels=c('AD-associated', 'not AD-associated'),
                     values=c('blue','red')) +
  theme(legend.position='bottom',legend.text=element_text(size=12)) +
  guides(colour=guide_legend(nrow=1,override.aes=list(size=5)),size='none')


setwd('/home/lorincn/beegfs/lorincn/data/phenotype/AD')
ad=fread('AD_Bellenguez_2022N_b37.tsv.gz') %>% filter(chr==1)
setwd('/home/lorincn/beegfs/lorincn/data/gene/eQTL/GTExV8/brain_cortex/')
eqtl=fread('Brain_Cortex.v8.EUR.allpairs.chr1.txt.gz')
gp=fread('/home/lorincn/isilon/Cheng-Noah/manuscripts/druggable_genes/data/genepositions_hg19.txt')
ensembl_ids=sapply(eqtl$phenotype_id,function(h) unlist(strsplit(h,'[.]'))[1])
ix=which(ensembl_ids %in% gp$ensembl_ids[gp$symbol=='NTRK1'])
eqtl=as.data.frame(eqtl)[ix,]
eqtl$symbol='NTRK1'
merged=inner_join(
  ad %>% select(rsid,position,ad_beta=beta,ad_se=se,ea1=effect_allele),
  eqtl %>% 
    filter(symbol=='NTRK1') %>%
    select(rsid=rsID,eqtl_beta=slope,eqtl_se=slope_se,ea2=ALT),
  by='rsid'
) %>%
  mutate(eqtl_beta=ifelse(ea1==ea2,eqtl_beta,-eqtl_beta))
bim=fread('/home/lorincn/beegfs/lorincn/data/reference_panels/1kg.v3/EUR.bim') %>% 
  filter(V1==1)
merged=inner_join(merged,
                  bim %>% select(rsid=V2,a1=V5),by='rsid')
mergedout=merged %>%
  mutate(zx=ad_beta/ad_se,
         zy=eqtl_beta/eqtl_se) %>%
  mutate(zx=ifelse(ea1==a1,zx,-zx),
         z2=ifelse(ea2==a1,zy,-zy)) %>%
  arrange(position)
# get LD
source('/home/lorincn/isilon/Cheng-Noah/software/corefunctions/functions.R')
ld=get_ld(mergedout$rsid,'/home/lorincn/beegfs/lorincn/data/reference_panels/1kg.v3/EUR')
ld=0.91*ld+(0.09*diag(nrow(ld)))
# save and load locally
out=list(merged=mergedout,ld=ld)
saveRDS(out,'/home/lorincn/isilon/Cheng-Noah/temp.Rds')
#### locally
############ you can load from my Documents
load(file='/Users/lorincn/Documents/working_space/NTRK1_mr.Rdata')
library(dplyr)
source('/Volumes/chengflab/Cheng-Noah/software/corefunctions/functions.R')
data=readRDS('/Volumes/chengflab/Cheng-Noah/temp.Rds')
merged=data$merged
ld=data$ld

d=eigen(ld)$values
meff=sum(+(d>=1)+d*(d<1))
q0m=sqrt(qchisq(1-0.05/meff,1))
q0=1.96
plot(merged$eqtl_beta,merged$ad_beta,type='n',
     ylim=c(-0.1,0.1),xlim=c(-1,1.5))
issig=merged$zx^2>q0m & merged$zy^2>q0m
for(i in 1:nrow(merged)) {
  # bx
  y=rep(merged$ad_beta[i],2)
  x=merged$eqtl_beta[i]+c(1,-1)*q0*merged$eqtl_se[i]
  if(issig[i]) {
    lines(x,y,col='blue')
  } else {
    lines(x,y,col='#0AA5FF35')
  }
  # by
  x=rep(merged$eqtl_beta[i],2)
  y=merged$ad_beta[i]+c(1,-1)*q0*merged$ad_se[i]
  if(issig[i]) {
    lines(x,y,col='blue')
  } else {
    lines(x,y,col='#0AA5FF35')
  }
}
points(merged$eqtl_beta[issig],merged$ad_beta[issig],col='blue',pch=19)
# do again
for(i in 1:nrow(merged)) {
  if(!issig[i]) next
  # bx
  y=rep(merged$ad_beta[i],2)
  x=merged$eqtl_beta[i]+c(1,-1)*q0*merged$eqtl_se[i]
  lines(x,y,col='blue')
  # by
  x=rep(merged$eqtl_beta[i],2)
  y=merged$ad_beta[i]+c(1,-1)*q0*merged$ad_se[i]
  lines(x,y,col='blue')
}
abline(h=0,v=0)
ellpar1=ell(c(0,0),cov(cbind(merged$eqtl_beta,merged$ad_beta)),level=1-0.05)
ellpar2=ell(c(0,0),cov(cbind(merged$eqtl_beta,merged$ad_beta)),level=1-0.05/meff)
lines(ellpar1[,1],ellpar1[,2],col='black')
lines(ellpar2[,1],ellpar2[,2],col='black',lty=2)
# MR
library(MRBEE)
by=merged$ad_beta
byse=merged$ad_se
bx=merged$eqtl_beta
bxse=merged$eqtl_se
ivs=clumpfun(merged$rsid,
             merged$position,
             pchisq(bx^2/bxse^2,1,lower.tail=FALSE),
             ld,
             clump_p=5e-5,clump_kb=500,clump_r2=0.1)
ix=which(merged$rsid %in% names(ivs))
bx=bx[ix]
by=by[ix]
bxse=bxse[ix]
byse=byse[ix]
res=MRBEE.IMRP.UV(by=by,bx=bx,byse=byse,bxse=bxse,Rxy=matrix(0,1,2))
pchisq(res$theta^2/res$vartheta,1,lower.tail=FALSE)
abline(a=0,b=res$theta,col='black')














