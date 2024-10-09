rm(list=ls(all=TRUE))
library(data.table);library(dplyr)
source('/home/lorincn/isilon/Cheng-Noah/software/corefunctions/functions.R')
#####
#####
tr=function(x)sum(diag(x))
xvegas=function(z,R,Z_xqtls) {
  z=as.matrix(z);Z_xqtls=as.matrix(Z_xqtls)
  m=nrow(R);p=ncol(Z_xqtls)
  L=matrix(0,nrow=nrow(Z_xqtls),ncol=nrow(Z_xqtls))
  for(i in 1:p) L=L+Z_xqtls[,i]%*%t(Z_xqtls[,i])
  L=L/(m*p)
  mu=tr(L%*%R)
  variance=2*tr(L%*%R%*%L%*%R)
  beta=mu/variance
  alpha=mu*beta
  stat=t(z)%*%L%*%z
  pval=pgamma(c(stat),shape=alpha,rate=beta,lower.tail=FALSE)
  list(pval=pval,stat=stat,alpha=alpha,beta=beta,mu=mu,variance=variance)
}
########################################################################################
ldref='/home/lorincn/beegfs/lorincn/data/reference_panels/1kg.v3/EUR'
gwasfp='/home/lorincn/beegfs/lorincn/data/phenotype/AD/AD_Bellenguez_2022N_b37.tsv.gz'
########################################################################################
# bim
bim=fread(paste0(ldref,'.bim')) %>%
  select(rsid=V2,a1=V5)
dupes=bim %>% group_by(rsid) %>% summarise(n=n()) %>% filter(n>1) %>% pull(rsid)
bim=bim %>% filter(!(rsid %in% dupes))
# gene positions
gpos=fread('/home/lorincn/isilon/Cheng-Noah/manuscripts/druggable_genes/data/genepositions_hg19_long.txt')
# AD data
gwas=fread(gwasfp) %>% filter(rsid %in% bim$rsid)
## genome-wide xGenT
chrs=1:22
for(cc in chrs) {
  cat('Starting chr',cc,'\n')
  setwd('/home/lorincn/beegfs/lorincn/data/gene/pQTL/MetaMayo')
  pqtl_fpin=paste0('chr',cc,'_raw_Cortex_MetaMayo_eQTL.Rds')
  pqtl_chr=readRDS(pqtl_fpin) %>%
    rename(z=t_statistic)
  gwas_chr=gwas %>% 
    filter(chr==cc) %>%
    inner_join(bim,by='rsid') %>%
    mutate(beta=ifelse(effect_allele==a1,beta,-beta))
  gene_chr=gpos %>% filter(chr==cc)
  genes1=(gene_chr %>% pull(symbol) %>% unique())
  genes2=(pqtl_chr %>% pull(symbol) %>% unique())
  genes=intersect(genes1,genes2)
  chrdf=data.frame()
  for(i in 1:length(genes)) {
    cat(' Starting gene',genes[i],'\n')
    gdfi=gene_chr %>% filter(symbol==genes[i]) 
    starti=(gdfi %>% pull(start))[1]
    endi=(gdfi %>% pull(end))[1]
    if((endi-starti)>1e6) {cat('  Skipping gene',genes[i],' (gene too big)\n');next}
    mergedi=inner_join(
      gwas_chr %>% 
        filter(position>(starti-5e4),position<(endi+5e4)) %>%
        mutate(gwas_z=beta/se) %>%
        select(rsid,position,ea1=effect_allele,gwas_z),
      pqtl_chr %>% 
        filter(symbol==genes[i]) %>%
        select(rsid,ea2=effect_allele,pqtl_z=z),
      by='rsid') %>%
      na.omit() %>%
      arrange(position)
    dupes=mergedi %>% group_by(rsid) %>% summarise(n=n()) %>% filter(n>1) %>% pull(rsid)
    mergedi=mergedi %>% filter(!(rsid %in% dupes))
    if(nrow(mergedi)==0) {cat('  Skipping gene',genes[i],'\n');next}
    # LD matrix
    ld=get_ld(mergedi$rsid,ldref)
    ix=which(apply(ld,2,function(h) all(is.na(h))))
    if(length(ix)>0) {
      ld=ld[-ix,-ix]
      mergedi=mergedi[-ix,]
    }
    # test
    Z_xqtls=mergedi %>% pull(pqtl_z) %>% as.matrix()
    z=mergedi %>% pull(gwas_z)
    testres=xvegas(z,ld,Z_xqtls)
    toadd=as.data.frame(testres) %>%
      mutate(gene=genes[i]) %>%
      select(gene,pval,stat,null_shape=alpha,null_rate=beta,null_mean=mu,null_variance=variance)
    chrdf=rbind(chrdf,toadd)
  }
  # save
  setwd('/home/lorincn/isilon/Cheng-Noah/manuscripts/druggable_genes/data/xGenT_results')
  fpout=paste0('MetaMayo_cortex_pQTL_chr',cc,'.Rds')
  saveRDS(chrdf,fpout)
}



