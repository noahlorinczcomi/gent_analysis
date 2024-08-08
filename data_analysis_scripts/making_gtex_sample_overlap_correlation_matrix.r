rm(list=ls(all=TRUE))
library(data.table);library(dplyr)
##
load_nonsig_gtex_gwas=function(chri=1,usetissues='all',verbose=TRUE) {
  bim=fread('/home/lorincn/beegfs/lorincn/data/reference_panels/1kg.v3/TRANS.bim') %>% select(chr=V1,rsid=V2,a1=V5) %>% filter(chr==chri) %>% select(-chr)
  setwd('/home/lorincn/beegfs/lorincn/data/gene/eQTL/GTExV8')
  tissues=dir()
  if(!all(usetissues=='all')) tissues=usetissues
  xdf=data.frame()
  for(i in 1:length(tissues)) {
    setwd('/home/lorincn/beegfs/lorincn/data/gene/eQTL/GTExV8')
    setwd(tissues[i])
    tissuei=unlist(strsplit(tissues[i],'_'))
    tissuei=stringr::str_to_title(tissuei)
    tissuei=paste(tissuei,collapse=' ')
    if(verbose) cat(' Loading GTEx',tissuei,'\n')
    fp=sapply(dir(),function(h) grepl(paste0('chr',chri),h))
    fp=dir()[which(fp)]
    if(chri==1) fp=fp[which(sapply(fp,function(h) grepl('chr1.',h,fixed=TRUE)))]
    if(chri==2) fp=fp[which(sapply(fp,function(h) grepl('chr2.',h,fixed=TRUE)))]
    dati=fread(fp) %>% as.data.frame() %>% filter(CHR==chri)
    dati$phenotype_id=sapply(dati$phenotype_id,function(h) unlist(strsplit(h,'[.]'))[1])
    dati=dati %>%
      select(ensembl_id=phenotype_id,beta=slope,se=slope_se,rsid=rsID,position=POS_b37,effect_allele=ALT,other_allele=REF,maf=maf) %>%
      inner_join(bim %>% select(rsid,a1),by='rsid') %>% 
      mutate(beta=ifelse(effect_allele==a1,beta,-beta)) %>%
      mutate(pval=pchisq((beta/se)^2,1,lower.tail=FALSE)) %>%
      select(-a1)
    xdf=rbind(xdf, dati %>% mutate(tissue=tissuei))
  }
  return(xdf %>% distinct())
}
chr1=load_nonsig_gtex_gwas()
# non-sig SNPs are already selected
# filter to SNPs present for all genes
alltissues=unique(chr1$tissue)
all_tissue_df=chr1 %>% 
    group_by(ensembl_id,rsid) %>% 
    summarise(n=n()) %>% 
    ungroup() %>% 
    filter(n==length(alltissues))
# rm(chr1)
# randomly select 2 SNPs from each SNP-gene pair and make wide
ransnps=all_tissue_df %>% 
    as_tibble() %>%
    group_by(ensembl_id) %>% 
    mutate(nsnps=length(rsid)) %>%
    filter(nsnps>=10) %>%
    mutate(keep=rsid %in% sample(rsid,10,replace=F)) %>% 
    filter(keep) %>%
    mutate(snp_gene_pair=paste0(ensembl_id,':',rsid))
ransnps=ransnps %>% 
    group_by(rsid) %>% 
    summarise(n=n(),snp_gene_pair=snp_gene_pair[1]) %>%
    ungroup() %>%
    filter(n==1)
df=chr1 %>%
    as_tibble() %>%
    mutate(snp_gene_pair=paste0(ensembl_id,':',rsid)) %>%
    filter(snp_gene_pair %in% ransnps$snp_gene_pair) %>%
    mutate(z=beta/se) %>%
    select(rsid,z,tissue) %>%
    distinct() %>%
    tidyr::pivot_wider(values_from='z',names_from='tissue') %>%
    tidyr::drop_na()
Zmat=df %>% select(-rsid) %>% as.matrix()
keepix=apply(Zmat,1,function(h) all(h^2<qchisq(1-0.01,1)))
R=cor(Zmat[which(keepix),],method='spearman')
rn=rownames(R)
setwd('/home/lorincn/isilon/Cheng-Noah/manuscripts/druggable_genes/data/')
saveRDS(R,'gtex_sample_overlap_correlations.Rds')

