library(data.table);library(dplyr);library(ggplot2)
setwd('/home/lorincn/isilon/Cheng-Noah/software/ldsc/nlc_ldscores/EUR')
lddf=fread('EUR.l2.ldscore.gz') %>%
  rename(rsid=SNP,chr=CHR,position=BP,ldscore=L2) %>%
  as_tibble()
setwd('/home/lorincn/beegfs/lorincn/data')
bim=fread('reference_panels/1kg.v3/EUR.bim') %>% select(CHR=V1,rsid=V2)
gwas=fread('phenotype/T2D/1kgSNPs_Suzuki_EUR_Metal_LDSC-CORR_Neff.v2.txt.gz')
gwas=gwas %>% filter(rsid %in% bim$rsid)
rm(bim)
gwas=gwas %>% filter(rsid %in% lddf$rsid)
gwasn=746120.8
gwas=inner_join(gwas,lddf %>% select(rsid,ldscore),by='rsid') %>% select(rsid,chr,position,beta,se,ldscore)
gp=fread('genepositions_hg19.txt')
sizes=c(10,50,100) # in Kb
resdf=data.frame()
for(cc in 1:22) {
  cat('chr',cc,'\n')
  gwaschr=gwas %>% filter(chr==cc)
  gpchr=gp %>% filter(chr==cc)
  genes=unique(gpchr$symbol)
  for(i in 1:length(genes)) {
    starti=gpchr$start[i]
    endi=gpchr$end[i]
    genei=gpchr$symbol[i]
    if(i%%round(length(genes)*0.1)==0) cat(' ',round(100*i/length(genes)),'% complete\n',sep='')
    for(j in 1:length(sizes)) {
      genedf=gwaschr %>% filter(position>(starti-1e3*sizes[j]),position<(endi+1e3*sizes[j]))
      nsnps=nrow(genedf)
      if(nsnps>20) {
        stat=(genedf$beta/genedf$se)^2
        ldscores=genedf$ldscore
        fit=lm.fit(cbind(1,gwasn*ldscores),stat)
        h2m=coef(fit)[2]
      } else {
        h2m=NA
      }
      toadd=data.frame(gene=genei,kbwindowsize=sizes[j],nsnps=nsnps,h2m=h2m,chr=gpchr$chr[i],start=starti,end=endi)
      resdf=rbind(resdf,toadd)
    }
  }
}
saveRDS(resdf,'different_SNP_window_sizes_number_SNPs.Rds')
resdf=readRDS('different_SNP_window_sizes_number_SNPs.Rds')
resdf %>%
  group_by(gene) %>%
  summarise(x=max(nsnps)/min(nsnps),
            y=max(nsnps)-min(nsnps)) %>%
  na.omit()

sizes=sort(unique(resdf$kbwindowsize))
resdf=resdf %>%
  mutate(kbwindowsize=factor(kbwindowsize),
         kbwindowsize=recode_factor(kbwindowsize,
                                    '10'='10Kb window',
                                    '50'='50Kb window',
                                    '100'='100Kb window'))
# SNP counts
resdf %>%
  filter(nsnps<1e4) %>%
  ggplot(aes(x=nsnps)) +
  geom_histogram(bins=100,fill='dodgerblue') +
  geom_vline(aes(xintercept=x),
             data=. %>% group_by(kbwindowsize) %>% summarise(x=median(nsnps)) %>% ungroup(),color='firebrick') +
  geom_text(aes(x=x+100,y=5000,label=round(x,2)),
            angle=270,color='firebrick',
            data=. %>% group_by(kbwindowsize) %>% summarise(x=median(nsnps)) %>% ungroup(),color='firebrick') +
  facet_wrap(~kbwindowsize,ncol=1) +
  coord_cartesian(xlim=c(0,3e3)) +
  labs(x='number of SNPs in window',y='count')
# log h2
resdf %>%
  filter(nsnps<1e4) %>%
  mutate(h2m=ifelse(h2m<0,0,h2m)) %>%
  ggplot(aes(x=log(h2m))) +
  geom_histogram(bins=100,fill='dodgerblue') +
  geom_vline(aes(xintercept=x),
             data=. %>% group_by(kbwindowsize) %>% summarise(x=median(log(h2m),na.rm=T)) %>% ungroup(),color='firebrick') +
  geom_text(aes(x=x-1/2,y=300,label=round(x,2)),
            angle=270,color='firebrick',
            data=. %>% group_by(kbwindowsize) %>% summarise(x=median(log(h2m),na.rm=T)) %>% ungroup(),color='firebrick') +
  facet_wrap(~kbwindowsize,ncol=1) +
  labs(x='log of estimated scaled SNP heritability',y='count')
# SNP-averaged log h2
resdf %>%
  filter(nsnps<1e4) %>%
  mutate(h2m=ifelse(h2m<0,0,h2m)) %>%
  ggplot(aes(x=log(h2m/nsnps))) +
  geom_histogram(bins=100,fill='dodgerblue') +
  geom_vline(aes(xintercept=x),
             data=. %>% group_by(kbwindowsize) %>% summarise(x=median(log(h2m/nsnps),na.rm=T)) %>% ungroup(),color='firebrick') +
  geom_text(aes(x=x-1/2,y=300,label=round(x,2)),
            angle=270,color='firebrick',
            data=. %>% group_by(kbwindowsize) %>% summarise(x=median(log(h2m/nsnps),na.rm=T)) %>% ungroup(),color='firebrick') +
  facet_wrap(~kbwindowsize,ncol=1) +
  labs(x='log of estimated SNP-averaged scaled SNP heritability',y='count')