library(data.table);library(dplyr)
library(mvnfast,lib='~/isilon/Cheng-Noah/Rpkgs')
population='EUR'
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# functions
find_n=function(h2,M,R,target=0.5,alpha=0.05/18160) {
  m=nrow(R) # number of tested SNPs
  ix=round(m/2) # index of causal SNP
  mu0=rep(0,m) # empty mean vector
  mu0[ix]=1 # this will be scaled by sqrt(N*h2/M) later
  E0=m # null expectation
  V0=sum(R^2) # null variance
  rate0=E0/V0 # null Gamma rate parameter
  shape0=E0*rate0 # null Gamma shape parameter
  q0=qgamma(1-alpha,shape=shape0,rate=rate0) # critical value
  powerfun=function(N) { # function to compute power given N
    muN=sqrt(N*h2/M)*mu0
    E1=E0+sum(muN^2)
    V1=V0+4*(t(muN)%*%R%*%muN)
    rate1=E1/V1
    shape1=E1*rate1
    pgamma(q0,shape=shape1,rate=rate1,lower.tail=FALSE)
  }
  n0=2e5 # starting value
  k=0; eps=0.01; error=eps+1
  while (k<15 & error>eps) {
    k=k+1
    ns=seq(n0*(1-2/k),n0*(1+5/k),length.out=100)
    ns=ns[ns>1e3 & ns<2e6]
    ps=sapply(ns,powerfun)
    ix=which.min(abs(target-ps))
    n0=ns[ix]
    error=abs(ps[ix]-target)
  }
  n0  
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
### STEPS
# GenT vs MAGMA Power
# Steps:
## 1) (R) Generate simulated GWAS summary statistics under alternative hypothesis
##  a) (R) Identify all SNP-gene pairs
##  b) (R) Extract LD matrices for each gene
##  c) (R) Draw effect sizes from multivariate normal distribution with known LD
##  d) (R) Save file of GWAS summary statistics
## 2) (HPC) Perform MAGMA and GenT
## 3) (R) Read in results and plot etc.
### NOTES
# all causal SNPs will have the same effect size
# all genes will have one causal SNP
# (all genes are causal)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

## 1) (R) Generate simulated GWAS summary statistics under alternative hypothesis
h2=0.2 # total SNP heritability of trait
M=18160 # number of causal SNPs
target=0.5 # target power to achieve
alpha=0.05/18160 # Type I error rate

# identify positions of gene windows
setwd('~/isilon/Cheng-Noah/manuscripts/druggable_genes')
gp=fread('data/genepositions_hg19.txt') %>% 
  as_tibble() %>% 
  filter((end-start)<1e6) %>%
  mutate(window_start=start-5e4,window_end=end+5e4)

# identify all SNPs in each gene window
setwd('~/isilon/Cheng-Noah/reference_panels/1kg.v3')
bim=fread(paste0(population,'.bim')) %>% 
  as_tibble() %>%
  rename(chr=V1,rsid=V2,x=V3,position=V4,a1=V5,a2=V6)

# make vector of gene-specific LD matrix filenames from reference data
magma_df=gent_df=data.frame()
for(cc in 1:22) {
  usedsnps=c()
  cat('chr',cc,'\n')
  gene_chr=gp %>% filter(chr==cc)
  bim_chr=bim %>% filter(chr==cc)
  genes=gene_chr %>% pull(symbol) %>% unique()
  for(i in 1:length(genes)) {
    if(i==length(genes)) cat(genes[i],'\n',sep='') else cat(genes[i],', ',sep='')
    tryCatch(
      {
        # identify SNPs for this gene and their LD matrix (pre-computed)
        setwd(paste0('~/isilon/Cheng-Noah/reference_data/ld_matrices/',population))
        setwd(paste0('chr',cc))
        gene_fps=dir()
        ix=which(gene_fps==paste0(genes[i],'.Rds'))
        if(length(ix)==0) next
        Ri=readRDS(gene_fps[ix])
        gene_chri=gene_chr %>% filter(symbol==genes[i]) %>% head(.,1)
        snps_around=bim_chr %>% 
          filter(position>gene_chri$window_start & position<gene_chri$window_end) %>% 
          arrange(position) %>%
          pull(rsid)
        rn=rownames(Ri)
        rndf=sapply(rn,\(.) unlist(strsplit(.,'_'))) %>% t() %>% as.data.frame() %>% rename(rsid=V1,a1=V2)
        rn=sapply(rn,\(.) unlist(strsplit(.,'_'))[1])
        rn=unname(rn)
        rownames(Ri)=colnames(Ri)=rn
        snps_around=intersect(snps_around,rn)
        snps_around=snps_around[!(snps_around %in% usedsnps)]
        usedsnps=unique(c(usedsnps,snps_around))
        Ri=Ri[snps_around,snps_around]
        mi=nrow(Ri)
        # need to make sure Ri is positive definite
        Ri=0.95*Ri+0.05*diag(nrow(Ri))
        # randomly draw GWAS summary statistics for this gene
        # just going to draw once because 18000 genes are like 18000 simulation replicates
        Ni=find_n(h2,M,Ri,target=target,alpha=alpha)
        b=rep(0,mi) # true effect sizes / sqrt(n)
        b[floor(mi/2)]=sqrt(Ni*h2/M)
        zi=rmvn(1,b,Ri) %>% c()
        # save (generic format)
        toadd=data.frame(rsid=snps_around,z=zi,p=pchisq(zi^2,1,lower.tail=FALSE)) %>%
          left_join(bim_chr %>% filter(rsid %in% snps_around) %>% select(rsid,chr,position,a1),by='rsid')
        # MAGMA format
        magma_dfi=toadd %>% 
          mutate(N=50000) %>% 
          select(SNP=rsid,P=p,N,Z=z,CHR=chr,BP=position,effect_allele=a1) %>% 
          mutate(gene=genes[i])
        # GenT format
        gent_dfi=toadd %>% 
          mutate(beta=z,se=1) %>% # making generic (beta,se) in case GenT expects it
          select(rsid,chr,position,effect_allele=a1,beta,se,z,pval=p) %>%
          mutate(gene=genes[i])
        # bind to running data frames
        magma_df=rbind(magma_df,magma_dfi) %>% distinct(SNP,.keep_all=TRUE)
        gent_df=rbind(gent_df,gent_dfi) %>% distinct(rsid,.keep_all=TRUE)
      },
      error=function(x) NA
    )
  }
}

# save data sets
setwd('/home/lorincn/isilon/Cheng-Noah/manuscripts/druggable_genes/MAGMA_simulations/GWAS')
write.table(magma_df,'MAGMA/power.txt',row.names=FALSE,quote=FALSE,sep='\t')
write.table(gent_df,'GenT/power.txt',row.names=FALSE,quote=FALSE,sep='\t')







