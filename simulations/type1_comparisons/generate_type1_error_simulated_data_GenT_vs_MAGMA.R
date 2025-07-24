
# GenT vs MAGMA Type I error
# Steps:
## 1) (R) Generate simulated GWAS summary statistics under null hypothesis
##  a) (R) Identify all SNP-gene pairs
##  b) (R) Extract LD matrices for each gene
##  c) (R) Draw effect sizes from multivariate normal distribution with known LD
##  d) (R) Save file of GWAS summary statistics
## 2) (HPC) Perform MAGMA and GenT
## 3) (R) Read in results and plot etc.

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## 1) (R) Generate simulated GWAS summary statistics under null hypothesis
library(data.table);library(dplyr)
library(mvnfast,lib='~/isilon/Cheng-Noah/Rpkgs')
population='EUR'

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
usedsnps=c()
for(cc in 1:22) {
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
                usedsnps=c(usedsnps,snps_around)
                Ri=Ri[snps_around,snps_around]
                mi=nrow(Ri)
                # need to make sure Ri is positive definite
                lambda=eigen(Ri)$values
                k=0;Ii=diag(mi)
                while(min(lambda)<1e-4) {
                    k=k+1
                    Ri=Ri*(1-k/20)+(k/20)*Ii
                    lambda=eigen(Ri)$values
                }
                # randomly draw GWAS summary statistics for this gene
                # just going to draw once because 18000 genes are like 18000 simulation replicates
                zi=rmvn(1,rep(0,mi),Ri) %>% c()
                # save (generic format)
                toadd=data.frame(rsid=snps_around,z=zi,p=pchisq(zi^2,1,lower.tail=FALSE)) %>%
                    left_join(bim_chr %>% filter(rsid %in% snps_around) %>% select(rsid,chr,position),by='rsid')
                # MAGMA format
                magma_dfi=toadd %>% mutate(N=50000) %>% select(SNP=rsid,P=p,N,Z=z,CHR=chr,BP=position) %>% mutate(gene=genes[i])
                # GenT format
                gent_dfi=toadd %>% 
                    left_join(rndf,by='rsid') %>%
                    mutate(effect_allele=a1) %>%
                    mutate(beta=z,se=1) %>% # making generic (beta,se) in case GenT expects it
                    select(rsid,chr,position,effect_allele,beta,se,z,pval=p) %>%
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
write.table(magma_df,'MAGMA/type1.txt',row.names=FALSE,quote=FALSE,sep='\t')
write.table(gent_df,'GenT/type1.txt',row.names=FALSE,quote=FALSE,sep='\t')
