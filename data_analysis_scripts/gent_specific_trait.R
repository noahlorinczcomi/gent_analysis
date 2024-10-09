rm(list=ls(all=TRUE))
library(data.table);library(dplyr);library(stringr)
source('/home/lorincn/isilon/Cheng-Noah/software/corefunctions/functions.R')
##############################################################################
<GWAS filepath (w/ extension) here>
<LD reference (w/o extension) here>
savedir='/home/lorincn/isilon/Cheng-Noah/GenT_results'
###########
genepositions=fread('/home/lorincn/isilon/Cheng-Noah/manuscripts/druggable_genes/data/genepositions_hg19_long.txt')
# analysis
trait_=unlist(strsplit(gwas_fp,'/'))
trait=tail(trait_,2)[1]
if(trait=='AI_derived_aging') trait=unlist(strsplit(tail(trait_,1),'[.]'))[1]
if(trait=='FunctionalConnectivity') {
  trait_=unlist(strsplit(gwas_fp,'/'))
  trait=tail(trait_,1)
  trait=unlist(strsplit(trait,'[.]'))[1]
}
if(trait=='COVID') trait='COVIDseverity'
pop=unlist(strsplit(ldref,'/'))
pop=tail(pop,1)
trait=paste0(trait,'_',pop)
gwas=fread(gwas_fp) %>% 
  as.data.frame() %>%
  mutate(position=as.numeric(position),
         beta=as.numeric(beta),
         se=as.numeric(se)) %>%
  na.omit()
genepositions=genepositions %>% filter((end-start)<1e6)
gw_gent(gwas,genepositions,BPwindow=5e4,gene_colname='symbol',ldref,save_dir=savedir,save_name=trait)

