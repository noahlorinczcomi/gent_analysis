rm(list=ls(all=TRUE))
library(data.table);library(dplyr);library(ggplot2);library(magrittr);library(tidyr)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# WILL NEED TO RE-DO AFTER NEW GENT RESULTS ARE FINISHED. ISILON WAS INTERRUPTED
# WHEN I WAS RUNNING IT AND IT ONLY DID 7K OF THE 18K GENES.
# NEW RESULTS SAVED TO:
#  ~/isilon/Cheng-Noah/Gent_results/NEWNEW_GenT_TypeI.Rds
#  or
#  ~/isilon/Cheng-Noah/Gent_results/AD_2025-07-23.Rds
# slurm job: 2980750
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
setwd('/mnt/isilon/w_gmi/chengflab/Cheng-Noah/manuscripts/druggable_genes/MAGMA_simulations')
magma_df=fread('output/type1_error/magma.genes.out') %>%
  mutate(chr=as.numeric(CHR)) %>%
  as_tibble() %>%
  mutate(test='MAGMA')
gent_df=readRDS('output/type1_error/NEWNEW_GenT_TypeI_2025-07-23.Rds') %>%
  na.omit() %>%
  mutate(chr=as.numeric(chr)) %>%
  as_tibble() %>%
  mutate(test='GenT')
merged=bind_rows(
  magma_df %>% select(gene=SYMBOL,chr=CHR,start=START,m=NSNPS,pval=P,test),
  gent_df %>% select(gene,chr,start=gene_start,m,pval,test)
)

# inflation
infl=function(data,type='MAGMA') {
  if(type=='MAGMA') {
    stat=data$ZSTAT^2
    q0=qchisq(0.5,1)
    median(stat)/q0
  } else {
    stat=data %$% qgamma(pval,shape=shape,rate=rate,lower.tail=FALSE)
    q0=data %$% qgamma(1/2,shape=shape,rate=rate)
    median(stat/q0)
  }
}
(magma_infl=infl(magma_df,'MAGMA'))
(gent_infl=infl(gent_df,'GenT'))

# QQ plot
labdf=data.frame(
  test=c('GenT','MAGMA'),
  x=c(1,1),
  y=c(4,4),
  label=c(
    paste0('inflation=',format(round(gent_infl,2),nsmall=2)),
    paste0('inflation=',format(round(magma_infl,2),nsmall=2)))
)

merged %>%
  na.omit() %>%
  group_by(test) %>%
  arrange(pval) %>%
  mutate(pp=ppoints(n())) %>%
  ungroup() %>%
  mutate(tpp=-log10(pp),tp=-log10(pval)) %>%
  ggplot(aes(tpp,tp,color=test)) +
  geom_abline(intercept=0,slope=1) +
  facet_wrap(~test) +
  geom_point() +
  scale_color_manual(NULL,values=c('#00B686','#234B84')) +
  geom_text(aes(x,y,label=label),data=labdf,inherit.aes=FALSE) +
  theme_bw() +
  theme(legend.position='none',
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank()) +
  guides(colour=guide_legend(override.aes=list(size=2))) +
  labs(x=expression('expected -log'[10]*'(P-value)'),
       y=expression('observed -log'[10]*'(P-value)'))
setwd('/mnt/isilon/w_gmi/chengflab/Cheng-Noah/manuscripts/druggable_genes/MAGMA_simulations/plots')
ggsave('Type1_QQplot.pdf',width=5.55,height=3.31)
ggsave('Type1_QQplot.png',width=5.55,height=3.31,units='in',dpi=400)

# Manhattan plot
chrdf=merged %>%
  na.omit() %>%
  group_by(test,chr) %>%
  mutate(newbp=chr+start/max(start)) %>%
  ungroup() %>%
  mutate(chrcol=chr %in% seq(1,22,2))
axisdf=chrdf %>% group_by(chr) %>% summarise(newbp=median(newbp)) %>% ungroup()
chrdf %>%
  mutate(tp=ifelse(test=='MAGMA',log10(pval),-log10(pval))) %>%
  ggplot(aes(newbp,tp,color=chrcol)) +
  geom_point() +
  geom_hline(yintercept=0,lwd=1/3) +
  geom_hline(yintercept=log10(0.05/12727),lwd=1/3,linetype='dashed') +
  geom_hline(yintercept=-log10(0.05/12727),lwd=1/3,linetype='dashed') +
  scale_color_manual(values=c('skyblue','cornflowerblue')) +
  theme_bw() +
  theme(legend.position='none',
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank()) +
  scale_x_continuous(breaks=axisdf$newbp,labels=axisdf$chr) +
  scale_y_continuous(breaks=seq(-15,15,5),labels=abs(seq(-15,15,5))) +
  coord_cartesian(ylim=c(-15,15)) +
  annotate('text',x=11,y=15,label='GenT',fontface='bold') +
  annotate('text',x=11,y=-15,label='MAGMA',fontface='bold') +
  labs(x='chromosome',
       y=expression('-log'[10]*'(gene-based test P-value)'))
setwd('/mnt/isilon/w_gmi/chengflab/Cheng-Noah/manuscripts/druggable_genes/MAGMA_simulations/plots')
ggsave('Type1_Manhattan.pdf',width=5.75,height=3.31)
ggsave('Type1_Manhattan.png',width=5.75,height=3.31,units='in',dpi=400)

# actual false positive rates
levels=seq(0.05/12727,0.05,length.out=50)
res=data.frame()
for(i in 1:length(levels)) {
  toadd=merged %>%
    group_by(test) %>%
    summarise(power=mean(pval<levels[i]),
              n=n()) %>%
    ungroup() %>%
    pivot_wider(values_from=c(power,n),names_from=test) %>%
    mutate(level=levels[i])
  res=rbind(res,toadd)
}
rr=range(c(res$power_GenT,res$power_MAGMA,res$level))
reslong=res %>%
  mutate(se_GenT=sqrt(power_MAGMA*(1-power_MAGMA)/n_MAGMA),
         se_MAGMA=sqrt(power_MAGMA*(1-power_MAGMA)/n_MAGMA),
         low_GenT=power_GenT-2*sqrt(power_GenT*(1-power_GenT)/n_GenT),
         high_GenT=power_GenT+2*sqrt(power_GenT*(1-power_GenT)/n_GenT),
         low_MAGMA=power_MAGMA-2*se_MAGMA,
         high_MAGMA=power_MAGMA+2*se_GenT) %>%
  pivot_longer(cols=c(power_GenT,power_MAGMA)) %>%
  rename(test=name) %>%
  select(test,level,value,contains(c('low','high')))
reslong %>%
  ggplot(aes(level,value,color=test)) +
  geom_abline(intercept=0,slope=1,lwd=1/2,linetype='dashed',color='gray50') +
  geom_ribbon(aes(ymin=low_GenT,ymax=high_GenT),fill='#00B686',alpha=1/5,color=NA) +
  geom_ribbon(aes(ymin=low_MAGMA,ymax=high_MAGMA),fill='#234B84',alpha=1/5,color=NA) +
  geom_path(lwd=1) +
  scale_color_manual(NULL,labels=c('GenT','MAGMA'),values=c('#00B686','#234B84')) +
  lims(x=rr,y=rr) +
  theme_bw() +
  theme(legend.position=c(0.305,0.862),
        legend.background=element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank()) +
  coord_cartesian(xlim=c(0,0.01),ylim=c(0,0.01)) +
  scale_x_continuous(breaks=c(0,0.005,0.01),labels=c('0','0.005','0.01')) +
  scale_y_continuous(breaks=c(0,0.005,0.01),labels=c('0','0.005','0.01')) +
  labs(x='target type 1 error rate',y='achieved type 1 error rate')
setwd('/mnt/isilon/w_gmi/chengflab/Cheng-Noah/manuscripts/druggable_genes/MAGMA_simulations/plots')
ggsave('Type1_level_maintext.pdf',width=2.3,height=2.3) # main text
ggsave('Type1_level_maintext.png',width=2.3,height=2.3,units='in',dpi=400) # main text
ggsave('Type1_level_supplement.pdf',width=3.3,height=3.3) # Supplement
ggsave('Type1_level_supplement.png',width=3.3,height=3.3,units='in',dpi=400) # Supplement
















