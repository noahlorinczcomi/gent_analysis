# power as more SNPs are included ####
ms=c(5:200,300,400,500,600,700,800,900,1000,1250,1500,1750,2000) # number of tested SNPs
alpha=0.05 # Type I error rate
m0=3 # number causal SNPs
h2=0.0005 # h2 explained by gene
ngwas=30000 # GWAS sample size
cortypes=c('CS','AR1') # CS, AR1
rhos=c(0,0.3,0.5,0.9)
POWER=array(dim=c(length(ms),length(rhos),length(cortypes)))
for(i in 1:length(ms)) {
  m=ms[i]
  causalix=round(seq(1,m,length.out=m0+2))[-c(1,m0+2)]
  b=rep(0,m)
  b[1:m0]=sqrt(h2/m0)
  for(j in 1:length(rhos)) {
    rho=rhos[j]
    for(k in 1:length(cortypes)) {
      if(k==1) R=CS(m,rho) else R=ar1(m,rho)
      adj=c(h2/t(b)%*%R%*%b)
      b=b*sqrt(adj)
      beta=R%*%b
      mu=sqrt(ngwas)*beta
      e0=m
      v0=2*tr(R%*%R)
      e1=m+sum(mu^2)
      v1=v0+4*t(mu)%*%R%*%mu
      q0=qgamma(1-alpha,shape=e0^2/v0,rate=e0/v0)
      POWER[i,j,k]=pgamma(q0,shape=e1^2/v1,rate=e1/v1,lower.tail=FALSE)
    }
  }
}
df=data.frame(
  x=c(POWER),
  m=rep(rep(ms,length(rhos)),length(cortypes)),
  rho=rep(rep(rhos,each=length(ms)),length(cortypes)),
  cortype=rep(cortypes,each=length(ms)*length(rhos))
) %>%
  mutate(cortype=factor(cortype),
         cortype=recode_factor(cortype,
                               'AR1'='First-order Autoregressive LD',
                               'CS'='Compound Symmetry LD'))
cols=paletteer::paletteer_d("ggsci::lanonc_lancet")[1:length(unique(df$rho))]
ggplot(df,aes(x=m,y=x,color=factor(rho))) +
  geom_path(aes(group=factor(rho))) +
  facet_wrap(~cortype) +
  lims(y=c(0,1)) +
  geom_hline(yintercept=0.05,lwd=1/3,linetype='dashed') +
  scale_color_manual('LD correlation parameter',values=cols) +
  theme(legend.position='bottom') +
  labs(x='number of tested SNPs (only 3 are cauasl)',
       y='power of gene-based association test',
       title='SNP window size and power of gene-based association test')