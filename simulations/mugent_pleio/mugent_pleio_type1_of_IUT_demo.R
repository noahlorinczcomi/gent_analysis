
# should be deflated for large p and small m
# should be inflated for large m and small p
alpha=0.05
at=function(alpha,meff,p) 1-(1-alpha^p)^meff
at(alpha=0.05,meff=30,p=5)
meff=1:100
ps=1:3
plot(meff,meff,type='n',ylim=c(0,1),xlab='effective number of independent SNPs',ylab='Type I error level',yaxt='n')
axis(side=2,at=c(0,0.25,0.5,0.75,1),labels=c(0,0.25,0.5,0.75,1))
cols=c('gold','orange','red')
for(i in ps) lines(meff,at(alpha,meff,ps[i]),col=cols[i],lwd=2)
abline(h=alpha,lty=3)
legend(x=75,y=0.7,title='# pops.',legend=ps,col=cols,lwd=rep(2,max(ps)),bg='white',cex=2/3)
text(x=90,y=0.08,label='0.05',cex=2/3)
# show corrected quantiles
at=function(alpha,meff,p) 1-(1-alpha^p)^meff
at(alpha=0.05,meff=30,p=5)
meff=1:100
ps=1:3
Qs=matrix(nr=length(meff),nc=length(ps))
for(i in ps) {
  atilde=at(alpha,meff,ps[i])
  Qs[,i]=qchisq(1-(1-(1-alpha)^(1/meff))^(1/ps[i]),1)
}
matplot(meff,Qs,col=cols,type='l',lwd=2,lty=1,xlab='effective number of independent SNPs',ylab='corrected quantile of IUT',yaxt='n')
axis(side=2,at=c(0,3,6,9,12),labels=c(0,3,6,9,12))
abline(h=qchisq(1-alpha,1),lty=3)
legend(x=75,y=11,title='# pops.',legend=ps,col=cols,lwd=rep(2,max(ps)),bg='white',cex=2/3)
text(x=70,y=qchisq(1-alpha,1)+1/3,label='nominal 1-0.05 quantile',cex=2/3)













