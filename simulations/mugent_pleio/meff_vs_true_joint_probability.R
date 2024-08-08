rm(list=ls(all=T))
library(mvnfast);library(mvtnorm)
m=30
rhos=c(0.5,0.75,0.9)
taus=exp(seq(log(1/1000000),log(1/10),length.out=20))
niter=1000
ACTUAL=EST1=EST2=EST3=matrix(nr=length(rhos),nc=length(taus))
for(i in 1:length(rhos)) {
  LD=ar1(m,rhos[i]) # AR(1)
  # LD=matrix(rhos[i],m,m);diag(LD)=1 # compound symmetry
  meff=meff_fun(LD)
  z=rmvn(niter,rep(0,m),LD)
  for(j in 1:length(taus)) {
    check=c()
    for(iter in 1:niter) {
      zi=z[iter,]
      check[iter]=all(zi^2>taus[j])
    }
    ACTUAL[i,j]=sum(check)/niter
    # using meff approximation
    EST1[i,j]=pchisq(taus[j],1,lower.tail=FALSE)^(2*meff)
    # using multivariate normal density
    EST2[i,j]=1-2*pmvnorm(lower=rep(taus[j],m),upper=rep(Inf,m),mean=rep(0,m),corr=LD)
    # using Monte Carlo
    EST3[i,j]=simfun(sqrt(taus[j]),rep(0,m),LD)
  }
}
## plot first
plot(0,1,type='n',pch=1/2,cex=1/2,xlim=c(0,1),ylim=c(0,1),xlab='true joint probability',ylab='approximated joint probability',xaxt='n',yaxt='n',main=expression('m'[eff.]*' approximation'))
axis(side=1,at=seq(0,1,0.2),labels=seq(0,1,0.2))
axis(side=2,at=seq(0,1,0.2),labels=seq(0,1,0.2))
abline(a=0,b=1,lty=2)
cols=c('gold','orange','brown')
for(i in 1:nrow(EST1)) points(ACTUAL[i,],EST1[i,],col=cols[i],type='b',pch=1/2,cex=2/3)
legend('bottomright',cex=2/3,title=expression('AR1('*rho*')'),legend=c(expression(rho*'=0.50'),expression(rho*'=0.75'),expression(rho*'=0.90')),pch=rep(1/2,3),col=cols)
## plot second
plot(0,1,type='n',pch=1/2,cex=1/2,xlim=c(0,1),ylim=c(0,1),xlab='true joint probability',ylab='approximated joint probability',xaxt='n',yaxt='n',main=expression('Probabilty density'))
axis(side=1,at=seq(0,1,0.2),labels=seq(0,1,0.2))
axis(side=2,at=seq(0,1,0.2),labels=seq(0,1,0.2))
abline(a=0,b=1,lty=2)
cols=c('#45EB40','deepskyblue','blue')
for(i in 1:nrow(EST2)) points(ACTUAL[i,],EST2[i,],col=cols[i],type='b',pch=1/2,cex=2/3)
legend('bottomright',cex=2/3,title=expression('AR1('*rho*')'),legend=c(expression(rho*'=0.50'),expression(rho*'=0.75'),expression(rho*'=0.90')),pch=rep(1/2,3),col=cols)
# plot third
plot(0,1,type='n',pch=1/2,cex=1/2,xlim=c(0,1),ylim=c(0,1),xlab='true joint probability',ylab='approximated joint probability',xaxt='n',yaxt='n',main=expression('Monte Carlo'))
axis(side=1,at=seq(0,1,0.2),labels=seq(0,1,0.2))
axis(side=2,at=seq(0,1,0.2),labels=seq(0,1,0.2))
abline(a=0,b=1,lty=2)
cols=c('#45EB40','deepskyblue','blue')
for(i in 1:nrow(EST3)) points(ACTUAL[i,],EST3[i,],col=cols[i],type='b',pch=1/2,cex=2/3)
legend('bottomright',cex=2/3,title=expression('AR1('*rho*')'),legend=c(expression(rho*'=0.50'),expression(rho*'=0.75'),expression(rho*'=0.90')),pch=rep(1/2,3),col=cols)






