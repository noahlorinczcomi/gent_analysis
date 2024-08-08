#### variance of vec(Z'Z)
rm(list=ls(all=TRUE))
library(mvnfast);library(corrplot)
source('simulations/mugent/functions.R')
####
m=100
p=3
niter=10000
LD=ar1(m,0.5)
ldlist=lapply(1:p,function(h) LD)
K=kronecker(diag(p),LD)
z=rmvn(niter,rep(0,m*p),K)
H=matrix(0,p^2,p^2)
for(iter in 1:niter) {
  Zi=c(z[iter,])
  Zi=matrix(Zi,m,p,byrow=FALSE)
  v=c(t(Zi)%*%Zi/m)
  H=H+v%*%t(v)
}
H=H/niter
g=expand.grid(1:p,1:p)
rn=paste0(g[,1],',',g[,2])
rownames(H)=colnames(H)=rn
# plot
vm=varmat(p,ldlist)
rownames(vm)=colnames(vm)=rn
corrplot(vm,
         main=paste0(p,' populations'),
         mar=c(0,0,2,0),
         tl.cex=4/5.5,
         tl.col='black',
         # cl.ratio=0.25,
         cl.ratio=0.35,
         addgrid.col='black',
         method='color',
         is.corr=FALSE,
         col=colorRampPalette(c("#FFF6BA","#3B9AB2","#F74F4F"))(50))
## covariance between (z_s^\top z_s) and (z_k^\top z_k) is 0:
# S=matrix(nr=niter,nc=p)
# for(iter in 1:niter) {
#   Zi=c(z[iter,])
#   Zi=matrix(Zi,m,p,byrow=FALSE)
#   S[iter,]=colSums(Zi^2)
# }
# cor(S)