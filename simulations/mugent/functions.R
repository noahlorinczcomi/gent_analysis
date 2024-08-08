ar1=function(n,rho=0.9) rho^toeplitz(0:(n-1))
tr=function(x) sum(diag(x))
## MuGenT
mugent=function(Z,ldlist) {
  # can actually be general for multiple populations
  Z=as.matrix(Z)
  p=ncol(Z);m=nrow(Z);j=rep(1,p)
  stat=c(t(j)%*%(t(Z)%*%Z/m)%*%j)
  EZ=t(j)%*%diag(ncol(Z))%*%j # under H0
  Kj=kronecker(t(j),t(j))
  VZ=tr(Kj%*%varmat(p,ldlist)%*%t(Kj))/m^2
  beta=EZ/VZ
  alpha=EZ*beta
  p=pgamma(stat,shape=alpha,rate=beta,lower.tail=FALSE)
  list(stat=stat,crit05=qgamma(0.95,shape=alpha,rate=beta),p=p)
}
# function to calculate the variance of vec[t(Z)%*%Z] under H0 of MuGenT
varmat=function(p,ldlist) {
  if(length(ldlist)!=p) stop(cat('The number of LD matrices in `ldlist` should be ',p,', not ',length(ldlist),'\n',sep=''))
  H=matrix(0,p^2,p^2)
  ix1=c(row(H[1:p,1:p]))
  ix2=c(col(H[1:p,1:p]))
  C=cbind(ix1,ix2)
  diags=seq(1,p^2,length.out=p) 
  for(i in 1:p^2) {
    for(j in 1:p^2) {
      left=C[i,]
      right=C[j,]
      both=c(left,right)
      boo1=length(unique(both))==2
      boo2=(length(unique(left))==2) & (length(unique(right))==2) 
      if(boo1 & boo2) H[i,j]=tr(ldlist[[left[1]]]%*%ldlist[[left[2]]])
      if(i==j & (i%in%diags)) H[i,j]=2*tr(ldlist[[left[1]]]%*%ldlist[[right[1]]])
    }
  }
  H
}