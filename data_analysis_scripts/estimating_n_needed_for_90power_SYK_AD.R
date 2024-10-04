beta0=-0.0439 # rs10512201 association with AD in Bellenguez et al (2022)
se0=0.0103
n0=450000
sigmaj=sqrt(n0)*se0
beta=beta0/sigmaj
se=1/sqrt(n0)
# now beta~N(a,1/n)
ns=seq(500000,1.5e6,1e4)
alpha=5e-8
q0=qchisq(1-5e-8,1)
power=c()
for(i in 1:length(ns)) {
	lambda=beta^2*ns[i]
	power[i]=pchisq(q0,1,lambda,lower.tail=FALSE)
}
plot(ns,power,type='l')
abline(h=0.9)
res=ns[which.min(abs(0.9-power))]
abline(v=res)
res
# function version
f=function(n) {
	alpha=5e-8
	q0=qchisq(1-5e-8,1)
	lambda=beta^2*n
	abs(0.9-pchisq(q0,1,lambda,lower.tail=FALSE))
}
optim(1e6,f,lower=5e5,upper=2e6,method='Brent')$par