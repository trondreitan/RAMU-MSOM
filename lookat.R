# Look at a single model:

#model="global"
#model="only_gammadyn"
#model="f"
#model="f_sf"
#model="f_sf_top"
#model="f_sf_newprior"
#model="f_sf_regocc_interact"
#model="newfull_dropint"
#model="colony"
model="smalltop"

source("read_data.R")
source(sprintf("model_%s.R",model))
hyper=hyper.this
hyper$N=N
hyper$N.sp=N.sp
hyper$Nf=Nf
hyper$N.forms=N.forms
hyper$sp=sp
hyper$genera=genera
hyper$genus.nr=genus.nr
hyper$genusnames=genusnames
hyper$unid.genus.nr=unid.genus.nr
hyper$N.un=N.un
hyper$inter=inter

source("make_laplace_wrapper.R")
library(coda)

par0=init.params.flat(hyper)

# Collate (if multiple runs)
N.start=1
N.runs=100

ll=rep(NA,N.runs)
# all
f=list()
for(i in 1:N.runs)
{
  show(i)
  load(sprintf("%s_%02d.Rdata",model,i+N.start-1))
  if(class(Fit3)=="demonoid")
  {
    if(i==1)
      res=mcmc.list(mcmc(Fit3$Monitor))
    if(i>1)
      res[[i]]=mcmc(Fit3$Monitor)
    f[[i]]=Fit3
    ll[i]=Fit3$LML
  }
  
  if(class(Fit3)=="demonoid.hpc")
  {
   start=1
   if(i==1)
   {
     res=mcmc.list(mcmc(Fit3[[1]]$Monitor))
     start=2
     f=Fit3
   }
   if(i>1)
   {
     for(j in 1:length(Fit3))
       {
         k=length(f)
         show(c(i,j,k))
 	 f[[k+1]]=Fit3[[j]]
       }
   }
   cc=length(Fit3)
   for(j in 1:cc)
    ll[j+(i-1)*cc]=Fit3[[j]]$LML
   k=length(res)
   for(j in start:cc)
     res[[j+k-start+1]]=mcmc(Fit3[[j]]$Monitor)
  }
}

source("make_laplace_wrapper.R")
data=make.data(par0,hyper,y)
d=Combine(f,Data=data)
d$LML
# global 5000: -3173.639
# global 20000: -3177.00
# global 100000: -3177.164
# only_gammadyn 5000: -3129.581
# only_gammadyn 20000: -3131.034
# f 5000: -2940.114
# f 20000: -2947.108
# f_sf 400000: NA
# f_sf_top 100000: NA
# f_sf_newprior 100000: NA
# newfull_dropint 100x200000: NA 
# colony 100x50000: NA
# smalltop 100x25000: NA

res2=mcmc(d$Monitor)




plot(res)

gelman.diag(res)
# global 5000: 1.03
# global 20000: 1.01
# global 100000: 1
# f 5000: 1.21
# f 20000: 1.07
# f_sf 400000: 1.38
# f_sf_top 100000: 1.05
# f_sf_top 400000: 1.03
# f_sf_newprior 100000: 1.02
# f_sf_regocc_interact 5000: 1.34
# newfull_dropint 100x200000: 1.03
# colony 100x50000: 1.11
# smalltop 100x100000: 1.08
# smalltop_wideprior: 100x25000: 4.54
# smalltop_widerprior: 100x75000: 6.85



# 4 worst parameters:
g=gelman.diag(res)
index=which(g$psrf[,1]>=sort(g$psrf[,1],decreasing=TRUE)[4])
#index=which(attributes(d$Monitor)$dimnames[[2]]=="log.sd.sf.bin1")
res.index=mcmc.list()
for(i in 1:length(res))
{
    res.index[[i]]=mcmc(res[[i]][,index])
}
plot(res.index)

gelman.diag(res.index)
#                   Point est. Upper C.I.
#log.const.lambda22       1.08       1.10
#f.lambda2                1.08       1.10
#sf.occu138               1.11       1.15
#sf.occu139               1.12       1.15

r=rep(0,4)
for(i in 1:4)
 r[i]=acf(res2[,index[i]],plot=F)$acf[2]
(1+r)/(1-r)
# 10.324567 23.532180 19.758560  7.224077 std prior
# 17.804152  2.754934  5.014976  3.967280
# 11.690535 13.853561  2.386205  2.147592



# Store estimates and uncertainty:

par.mean=apply(res2,2,mean)
par.median=apply(res2,2,median)
par.sd=apply(res2,2,sd)

save(par.mean,file=sprintf("par_mean_%s.RData",model))
save(par.sd,file=sprintf("par_sd_%s.RData",model))
save(par.median,file=sprintf("par_median_%s.RData",model))








# Look at top parameters:

top=c(1:9,53:54,65:66,109,560)
names(par0)[top]

m=c(hyper$log.kappa.mu, hyper$logit.mu, hyper$logsd.mu,hyper$log.lambda.mu,
  rep(hyper$logsd.mu,10), hyper$logit.mu)
s=c(hyper$log.kappa.sd, hyper$logit.sd, hyper$logsd.sd,hyper$log.lambda.sd,
  rep(hyper$logsd.sd,10), hyper$logit.sd)

for(i in 1:length(top))
{
  pname=names(par0)[top[i]]
  png(sprintf("par_%s_%s.png", pname, model),height=1200, width=1600)
  par(cex=3)
  hist(res2[,top[i]],freq=F,xlab=pname,main=sprintf("Histogram of %s",pname))
  x=seq(min(res2[,top[i]])-1,max(res2[,top[i]])+1,0.01)
  lines(x,dnorm(x, m[i], s[i]))
  dev.off()

  print.srcref(sprintf("%18s: m=%5.2f (%5.2f)  s=%5.2f (%5.2f)  S=%7.3f",
         pname, mean(res2[,top[i]]), m[i], sd(res2[,top[i]]), s[i],
	 s[i]/sd(res2[,top[i]])))
}
#         log.kappa: m= 0.02 ( 0.00)  s= 0.07 (20.00)  S=266.948
#     mu.const.occu: m= 0.48 ( 0.00)  s= 0.99 (10.00)  S= 10.090
#    lsd.const.occu: m= 1.07 ( 0.00)  s= 0.28 ( 6.11)  S= 21.873
#   mu.const.lambda: m=-3.40 (-3.00)  s= 0.57 (10.00)  S= 17.591
#  lsd.const.lambda: m= 0.34 ( 0.00)  s= 0.19 ( 6.11)  S= 32.111
#    mu.lsd.sf.occu: m= 0.64 ( 0.00)  s= 1.04 ( 6.11)  S=  5.900
#   lsd.lsd.sf.occu: m=-3.42 ( 0.00)  s= 2.56 ( 6.11)  S=  2.388
#  mu.lsd.sf.lambda: m=-0.08 ( 0.00)  s= 0.20 ( 6.11)  S= 30.316
# lsd.lsd.sf.lambda: m=-4.36 ( 0.00)  s= 2.47 ( 6.11)  S=  2.474
#    logit.gamma.mu: m= 1.50 ( 0.00)  s= 0.36 ( 6.11)  S= 16.944
#   logit.gamma.lsd: m=-0.05 ( 0.00)  s= 0.30 ( 6.11)  S= 20.057
#     log.sd.f.occu: m=-3.22 ( 0.00)  s= 2.70 ( 6.11)  S=  2.261
#   log.sd.f.lambda: m= 0.05 ( 0.00)  s= 0.32 ( 6.11)  S= 18.995
#log.sd.sf.lambda22: m=-4.70 ( 0.00)  s= 2.74 ( 6.11)  S=  2.234
#    logit.reg.occu: m= 5.38 ( 0.00)  s= 5.00 (10.00)  S=  2.000

ident=function(x) x
lower=function(x) quantile(x,0.025)
upper=function(x) quantile(x,0.975)
trans=c(exp,ident,exp,ident,exp,ident,exp,ident,exp,ident,exp,exp,exp,exp,ilogit)
for(i in 1:length(top))
{
  pname=names(par0)[top[i]]
  png(sprintf("par_%s_%s.png", pname, model),height=1200, width=1600)
 
  print.srcref(sprintf("%18s: m=%5.2f (%5.2f-%5.2f)",
         pname, mean(trans[[i]](res2[,top[i]])), lower(trans[[i]](res2[,top[i]])),
	 upper(trans[[i]](res2[,top[i]])) ))
}
# For regular prior:
#         log.kappa: m= 1.01 ( 0.88- 1.16)
#     mu.const.occu: m= 0.47 (-1.05- 2.10)
#    lsd.const.occu: m= 2.90 ( 1.67- 4.75)
#   mu.const.lambda: m=-3.33 (-4.33--2.40)
#  lsd.const.lambda: m= 1.40 ( 0.98- 2.01)
#    mu.lsd.sf.occu: m= 0.78 (-0.56- 1.26)
#   lsd.lsd.sf.occu: m= 0.33 ( 0.03- 1.19)
#  mu.lsd.sf.lambda: m=-0.09 (-0.38- 0.19)
# lsd.lsd.sf.lambda: m= 0.20 ( 0.02- 0.54)
#    logit.gamma.mu: m= 1.47 ( 0.75- 2.17)
#   logit.gamma.lsd: m= 1.00 ( 0.55- 1.77)
#     log.sd.f.occu: m= 0.47 ( 0.03- 1.50)
#   log.sd.f.lambda: m= 1.06 ( 0.55- 1.89)
#log.sd.sf.lambda22: m= 0.23 ( 0.02- 0.69)
#    logit.reg.occu: m= 0.90 ( 0.75- 0.99)



# log.kappa
i=1
hist(res2[,i],freq=F)
x=seq(min(res2[,i])-1,max(res2[,i])+1,0.01)
lines(x,dnorm(x, hyper$log.kappa.mu, hyper$log.kappa.sd))
c(pnorm(min(res2[,i]), hyper$log.kappa.mu, hyper$log.kappa.sd),
  1-pnorm(max(res2[,i]), hyper$log.kappa.mu, hyper$log.kappa.sd))

mean(exp(res2[,i]))
# 1.015191

sd(res2[,i])/hyper$log.kappa.sd
# 0.003746048


# mu.const.occu
i=2
hist(res2[,i],freq=F)
x=seq(min(res2[,i])-1,max(res2[,i])+1,0.01)
lines(x,dnorm(x, hyper$logit.mu, hyper$logit.sd))
c(pnorm(min(res2[,i]), hyper$logit.mu, hyper$logit.sd),
  1-pnorm(max(res2[,i]), hyper$logit.mu, hyper$logit.sd))
c(pnorm(quantile(res2[,i],0.025), hyper$logit.mu, hyper$logit.sd),
  1-pnorm(quantile(res2[,i],0.975), hyper$logit.mu, hyper$logit.sd))

sd(res2[,i])/hyper$logit.sd
# 0.09911241


# lsd.const.occu
i=3
hist(res2[,i],freq=F)
x=seq(min(res2[,i])-1,max(res2[,i])+1,0.01)
lines(x,dnorm(x, hyper$logsd.mu, hyper$logsd.sd))
c(pnorm(min(res2[,i]), hyper$logsd.mu, hyper$logsd.sd),
  1-pnorm(max(res2[,i]), hyper$logsd.mu, hyper$logsd.sd))
c(pnorm(quantile(res2[,i],0.025), hyper$logsd.mu, hyper$logsd.sd),
  1-pnorm(quantile(res2[,i],0.975), hyper$logsd.mu, hyper$logsd.sd))

sd(res2[,i])/hyper$logsd.sd
# 0.04571877


# mu.const.lambda
i=4
hist(res2[,i],freq=F)
x=seq(min(res2[,i])-1,max(res2[,i])+1,0.01)
lines(x,dnorm(x, hyper$log.lambda.mu, hyper$log.lambda.sd))
c(pnorm(min(res2[,i]), hyper$log.lambda.mu, hyper$log.lambda.sd),
  1-pnorm(max(res2[,i]), hyper$log.lambda.mu, hyper$log.lambda.sd))
c(pnorm(quantile(res2[,i],0.025), hyper$logsd.mu, hyper$logsd.sd),
  1-pnorm(quantile(res2[,i],0.975), hyper$logsd.mu, hyper$logsd.sd))

sd(res2[,i])/hyper$logsd.sd
# 0.09298217

# lsd.const.lambda
i=5
hist(res2[,i],freq=F)
x=seq(min(res2[,i])-1,max(res2[,i])+1,0.01)
lines(x,dnorm(x, hyper$logsd.mu, hyper$logsd.sd))
c(pnorm(min(res2[,i]), hyper$logsd.mu, hyper$logsd.sd),
  1-pnorm(max(res2[,i]), hyper$logsd.mu, hyper$logsd.sd))
c(pnorm(quantile(res2[,i],0.025), hyper$logsd.mu, hyper$logsd.sd),
  1-pnorm(quantile(res2[,i],0.975), hyper$logsd.mu, hyper$logsd.sd))

sd(res2[,i])/hyper$logsd.sd
# 0.03114168


# mu.lsd.sf.occu
i=6
hist(res2[,i],freq=F)
x=seq(min(res2[,i])-1,max(res2[,i])+1,0.01)
lines(x,dnorm(x, hyper$logsd.mu, hyper$logsd.sd))
c(pnorm(min(res2[,i]), hyper$logsd.mu, hyper$logsd.sd),
  1-pnorm(max(res2[,i]), hyper$logsd.mu, hyper$logsd.sd))
c(pnorm(quantile(res2[,i],0.025), hyper$logsd.mu, hyper$logsd.sd),
  1-pnorm(quantile(res2[,i],0.975), hyper$logsd.mu, hyper$logsd.sd))

sd(res2[,i])/hyper$logsd.sd
# 0.169494

# lsd.lsd.sf.occu
i=7
hist(res2[,i],freq=F)
x=seq(min(res2[,i])-1,max(res2[,i])+1,0.01)
lines(x,dnorm(x, hyper$logsd.mu, hyper$logsd.sd))
c(pnorm(min(res2[,i]), hyper$logsd.mu, hyper$logsd.sd),
  1-pnorm(max(res2[,i]), hyper$logsd.mu, hyper$logsd.sd))
c(pnorm(quantile(res2[,i],0.025), hyper$logsd.mu, hyper$logsd.sd),
  1-pnorm(quantile(res2[,i],0.975), hyper$logsd.mu, hyper$logsd.sd))
# This one seems almost just truncated even for wider prior
sd(res2[,i])/hyper$logsd.sd
# 0.42

# mu.lsd.sf.lambda
i=8
hist(res2[,i],freq=F)
x=seq(min(res2[,i])-1,max(res2[,i])+1,0.01)
lines(x,dnorm(x, hyper$logsd.mu, hyper$logsd.sd))
c(pnorm(min(res2[,i]), hyper$logsd.mu, hyper$logsd.sd),
  1-pnorm(max(res2[,i]), hyper$logsd.mu, hyper$logsd.sd))
c(pnorm(quantile(res2[,i], 0.025), hyper$logsd.mu, hyper$logsd.sd),
  1-pnorm(quantile(res2[,i],0.975), hyper$logsd.mu, hyper$logsd.sd))

sd(res2[,i])/hyper$logsd.sd
# 0.033

# lsd.lsd.sf.lambda
i=9
hist(res2[,i],freq=F)
x=seq(min(res2[,i])-1,max(res2[,i])+1,0.01)
lines(x,dnorm(x, hyper$logsd.mu, hyper$logsd.sd))
c(pnorm(min(res2[,i]), hyper$logsd.mu, hyper$logsd.sd),
  1-pnorm(max(res2[,i]), hyper$logsd.mu, hyper$logsd.sd))
c(pnorm(quantile(res2[,i], 0.025), hyper$logsd.mu, hyper$logsd.sd),
  1-pnorm(quantile(res2[,i],0.975), hyper$logsd.mu, hyper$logsd.sd))

sd(res2[,i])/hyper$logsd.sd
# 0.40



# log.sd.f.occu
i=65
hist(res2[,i],freq=F, xlab=names(par0)[i], main=sprintf("Histrogram of %s",names(par0)[i]))
x=seq(min(res2[,i])-1,max(res2[,i])+1,0.01)
lines(x,dnorm(x, hyper$logsd.mu, hyper$logsd.sd))
c(pnorm(min(res2[,i]), hyper$logsd.mu, hyper$logsd.sd),
  1-pnorm(max(res2[,i]), hyper$logsd.mu, hyper$logsd.sd))
c(pnorm(quantile(res2[,i], 0.025), hyper$logsd.mu, hyper$logsd.sd),
  1-pnorm(quantile(res2[,i],0.975), hyper$logsd.mu, hyper$logsd.sd))

sd(res2[,i])/hyper$logsd.sd
# 0.4422628

# log.sd.f.lambda
i=66
hist(res2[,i],freq=F)
x=seq(min(res2[,i])-1,max(res2[,i])+1,0.01)
lines(x,dnorm(x, hyper$logsd.mu, hyper$logsd.sd))
c(pnorm(min(res2[,i]), hyper$logsd.mu, hyper$logsd.sd),
  1-pnorm(max(res2[,i]), hyper$logsd.mu, hyper$logsd.sd))
c(pnorm(quantile(res2[,i], 0.025), hyper$logsd.mu, hyper$logsd.sd),
  1-pnorm(quantile(res2[,i],0.975), hyper$logsd.mu, hyper$logsd.sd))

sd(res2[,i])/hyper$logsd.sd
# 0.05264587


# log.sd.sf.lambda22 (superspecies)
i=109
hist(res2[,i],freq=F)
x=seq(min(res2[,i])-1,max(res2[,i])+1,0.01)
lines(x,dnorm(x, hyper$logsd.mu, hyper$logsd.sd))
c(pnorm(min(res2[,i]), hyper$logsd.mu, hyper$logsd.sd),
  1-pnorm(max(res2[,i]), hyper$logsd.mu, hyper$logsd.sd))
c(pnorm(quantile(res2[,i], 0.025), hyper$logsd.mu, hyper$logsd.sd),
  1-pnorm(quantile(res2[,i],0.975), hyper$logsd.mu, hyper$logsd.sd))

sd(res2[,i])/hyper$logsd.sd
# 0.45

# logit.gamma.mu
i=53
hist(res2[,i],freq=F)
x=seq(min(res2[,i])-1,max(res2[,i])+1,0.01)
lines(x,dnorm(x, hyper$logsd.mu, hyper$logsd.sd))
c(pnorm(min(res2[,i]), hyper$logsd.mu, hyper$logsd.sd),
  1-pnorm(max(res2[,i]), hyper$logsd.mu, hyper$logsd.sd))
c(pnorm(quantile(res2[,i], 0.025), hyper$logsd.mu, hyper$logsd.sd),
  1-pnorm(quantile(res2[,i],0.975), hyper$logsd.mu, hyper$logsd.sd))

sd(res2[,i])/hyper$logsd.sd
# 0.059




mean(res2[,i])
# 1.467554

mean(ilogit(res2[,i]))
# 0.8069108

# Compare to data



# logit.gamma.log.sd
i=54
hist(res2[,i],freq=F)
x=seq(min(res2[,i])-1,max(res2[,i])+1,0.01)
lines(x,dnorm(x, hyper$logit.mu, hyper$logit.sd))
c(pnorm(min(res2[,i]), hyper$logit.mu, hyper$logit.sd),
  1-pnorm(max(res2[,i]), hyper$logit.mu, hyper$logit.sd))
c(pnorm(quantile(res2[,i], 0.025), hyper$logit.mu, hyper$logit.sd),
  1-pnorm(quantile(res2[,i],0.975), hyper$logit.mu, hyper$logit.sd))

mean(exp(res2[,i]))
# 0.9994813

x=ilogit(rnorm(50000,res2[,53],exp(res2[,54])))

sd(res2[,i])/hyper$logsd.sd
# 0.04985885





# logit.reg.occu
i=560
hist(res2[,i],freq=F)
x=seq(min(res2[,i])-1,max(res2[,i])+1,0.01)
lines(x,dnorm(x, hyper$logit.mu, hyper$logit.sd))
c(pnorm(min(res2[,i]), hyper$logit.mu, hyper$logit.sd),
  1-pnorm(max(res2[,i]), hyper$logit.mu, hyper$logit.sd))
c(pnorm(quantile(res2[,i], 0.025), hyper$logit.mu, hyper$logit.sd),
  1-pnorm(quantile(res2[,i],0.975), hyper$logit.mu, hyper$logit.sd))


mean(ilogit(res2[,i]))
# 0.8982944 # for standard prior

# Compare to raw estimate:
sum(no.presence[,1:(N.sp-1)]==0 | inter[,1:(N.sp-1)]==1)/length(c(inter[,1:(N.sp-1)]))
# 0.6571429


sd(res2[,i])/hyper$logsd.sd
# 0.8178026

r=ilogit(res2[,i])
hist(r,freq=F)
mean(r<1-1/Nf/N.sp)
# 0.6744067 (wider prior)


# Look at regional occupancy states:

reg.occu.index=which(substr(names(par0),1,8)=="reg.occu")

inter.active=which(inter[,1:(N.sp-1)]==1,arr.ind=TRUE)
inter.inactive=which(inter[,1:(N.sp-1)]==0,arr.ind=TRUE)


for(i in 1:length(reg.occu.index))
{
  x=res2[,reg.occu.index[i]]
  s=inter.inactive[i,2]
  f=inter.inactive[i,1]

  print.srcref(sprintf("s=%d,f=%d: m=%f, s=%f, Pr(state=0)=%f",
    s,f,mean(x),sd(x),mean(x<0)))
}

