
for(sim in 1:length(sim.nrs))
 for(err in 1:length(err.nrs))
{

  sim.nr=sim.nrs[sim]
  err.nr=err.nrs[err]

is.o18bin=FALSE
if(length(grep("o18bin",model))>0)
  is.o18bin=TRUE

is.cabin=FALSE
if(length(grep("cabin",model))>0)
  is.cabin=TRUE

source(sprintf("%s_read_data.R",model))
source("betabin.R")
source(sprintf("model_%s.R",model))
hyper=hyper.this
hyper$N=N
hyper$N.sp=N.sp
hyper$Nf=Nf
hyper$N.forms=N.forms

library(coda)




# Collate (if multiple runs)
N.start=runs[1]
N.runs=length(runs)

ll=rep(NA,N.runs)
# all
f=list()
for(i in 1:N.runs)
{
  load(sprintf("%s_s%03d_e%02d_%03d.Rdata",model,sim.nr,err.nr,i+N.start-1))
  if(i==1)
  {
    res=mcmc.list(mcmc(Fit3$Monitor))
    f=list(f1=Fit3)
  }
  if(i>1)
  {
     res[[i]]=mcmc(Fit3$Monitor)
     f[[i]]=Fit3
  }
  ll[i]=Fit3$LML
}

source("make_laplace_wrapper.R")
par0=init.params.flat(hyper)
data=make.data(par0,hyper,y)


n1=dim(res[[1]])[2]
n2=0
for(i in 1:length(res))
{
  n2=n2+dim(res[[i]])[1]
}
res2=array(NA,c(n2,n1))
j=1
for(i in 1:length(res))
{
  nn=dim(res[[i]])[1]
  res2[j:(j+nn-1),]=res[[i]]
  j=j+nn
}

NN=dim(res2)[1]

pname=names(par0)

const.occu.index=which(substr(pname,1,16)=="logit.const.occu")
const.logit.occu=res2[,const.occu.index]

const.lambda.index=which(substr(pname,1,16)=="log.const.lambda")
const.log.lambda=res2[,const.lambda.index]

f.occu.index=which(substr(pname,1,6)=="f.occu")
f.occu=res2[,f.occu.index]

f.lambda.index=which(substr(pname,1,8)=="f.lambda")
f.lambda=res2[,f.lambda.index]

sf.occu.index=which(substr(pname,1,7)=="sf.occu")
sf.occu=array(NA,c(NN,Nf,N.sp-1))
for(i in 1:NN)
  sf.occu[i,,]=matrix(res2[i,sf.occu.index], ncol=N.sp-1)

sf.lambda.index=which(substr(pname,1,9)=="sf.lambda")
sf.lambda=array(NA,c(NN,Nf,N.sp))
for(i in 1:NN)
  sf.lambda[i,,]=matrix(res2[i,sf.lambda.index], ncol=N.sp)



logit.occu.formation.species=array(NA,c(NN,Nf,N.sp-1))
for(i in 1:Nf)
 logit.occu.formation.species[,i,]=const.logit.occu+f.occu[,i]%*%t(rep(1,N.sp-1))+sf.occu[,i,]
occu.formation.species=ilogit(logit.occu.formation.species)


d=dim(occu.formation.species)
occu.formation.species2=array(NA,c(d[1],d[2],d[3]+1))
occu.formation.species2[,,1:d[3]]=occu.formation.species
occu.formation.species2[,,d[3]+1]=1



log.lambda.formation.species=array(NA,c(NN,Nf,N.sp))
for(i in 1:Nf)
 log.lambda.formation.species[,i,]=const.log.lambda+f.lambda[,i]%*%t(rep(1,N.sp))+sf.lambda[,i,]
lambda.formation.species=exp(log.lambda.formation.species)

bin.formation.species=1-exp(-lambda.formation.species)


lower=function(x) quantile(x, 0.025)
upper=function(x) quantile(x, 0.975)

forms=sort(unique(y$form.nr))
Nf2=length(forms) # check that this is the same as Nf
Nf3=max(forms) # check that this is also the same as Nf

#K1=rep(NA,Nf)
#K2=rep(NA,Nf)
#for(i in 1:Nf)
#{
#  K1[i]=y$time.start[y$form.nr==i][1]
#  K2[i]=y$time.end[y$form.nr==i][1]
#}



#K.mid=(K1+K2)/2

#K=K.mid
#Ks=K1
#Ke=K2

K=1:Nf

 
sp=names(y)[substr(names(y),1,8)=="Species_"]
#sp=substr(sp,9,nchar(sp))




#####################
# Summary plots:
#####################


# Plot raw detection rate data:

raw.rel.abundance=array(NA,c(Nf,N.sp))
for(s in 1:N.sp)
 for(f in 1:Nf)
   raw.rel.abundance[f,s]=mean(y[y$form.nr==f,names(y)==sp[s]]/y$Total[y$form.nr==f])
  



# Now make relative abundance, formation by formation:
# abundance.given.occupancy=bin.formation.species/(1-bin.formation.species)
#abundance.given.occupancy=-log(1-bin.formation.species)
abundance.given.occupancy=lambda.formation.species
abundance=abundance.given.occupancy*occu.formation.species2


rel.abundance=abundance
for(i in 1:Nf)
  rel.abundance[,i,]=abundance[,i,]/apply(abundance[,i,],1,sum)

mean.rel.abundance=array(NA,c(Nf,N.sp))
median.rel.abundance=array(NA,c(Nf,N.sp))
lower.rel.abundance=array(NA,c(Nf,N.sp))
upper.rel.abundance=array(NA,c(Nf,N.sp))
for(f in 1:Nf)
  for(s in 1:N.sp)
  {
    mean.rel.abundance[f,s]=mean(rel.abundance[,f,s])
    median.rel.abundance[f,s]=median(rel.abundance[,f,s])
    lower.rel.abundance[f,s]=lower(rel.abundance[,f,s])
    upper.rel.abundance[f,s]=upper(rel.abundance[,f,s])
  }


rel.abundance.shell.mod[sim,err,,]=mean.rel.abundance
rel.abundance.shell.raw[sim,err,,]=raw.rel.abundance
rel.abundance.shell.upper[sim,err,,]=upper.rel.abundance
rel.abundance.shell.lower[sim,err,,]=lower.rel.abundance


png(sprintf("%s_rel_abundance_s%03d_e%02d.png",model,sim.nr,err.nr),height=1200,width=1600)
par(cex=3,cex.lab=2.9,cex.sub=1.9,cex.axis=2.2,cex.main=1.9,mar=c(10,5,6,1))
cols=rainbow(N.sp)
par(mfrow=c(1,1))
plot(K, mean.rel.abundance[,1], type="b", ylim=c(0.001,1),
  xlab="time", ylab="Rel. Abundance",col=cols[1],lwd=4,log="y")
for(s in 2:N.sp)
  lines(K, mean.rel.abundance[,s], type="b",col=cols[s],lwd=4)
for(s in 1:N.sp)
   lines(K, rel.abundance.real[,s], type="b",col=cols[s],lwd=1,lty=1)
for(s in 1:N.sp)
   lines(K, raw.rel.abundance[,s], type="b",col=cols[s],lwd=2,lty=3)
dev.off()


rmse.mod[sim,err]=sqrt(mean((mean.rel.abundance-rel.abundance.real)^2))
rmse.raw[sim,err]=sqrt(mean((raw.rel.abundance-rel.abundance.real)^2))

c(rmse.mod[sim,err], rmse.raw[sim,err])

include=which(raw.rel.abundance>0,arr.ind=T)
log.rmse.mod[sim,err]=sqrt(mean((log(mean.rel.abundance[include])-log(rel.abundance.real[include]))^2))
log.rmse.raw[sim,err]=sqrt(mean((log(raw.rel.abundance[include])-log(rel.abundance.real[include]))^2))
c(log.rmse.mod[sim,err], log.rmse.raw[sim,err])

}

