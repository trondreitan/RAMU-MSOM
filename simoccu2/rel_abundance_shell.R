##############################
# Simulation parameters:
##############################

# From simoccu_save.R: 

logit=function(x) log(x/(1-x))
ilogit=function(x) 1/(1+exp(-x))

# Important settings:
Nf=10
N.focal=12
N.sp=N.focal+1
sp=c(sprintf("sp%02d",1:N.sp),"Super")

N.site=10
N.shell=100



# Constant variants:
lambda.high=0.3
lambda.mod=0.1
lambda.low=0.02
psi.high=0.8
psi.mod=0.6
psi.low=0.4


# Set species constants:
lambda.sp=rep(0,N.sp)
psi.sp=rep(0,N.sp)

lambdas=c(lambda.high,lambda.mod,lambda.low)
psis=c(psi.low,psi.high, psi.mod)
for(i in 1:3)
{
  lambda.sp[(i-1)*4+1:4]=lambdas[i]
  psi.sp[(i-1)*4+1:4]=psis[i]
}

# Superspecies:
lambda.sp[N.sp]=lambda.high
psi.sp[N.sp]=1



# Dynamics:
lambda=array(0,c(Nf,N.sp))
psi=array(0,c(Nf,N.sp))

lambda.maxvar=log(3)
occu.maxvar=1.5
for(s in 1:N.sp)
{
 for(f in 1:Nf)
 {
   if(s%%4==1)
   {
     lambda[f,s]=
       exp(log(lambda.sp[s])+
               (f-Nf/2-.5)/(Nf/2-.5)*lambda.maxvar)
     if(s<N.sp)
       psi[f,s]=psi.sp[s]
   }

   if(s%%4==2)
   {
     lambda[f,s]=lambda.sp[s]
     if(s<N.sp)
       psi[f,s]=ilogit(logit(psi.sp[s])+
                       (f-Nf/2-.5)/(Nf/2-.5)*occu.maxvar)
   }

  if(s%%4==3)
   {
     lambda[f,s]=
       exp(log(lambda.sp[s])-
               (f-Nf/2-0.5)/(Nf/2-0.5)*lambda.maxvar)
     if(s<N.sp)
       psi[f,s]=ilogit(logit(psi.sp[s])-
                       (f-Nf/2-.5)/(Nf/2-.5)*occu.maxvar)
   }

   if(s%%4==0)
   {
     lambda[f,s]=
       exp(log(lambda.sp[s])-
               (f-Nf/2-.5)/(Nf/2-.5)*lambda.maxvar)
     if(s<N.sp)
       psi[f,s]=ilogit(logit(psi.sp[s])+
                           (f-Nf/2-.5)/(Nf/2-.5)*occu.maxvar)
   }

   if(s==N.sp)
   {
     lambda[f,s]=exp(lambda.maxvar)*lambda.high
     psi[f,s]=1
   }
 } 
}

psilambda=psi*lambda

rel.abundance.real=psilambda
for(i in 1:Nf)
  rel.abundance.real[i,]=rel.abundance.real[i,]/sum(rel.abundance.real[i,])






##############################
# Analysis results:
##############################


model="shell"

runs=5:8


rmse.mod=rep(0,100)
rmse.raw=rep(0,100)
log.rmse.mod=rep(0,100)
log.rmse.raw=rep(0,100)


rel.abundance.shell.mod=array(0,c(100,Nf,N.sp))
rel.abundance.shell.raw=array(0,c(100,Nf,N.sp))
rel.abundance.shell.upper=array(0,c(100,Nf,N.sp))
rel.abundance.shell.lower=array(0,c(100,Nf,N.sp))

for(sim.nr in 1:100)
{

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
  load(sprintf("%s_s%03d_%03d.Rdata",model,sim.nr,i+N.start-1))
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


rel.abundance.shell.mod[sim.nr,,]=mean.rel.abundance
rel.abundance.shell.raw[sim.nr,,]=raw.rel.abundance
rel.abundance.shell.upper[sim.nr,,]=upper.rel.abundance
rel.abundance.shell.lower[sim.nr,,]=lower.rel.abundance


png(sprintf("%s_rel_abundance_s%03d.png",model,sim.nr),height=1200,width=1600)
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


rmse.mod[sim.nr]=sqrt(mean((mean.rel.abundance-rel.abundance.real)^2))
rmse.raw[sim.nr]=sqrt(mean((raw.rel.abundance-rel.abundance.real)^2))

c(rmse.mod[sim.nr], rmse.raw[sim.nr])

include=which(raw.rel.abundance>0,arr.ind=T)
log.rmse.mod[sim.nr]=sqrt(mean((log(mean.rel.abundance[include])-log(rel.abundance.real[include]))^2))
log.rmse.raw[sim.nr]=sqrt(mean((log(raw.rel.abundance[include])-log(rel.abundance.real[include]))^2))
c(log.rmse.mod[sim.nr], log.rmse.raw[sim.nr])

}


c(mean(rmse.mod),mean(rmse.raw), 1-mean(rmse.mod)/mean(rmse.raw))
# 0.03488984 0.08465984 0.58788203
# 0.03490922 0.08465984 0.58765308

c(mean(log.rmse.mod),mean(log.rmse.raw),1-mean(log.rmse.mod)/mean(log.rmse.raw))
# 0.5458062 1.0787304 0.4940291
# 0.5502492 1.0787304 0.4899104






save(rel.abundance.shell.mod, file="rel_abundance_shell_mod.Rdata")
save(rel.abundance.shell.raw, file="rel_abundance_shell_raw.Rdata")
save(rel.abundance.shell.upper, file="rel_abundance_shell_upper.Rdata")
save(rel.abundance.shell.lower, file="rel_abundance_shell_lower.Rdata")

save(rel.abundance.real,  file="rel_abundance_real.Rdata")
