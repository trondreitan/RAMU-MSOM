
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

source("colony_read_data.R")
source("betabin.R")
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
hyper$no.presence=no.presence
hyper$n.no.presence=n.no.presence

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
num.indep=length(res)/4
rm(res)


par0=init.params.flat(hyper)
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


reg.occu.state.index=which(substr(pname,1,14)=="reg.occu.state")
reg.occu.state=array(1,c(NN,Nf,N.sp-1))
inter.inactive=which(inter[,1:(N.sp-1)]==0,arr.ind=TRUE)
inter.active=which(inter[,1:(N.sp-1)]==1,arr.ind=TRUE)
if(length(reg.occu.state.index)>0)
 for(i in 1:NN)
{
  rbuff=array(1, c(Nf,N.sp-1))
  if(length(grep("dropint",model))>0)
    rbuff[inter.inactive]=
    1*(res2[i,reg.occu.state.index]>0)      
  if(length(grep("dropint",model))==0)
    rbuff=matrix(1*(res2[i,reg.occu.state.index]>0), ncol=N.sp-1)
  rbuff[inter.active]=1
  reg.occu.state[i,,]=rbuff
}



d=dim(reg.occu.state)
reg.occu.state2=array(NA,c(d[1],d[2],d[3]+1))
reg.occu.state2[,,1:d[3]]=reg.occu.state
reg.occu.state2[,,d[3]+1]=1


# Sanity check this:
counts=as.matrix(y[,which(substr(names(y),1,7)=="Species")])
# Check if formation absence is set for a species-formation combination with
# nonzero counts:
m.max=0
f.max=0
s.max=0
b=0
for(f in 1:Nf)
 for(s in 1:(N.sp-1))
  #if(sum(counts[y$form.nr==f,s])>0)
{
  m=mean(reg.occu.state[,f,s]==0)

  #if(m>0)
  #  print.srcref(sprintf("Formation=%s(%d) species=%s(%d): counts=%d reg.absence=%4.1f%%",
  #      y$Formation_name[y$form.nr==f][1],f,sp[s],s, sum(counts[y$form.nr==f,s]), m*100))

  if(m>0.3)
    b=b+1

  if(m>m.max)
  {
    m.max=m
    s.max=s
    f.max=f
  }
}

#print.srcref(sprintf("Maximum absence probability: %s - %s: %f", sp[s.max],
#   y$Formation_name[y$form.nr==f.max][1],m.max))
   



#o18.formation=rep(0,Nf)
#ca.formation=rep(0,Nf)
#for(f in 1:Nf)
#{
#  o18.formation[f]=y$O18.mean[y$form.nr==f][1]
#  ca.formation[f]=y$ca.mean[y$form.nr==f][1]	
#}
# o18 and ca not in yet


logit.occu.formation.species=array(NA,c(NN,Nf,N.sp-1))
for(i in 1:Nf)
 logit.occu.formation.species[,i,]=const.logit.occu+f.occu[,i]%*%t(rep(1,N.sp-1))+sf.occu[,i,]
occu.formation.species=ilogit(logit.occu.formation.species)


d=dim(occu.formation.species)
occu.formation.species2=array(NA,c(d[1],d[2],d[3]+1))
occu.formation.species2[,,1:d[3]]=occu.formation.species
occu.formation.species2[,,d[3]+1]=1


#o18bin=array(0,c(NN,Nf,N.sp))
#if(is.o18bin)
#{
#  for(f in 1:Nf)
#   for(s in 1:N.sp)
#   {
#     o18.o.ind=which(names(par0)==sprintf("o18.optim.bin%d",s))
#     o18.w.ind=which(names(par0)==sprintf("o18.lwidth.bin%d",s))
#
#     o18bin[,f,s]=-(o18.formation[f]-res2[,o18.o.ind])^2/exp(res2[,o18.w.ind])^2
#   }
#}
#
#cabin=array(0,c(NN,Nf,N.sp))
#if(is.cabin)
#{
#  for(f in 1:Nf)
#   for(s in 1:N.sp)
#   {
#     ca.o.ind=which(names(par0)==sprintf("ca.optim.bin%d",s))
#     ca.w.ind=which(names(par0)==sprintf("ca.lwidth.bin%d",s))
#
#     cabin[,f,s]=-(ca.formation[f]-res2[,ca.o.ind])^2/exp(res2[,ca.w.ind])^2
#   }
#}

log.lambda.formation.species=array(NA,c(NN,Nf,N.sp))
for(i in 1:Nf)
 log.lambda.formation.species[,i,]=
    const.log.lambda+
    f.lambda[,i]%*%t(rep(1,N.sp))+
    sf.lambda[,i,] # +o18bin[,i,]+cabin[,i,]


#for(f in 1:Nf)
# for(s in 1:N.sp)
#  print.srcref(sprintf("%s-%d: %f", 
#      sp[s],f,mean(log.lambda.formation.species[,f,s])))

lambda.formation.species=exp(log.lambda.formation.species)


#for(f in 1:Nf)
# for(s in 1:N.sp)
#  print.srcref(sprintf("%s-%d: %f %f", 
#      sp[s],f,mean(lambda.formation.species[,f,s]),
#       median(lambda.formation.species[,f,s])))


bin.formation.species=1-exp(-lambda.formation.species)


lower=function(x) quantile(x, 0.025)
upper=function(x) quantile(x, 0.975)

forms=sort(unique(y$form.nr))
Nf2=length(forms) # check that this is the same as Nf
Nf3=max(forms) # check that this is also the same as Nf

# Formation-age info not yet in
#K1=rep(NA,Nf)
#K2=rep(NA,Nf)
#for(i in 1:Nf)
#{
#  K1[i]=y$time.start[y$form.nr==i][1]
#  K2[i]=y$time.end[y$form.nr==i][1]
#}
#K.mid=(K1+K2)/2
#
#K=K.mid
#Ks=K1
#Ke=K2

K=1:Nf



# Now make relative abundance, formation by formation:
# abundance.given.occupancy = 
#    bin.formation.species/(1-bin.formation.species)
#abundance.given.occupancy=-log(1-bin.formation.species)


log.lambda.corrected=array(NA,c(NN,Nf,N.sp))
for(i in 1:Nf)
 log.lambda.corrected[,i,]=
    const.log.lambda+
    # f.lambda[,i]%*%t(rep(1,N.sp))+
    sf.lambda[,i,] # +o18bin[,i,]+cabin[,i,]

#for(f in 1:Nf)
# for(s in 1:N.sp)
#  print.srcref(sprintf("%s-%d: %f %f", 
#      sp[s],f,mean(log.lambda.corrected[,f,s]), mean(log.lambda.formation.species[,f,s])))

lambda.corrected=exp(log.lambda.corrected)

abundance.given.occupancy=lambda.corrected
rm(lambda.corrected)
rm(log.lambda.corrected)
rm(lambda.formation.species)

#for(f in 1:Nf)
# for(s in 1:N.sp)
#  print.srcref(sprintf("%s-%d: %f", 
#      sp[s],f,mean(abundance.given.occupancy[,f,s])))




abundance=abundance.given.occupancy*occu.formation.species2*reg.occu.state2
abundance.corrected=abundance

gam.index=which(substr(names(par0),1,11)=="logit.gamma")
gam.raw=ilogit(res2[,gam.index])
gam=array(NA,c(dim(res2)[1],Nf,N.un))
if(N.un>0)
 for(i in 1:N.un)
  gam[,,i]=gam.raw[,(i-1)*Nf+1:Nf]

for(i in 1:length(unid.genus.nr))
 for(s in which(genus.nr==unid.genus.nr[i]))
  for(f in 1:Nf)
    abundance[,f,s]=abundance.corrected[,f,s]*gam[,f,i]


#for(f in 1:Nf)
# for(s in 1:N.sp)
#  print.srcref(sprintf("%s-%d: %f %f", 
#      sp[s],f,mean(abundance.corrected[,f,s]),
#      mean(abundance[,f,s])))


rel.abundance.raw=abundance
for(i in 1:Nf)
 rel.abundance.raw[,i,]=abundance[,i,]/
         apply(abundance[,i,],1,sum)

rel.abundance=abundance.corrected
for(i in 1:Nf)
 rel.abundance[,i,]=abundance.corrected[,i,]/
         apply(abundance.corrected[,i,],1,sum)



mean.rel.abundance=array(NA,c(Nf,N.sp))
median.rel.abundance=array(NA,c(Nf,N.sp))
lower.rel.abundance=array(NA,c(Nf,N.sp))
upper.rel.abundance=array(NA,c(Nf,N.sp))
mean.rel.abundance.raw=array(NA,c(Nf,N.sp))
median.rel.abundance.raw=array(NA,c(Nf,N.sp))
lower.rel.abundance.raw=array(NA,c(Nf,N.sp))
upper.rel.abundance.raw=array(NA,c(Nf,N.sp))

for(f in 1:Nf)
  for(s in 1:N.sp)
  {
    mean.rel.abundance[f,s]=mean(rel.abundance[,f,s])
    median.rel.abundance[f,s]=median(rel.abundance[,f,s])
    lower.rel.abundance[f,s]=lower(rel.abundance[,f,s])
    upper.rel.abundance[f,s]=upper(rel.abundance[,f,s])
    mean.rel.abundance.raw[f,s]=mean(rel.abundance.raw[,f,s])
    median.rel.abundance.raw[f,s]=
      median(rel.abundance.raw[,f,s])
   lower.rel.abundance.raw[f,s]=lower(rel.abundance.raw[,f,s])
   upper.rel.abundance.raw[f,s]=upper(rel.abundance.raw[,f,s])
  }
#save(mean.rel.abundance, file=sprintf("mean_ra_%s.Rdata",model))
#save(median.rel.abundance, file=sprintf("median_ra_%s.Rdata",model))
#save(lower.rel.abundance, file=sprintf("lower_ra_%s.Rdata",model))
#save(upper.rel.abundance, file=sprintf("upper_ra_%s.Rdata",model))


raw.rel.abundance=0*mean.rel.abundance
sp.index=which(substr(names(y),1,8)=="Species_")
for(f in 1:Nf)
 for(s in 1:N.sp)
  {
    raw.rel.abundance[f,s]=sum(y[y$form.nr==f,sp.index[s]])/sum(y[y$form.nr==f,sp.index])
  }  


rel.abundance.colony.mod[sim,err,,]=mean.rel.abundance
rel.abundance.colony.raw[sim,err,,]=raw.rel.abundance
rel.abundance.colony.upper[sim,err,,]=upper.rel.abundance
rel.abundance.colony.lower[sim,err,,]=lower.rel.abundance



cols=rainbow(N.sp)

#plot(rnorm(1000,1,0.1),rnorm(1000,0,0.1),
#  col=cols[1],xlim=c(0,12))
#for(s in 2:N.sp)
# points(rnorm(1000,s,0.1),rnorm(1000,0,0.1),col=cols[s])
#sp

png(sprintf("%s_rel_abundance_s%03d_e%02d.png",model,sim.nr,err.nr),height=1200,width=1600)
par(cex=3,cex.lab=2.9,cex.sub=1.9,cex.axis=2.2,cex.main=1.9,mar=c(10,5,6,1))
par(mfrow=c(1,1))
plot(K, mean.rel.abundance[,1], type="b", ylim=c(0.0001,1),
  xlab="time", ylab="Rel. Abundance",col=cols[1],lwd=4,log="y")
for(s in 2:N.sp)
  lines(K, mean.rel.abundance[,s], type="b",col=cols[s],lwd=4)
#for(s in 1:N.sp)
#  lines(K, mean.rel.abundance.raw[,s], type="b",col=cols[s],lwd=1,lty=2)
for(s in 1:N.sp)
   lines(K, rel.abundance.real[,s], type="b",col=cols[s],lwd=1,lty=1)
for(s in 1:N.sp)
   lines(K, raw.rel.abundance[,s], type="b",col=cols[s],lwd=2,lty=3)

dev.off()




rmse.mod[sim,err]=sqrt(mean((mean.rel.abundance-rel.abundance.real)^2))
rmse.raw[sim,err]=sqrt(mean((raw.rel.abundance-rel.abundance.real)^2))

c(rmse.mod[sim,err], rmse.raw[sim,err])
# 0.01723286 0.02190014

include=which(raw.rel.abundance>0,arr.ind=T)
log.rmse.mod[sim,err]=sqrt(mean((log(mean.rel.abundance[include])-log(rel.abundance.real[include]))^2))
log.rmse.raw[sim,err]=sqrt(mean((log(raw.rel.abundance[include])-log(rel.abundance.real[include]))^2))
c(log.rmse.mod[sim,err], log.rmse.raw[sim,err])
# 0.3766723 0.5787378

}


