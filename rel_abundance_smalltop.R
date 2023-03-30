model="smalltop"

source("read_data.R")
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

res2=mcmc(d$Monitor)

NN=dim(res2)[1]


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
reg.occu.state=array(1,c(NN,Nf,N.sp))

inter.active=which(inter==1,arr.ind=TRUE)
inter.inactive=which(inter==0,arr.ind=TRUE)
  
if(length(reg.occu.state.index)>0)
 for(i in 1:NN)
{
  rbuff=array(1, c(Nf,N.sp))
  rbuff[inter.inactive]=
    1*(res2[i,reg.occu.state.index]>0)      
  rbuff[inter.active]=1
  reg.occu.state[i,,]=rbuff
}
 

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

print.srcref(sprintf("Maximum absence probability: %s - %s: %f", sp[s.max],
  y$Formation_name[y$form.nr==f.max][1],m.max))
# wider prior:  Maximum absence probability: Calloporina_angustipora - Lower Kai-iwi Shellbed: 0.508407


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
 log.lambda.formation.species[,i,]=
    const.log.lambda+
    f.lambda[,i]%*%t(rep(1,N.sp))+
    sf.lambda[,i,] # +o18bin[,i,]+cabin[,i,]


for(f in 1:Nf)
 for(s in 1:N.sp)
  print.srcref(sprintf("%s-%d: %f", 
      sp[s],f,mean(log.lambda.formation.species[,f,s])))

lambda.formation.species=exp(log.lambda.formation.species)


lower=function(x) quantile(x, 0.025)
upper=function(x) quantile(x, 0.975)

forms=sort(unique(y$form.nr))
Nf2=length(forms) # check that this is the same as Nf
Nf3=max(forms) # check that this is also the same as Nf

K=1:Nf



# Now make relative abundance, formation by formation:

log.lambda.corrected=array(NA,c(NN,Nf,N.sp))
for(i in 1:Nf)
 log.lambda.corrected[,i,]=
    const.log.lambda+
    # f.lambda[,i]%*%t(rep(1,N.sp))+
    sf.lambda[,i,] 

lambda.corrected=exp(log.lambda.corrected)

abundance.given.occupancy=lambda.corrected
rm(lambda.corrected)
rm(log.lambda.corrected)
rm(lambda.formation.species)

for(f in 1:Nf)
 for(s in 1:N.sp)
  print.srcref(sprintf("%s-%d: %f", 
      sp[s],f,mean(abundance.given.occupancy[,f,s])))

abundance=abundance.given.occupancy*
   occu.formation.species2*reg.occu.state
abundance.corrected=abundance


gam.index=which(substr(names(par0),1,11)=="logit.gamma")
gam.from.mcmc=ilogit(res2[,gam.index])
gam=array(NA,c(dim(res2)[1],Nf,N.un))
gam.mean=array(NA,c(Nf,N.un))
if(N.un>0)
 for(i in 1:N.un)
 {
  gam[,,i]=gam.from.mcmc[,(i-1)*Nf+1:Nf]
  for(f in 1:Nf)
   gam.mean[f,i]=mean(gam[,f,i])
 }


lambda.mean=array(NA,c(Nf,N.sp))
psi.mean=array(NA,c(Nf,N.sp))
for(f in 1:Nf)
 for(s in 1:N.sp)
  {
   lambda.mean[f,s]=
        mean(abundance.given.occupancy[,f,s]*reg.occu.state[,f,s])
	
   psi.mean[f,s]=
       mean(occu.formation.species2[,f,s]*reg.occu.state[,f,s])
  }

hist(log(lambda.mean),freq=F)
xx=seq(-10,10,0.01)
lines(xx,dnorm(xx,hyper$log.lambda.mu, hyper$log.lambda.sd))

hist(logit(psi.mean[,1:(N.sp-1)]),freq=F)
xx=seq(-10,10,0.01)
lines(xx,dnorm(xx,hyper$logit.mu, hyper$logit.sd))


rel.abundance=abundance.corrected
for(i in 1:Nf)
 rel.abundance[,i,]=abundance.corrected[,i,]/
         apply(abundance.corrected[,i,],1,sum)


mean.rel.abundance=array(NA,c(Nf,N.sp))
median.rel.abundance=array(NA,c(Nf,N.sp))
lower.rel.abundance=array(NA,c(Nf,N.sp))
upper.rel.abundance=array(NA,c(Nf,N.sp))

for(f in 1:Nf)
{
  for(s in 1:N.sp)
  {
    mean.rel.abundance[f,s]=mean(rel.abundance[,f,s])
    median.rel.abundance[f,s]=median(rel.abundance[,f,s])
    lower.rel.abundance[f,s]=lower(rel.abundance[,f,s])
    upper.rel.abundance[f,s]=upper(rel.abundance[,f,s])
  }
}

save(mean.rel.abundance, file=sprintf("mean_ra_%s.Rdata",model))
save(median.rel.abundance, file=sprintf("median_ra_%s.Rdata",model))
save(lower.rel.abundance, file=sprintf("lower_ra_%s.Rdata",model))
save(upper.rel.abundance, file=sprintf("upper_ra_%s.Rdata",model))

raw.rel.abundance.uncorr=0*mean.rel.abundance
raw.rel.abundance=0*mean.rel.abundance
sp.index=which(substr(names(y),1,8)=="Species_")
for(f in 1:Nf)
 for(s in 1:N.sp)
  {
    raw.rel.abundance.uncorr[f,s]=
      sum(y[y$form.nr==f,sp.index[s]])/sum(y[y$form.nr==f,sp.index])
    if(s<N.sp)
      raw.gam=(0.5+sum(y[y$form.nr==f,sp.index[genus.nr==genus.nr[s]]]))/
	(1+sum(y[y$form.nr==f,names(y)==sprintf("Genus_%s",genus[s])]))
    if(s==N.sp)
      raw.gam=1
    raw.rel.abundance[f,s]=
      sum(y[y$form.nr==f,sp.index[s]])/sum(y[y$form.nr==f,sp.index])/raw.gam
  }  

save(raw.rel.abundance, file="raw_rel_abundance.Rdata")
save(raw.rel.abundance.uncorr, file="raw_rel_abundance_uncorr.Rdata")



#PS: abundance=abundance.given.occupancy*occu.formation.species2*reg.occu.state
mean.abundance.given.occupancy=array(NA,c(Nf,N.sp))
median.abundance.given.occupancy=array(NA,c(Nf,N.sp))
lower.abundance.given.occupancy=array(NA,c(Nf,N.sp))
upper.abundance.given.occupancy=array(NA,c(Nf,N.sp))
mean.site.occupancy=array(NA,c(Nf,N.sp))
median.site.occupancy=array(NA,c(Nf,N.sp))
lower.site.occupancy=array(NA,c(Nf,N.sp))
upper.site.occupancy=array(NA,c(Nf,N.sp))
mean.reg.occu=array(NA,c(Nf,N.sp))

mean.occupancy=array(NA,c(Nf,N.sp))
median.occupancy=array(NA,c(Nf,N.sp))
lower.occupancy=array(NA,c(Nf,N.sp))
upper.occupancy=array(NA,c(Nf,N.sp))

for(f in 1:Nf)
{
  for(s in 1:N.sp)
  {
    mean.abundance.given.occupancy[f,s]=mean(abundance.given.occupancy[,f,s]*reg.occu.state[,f,s])
    median.abundance.given.occupancy[f,s]=median(abundance.given.occupancy[,f,s]*reg.occu.state[,f,s])
    lower.abundance.given.occupancy[f,s]=lower(abundance.given.occupancy[,f,s]*reg.occu.state[,f,s])
    upper.abundance.given.occupancy[f,s]=upper(abundance.given.occupancy[,f,s]*reg.occu.state[,f,s])
    mean.site.occupancy[f,s]=mean(occu.formation.species2[,f,s])
    median.site.occupancy[f,s]=median(occu.formation.species2[,f,s])
    lower.site.occupancy[f,s]=lower(occu.formation.species2[,f,s])
    upper.site.occupancy[f,s]=upper(occu.formation.species2[,f,s])
    mean.reg.occu[f,s]=mean(reg.occu.state[,f,s])

    mean.occupancy[f,s]=mean(reg.occu.state[,f,s]*occu.formation.species2[,f,s])
    median.occupancy[f,s]=median(reg.occu.state[,f,s]*occu.formation.species2[,f,s])
    lower.occupancy[f,s]=lower(reg.occu.state[,f,s]*occu.formation.species2[,f,s])
    upper.occupancy[f,s]=upper(reg.occu.state[,f,s]*occu.formation.species2[,f,s])
  }
}
save(mean.abundance.given.occupancy, file="mean_abundance_given_occupancy.Rdata");
save(median.abundance.given.occupancy, file="median_abundance_given_occupancy.Rdata");
save(lower.abundance.given.occupancy, file="lower_abundance_given_occupancy.Rdata");
save(upper.abundance.given.occupancy, file="upper_abundance_given_occupancy.Rdata");
save(mean.site.occupancy, file="mean_site_occupancy.Rdata");
save(median.site.occupancy, file="median_site_occupancy.Rdata");
save(lower.site.occupancy, file="lower_site_occupancy.Rdata");
save(upper.site.occupancy, file="upper_site_occupancy.Rdata");
save(mean.occupancy, file="mean_occupancy.Rdata");
save(median.occupancy, file="median_occupancy.Rdata");
save(lower.occupancy, file="lower_occupancy.Rdata");
save(upper.occupancy, file="upper_occupancy.Rdata");
save(mean.reg.occu, file="mean_reg_occu.Rdata")






# Preliminary plots:

for(s in 1:N.sp)
{
  png(sprintf("rel_abundance_%s_s%02d.png",model,s), width=1600,height=1200)
  par(cex=2)
  plot(1:Nf,mean.rel.abundance[,s],xlab="Formations",
	ylab="Relative abundance", ylim=c(min(lower.rel.abundance[,s]),
	max(upper.rel.abundance[,s])),type="b",main=sp[s])

  lines(1:Nf, lower.rel.abundance[,s], col="red")
  lines(1:Nf, upper.rel.abundance[,s], col="red")
  dev.off()
}



