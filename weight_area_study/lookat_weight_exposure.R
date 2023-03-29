#model="global"
#model="only_gammadyn"
#model="f"
#model="f_sf"
#model="f_sf_top"
#model="f_sf_newprior"
#model="f_sf_regocc_interact"
model="weight_exposure"

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
N.runs=10

f=list()
for(i in 1:N.runs)
{
  show(i)
  load(sprintf("%s_%02d.Rdata",model,i+N.start-1))
  if(i==1)
  {
    res=mcmc.list(mcmc(Fit3$Monitor))
  }
  if(i>1)
    res[[i]]=mcmc(Fit3$Monitor)
  f[[i]]=Fit3
}



source("make_laplace_wrapper.R")
data=make.data(par0,hyper,y)
d=Combine(f,Data=data)
d$LML

gelman.diag(res)


sel=mcmc.list()
for(i in 1:length(res))
 sel[[i]]=res[[i]][,c(2,3,4)]

sel.shell=array(0,c(length(sel),dim(sel[[1]])[1]))
sel.vol=array(0,c(length(sel),dim(sel[[1]])[1]))
sel.area=array(0,c(length(sel),dim(sel[[1]])[1]))


for(i in 1:length(res))
{
  sel.shell[i,]=sel[[i]][,1]
  sel.vol[i,]=sel[[i]][,2]
  sel.area[i,]=sel[[i]][,3]
}
for(j in 1:dim(sel.area)[2])
  for(i in 1:length(res))
{
  sel=c(ilogit(sel.shell[i,j]),ilogit(sel.vol[i,j]),ilogit(sel.area[i,j]))
  sel=sel/sum(sel)
  sel.shell[i,j]=sel[1]
  sel.vol[i,j]=sel[2]
  sel.area[i,j]=sel[3]
}


matplot(t(sel.shell))
matplot(t(sel.vol))
matplot(t(sel.area))

mean(sel.area)
#[1] 0.4260445
mean(sel.vol)
#[1] 0.281479
mean(sel.shell)
#[1] 0.2924765
