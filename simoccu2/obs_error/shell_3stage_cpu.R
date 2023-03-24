
# Requires setting the variables:
# model
# prevmodel
# N.mcmc
# sim.nr
# run.nr
# thin 
# cpus

model="shell"

source("shell_read_data.R")
source("betabin.R")
source("model_shell.R")
hyper=hyper.this
hyper$N=N
hyper$N.sp=N.sp
hyper$Nf=Nf
hyper$N.forms=N.forms
hyper$forms=forms
hyper$K1=K1
hyper$K2=K2
hyper$span=K2-K1
hyper$K=K
hyper$inter=inter


source("init_rerun.R")

if(!is.null(prevmodel) & prevmodel!="")
{
  load(sprintf("par_median_%s.RData", prevmodel))
  load(sprintf("par_mean_%s.RData", prevmodel))
  load(sprintf("par_sd_%s.RData", prevmodel))
  par0=init.params.flat2.prevmodel(10000,hyper,par.median,par.sd)
}
if(is.null(prevmodel) | prevmodel=="")
  par0=init.params.flat2(10000,hyper)
  

source("make_laplace_wrapper.R")
data=make.data(par0,hyper,y)

source("find_best_par.R")

# Stage 1:
if(cpus>1)
{
 Fit1=LaplacesDemon.hpc(Model, Data=data, par0,
  Iterations=as.integer(ceiling(N.mcmc*3/10)), 
  Algorithm="UESS", Chains=cpus, CPUs=cpus,
  Specs=list(A=Inf, B=NULL, m=100, n=0)
 )
 }
if(cpus==1)
{
 Fit1=LaplacesDemon(Model, Data=data, par0,
  Iterations=as.integer(ceiling(N.mcmc*3/10)),
  Algorithm="UESS", Specs=list(A=Inf, B=NULL, m=100, n=0) )
}
# or 
# Fit1=LaplacesDemon(Model, Data=data, par0,Iterations=6000,Algorithm="UESS", Specs=list(A=Inf, B=NULL, m=100, n=0) )

best1=find.best.par(Fit1,data)
par0=best1$par
data=make.data(par0,hyper,y)
# logprior.flat(par0,hyper)+loglik.flat(y,hyper,par0)


# Stage 2:
if(cpus>1)
{
 Fit2=LaplacesDemon.hpc(Model, Data=data, par0,
  Iterations=as.integer(ceiling(N.mcmc*3/50)),
  Algorithm="AMWG",Chains=cpus, CPUs=cpus,
  Specs=list(B=NULL, n=0, Periodicity=50)
 )
}
if(cpus==1)
{
 Fit2=LaplacesDemon(Model, Data=data, par0,
  Iterations=as.integer(ceiling(N.mcmc*3/50)),
  Algorithm="AMWG",Specs=list(B=NULL, n=0, Periodicity=50) )
}
# or
# Fit2=LaplacesDemon(Model, Data=data, par0,Iterations=2000,Algorithm="AMWG",Specs=list(B=NULL, n=0, Periodicity=50) )

best2=find.best.par(Fit2,data)
par0=best2$par
data=make.data(par0,hyper,y)


# Stage 3:
if(cpus>1)
{
 Fit3=LaplacesDemon.hpc(Model, Data=data, par0,
  Iterations=N.mcmc,Thinning=thin,
  Chains=cpus, CPUs=cpus,Specs=list(alpha.star=0.44),
  Algorithm="CHARM",Covar=Fit2$Covar)
}
if(cpus==1)
{
 Fit3=LaplacesDemon(Model, Data=data, par0,
  Iterations=N.mcmc,Thinning=thin,
  Specs=list(alpha.star=0.44),Algorithm="CHARM",Covar=Fit2$Covar)
}
# or
# Fit3=LaplacesDemon(Model, Data=data, par0,Iterations=10000,Algorithm="CHARM",Covar=Fit2$Covar)

save(Fit3, file=
  sprintf("%s_s%03d_e%02d_%03d.Rdata",model,sim.nr,err.nr,run.nr))
 
