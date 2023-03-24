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


runs=1:4
sim.nrs=1:10
err.nrs=1:9

p.err=c(0,0.0001,0.0003,0.001,0.003,0.01,0.03,0.1,0.3)
N.err=length(p.err)


n.sim=length(sim.nrs)
n.err=length(err.nrs)

model="shell"



rmse.mod=array(0,c(n.sim,n.err))
rmse.raw=array(0,c(n.sim,n.err))
log.rmse.mod=array(0,c(n.sim,n.err))
log.rmse.raw=array(0,c(n.sim,n.err))


rel.abundance.shell.mod=array(0,c(n.sim,n.err,Nf,N.sp))
rel.abundance.shell.raw=array(0,c(n.sim,n.err,Nf,N.sp))
rel.abundance.shell.upper=array(0,c(n.sim,n.err,Nf,N.sp))
rel.abundance.shell.lower=array(0,c(n.sim,n.err,Nf,N.sp))

source("read_rel_abundance_shell.R")


rmse.mod.mean=apply(rmse.mod,2,mean)
rmse.raw.mean=apply(rmse.raw,2,mean)
log.rmse.mod.mean=apply(log.rmse.mod,2,mean)
log.rmse.raw.mean=apply(log.rmse.raw,2,mean)


for(err in 1:length(err.nrs))
  print.srcref(sprintf("Error %02d (p=%7.5f): rmse.mod=%f rmse.raw=%f",
    err.nrs[err],p.err[err.nrs[err]],
    rmse.mod.mean[err], rmse.raw.mean[err]))



png("RMSE_compare_shell.png",height=2000,width=3000)
par(cex=4)
plot(p.err,rmse.raw.mean,type="b",log="xy",
  xlab="Measurement error probability",ylab="RMSE",
  ylim=c(min(min(rmse.raw.mean),min(rmse.mod.mean)),
	 max(max(rmse.raw.mean),max(rmse.mod.mean))),
	 col="red")
lines(p.err,rmse.mod.mean,type="b",col="black")
dev.off()

png("logRMSE_compare_shell.png",height=2000,width=3000)
par(cex=4)
plot(p.err,log.rmse.raw.mean,type="b",col="red",log="xy",xlab="Measurement error probability",ylab="log-RMSE")
lines(p.err,log.rmse.mod.mean,type="b",col="black")
dev.off()

rmse.mod.shell=rmse.mod.mean
rmse.raw.shell=rmse.raw.mean
save(rmse.mod.shell, file="rmse_mod_shell.Rdata")
save(rmse.raw.shell, file="rmse_raw_shell.Rdata")

logrmse.mod.shell=log.rmse.mod.mean
logrmse.raw.shell=log.rmse.raw.mean
save(logrmse.mod.shell, file="logrmse_mod_shell.Rdata")
save(logrmse.raw.shell, file="logrmse_raw_shell.Rdata")



save(rel.abundance.shell.mod, file="rel_abundance_shell_mod.Rdata")
save(rel.abundance.shell.raw, file="rel_abundance_shell_raw.Rdata")
save(rel.abundance.shell.upper, file="rel_abundance_shell_upper.Rdata")
save(rel.abundance.shell.lower, file="rel_abundance_shell_lower.Rdata")

save(rel.abundance.real,  file="rel_abundance_real.Rdata")
