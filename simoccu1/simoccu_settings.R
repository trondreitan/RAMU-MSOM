

logit=function(x) log(x/(1-x))
ilogit=function(x) 1/(1+exp(-x))

# Important settings:
Nf=10

N.genera=6
genera.names=sprintf("G%02d",1:N.genera)

N.species.per.genera=4
genera.num.species=rep(N.species.per.genera,N.genera)

N.focal=sum(genera.num.species)
N.sp=N.focal+1

species.genera=rep(0,N.sp) # 0 means superspecies
k=1
for(g in 1:N.genera)
{
 for(s in 1:genera.num.species[g])
 {
   species.genera[k]=g
   k=k+1
 }
}

sp=species.names=c(sprintf("%s_sp%02d",genera.names[species.genera[1:N.focal]],1:N.focal),"Superspecies")

N.site=10
N.shell=100


log.lambda.genera=log(sort(rep(exp(seq(log(0.025),log(0.1),length.out=N.genera/2)),2),
 decreasing=TRUE))

# Constant variants:
log.lambda.mult.high=log(3)
log.lambda.mult.mod=log(0.8)
log.lambda.mult.low=log(0.2)
logit.psi.high=logit(0.8)
logit.psi.mod=logit(0.6)
logit.psi.low=logit(0.4)


# Set species constants:
log.lambda.sp=rep(0,N.sp)
logit.psi.sp=rep(0,N.sp)

s=1
for(g in 1:N.genera)
{
 for(k in 1:genera.num.species[g])
 {
  if(k==1)
  {
    log.lambda.sp[s]=log.lambda.genera[g]+log.lambda.mult.high
    logit.psi.sp[s]=logit.psi.mod
  }
  if(k==2 && k!=genera.num.species[g])
  {
    log.lambda.sp[s]=log.lambda.genera[g]+log.lambda.mult.mod
    logit.psi.sp[s]=logit.psi.low
  }
  if(k>=3 && k!=genera.num.species[g])
  {
    log.lambda.sp[s]=log.lambda.genera[g]+log.lambda.mult.mod
    logit.psi.sp[s]=logit.psi.high
  }
  if(k==genera.num.species[g])
  {
    log.lambda.sp[s]=log.lambda.genera[g]+log.lambda.mult.low
    logit.psi.sp[s]=logit.psi.mod
  }
  s=s+1
 }
}

psi.sp=ilogit(logit.psi.sp)


# Superspecies:
log.lambda.sp[N.sp]=log(0.5)
psi.sp[N.sp]=1

lambda.sp=exp(log.lambda.sp)


##################
# Dynamics:
##################

# Species-wise dynamics

log.lambda.sf=array(0,c(Nf,N.sp))
logit.psi.sf=array(0,c(Nf,N.sp))

log.lambda.maxvar=log(3)
logit.occu.maxvar=logit(0.75)

increase=(1:Nf-Nf/2-.5)/(Nf/2-.5)
decrease=-increase
top=sin(2*pi*(1:Nf-1)/(Nf-1))
top=-1+2/(max(top)-min(top))*(top-min(top))
bottom=-top


s=1
for(g in 1:N.genera)
{
 for(k in 1:genera.num.species[g])
 {
   if(k>=1 && k!=genera.num.species[g])
   {
     log.lambda.sf[,s]=top*log.lambda.maxvar/2
     logit.psi.sf[,s]=increase*logit.occu.maxvar*2
   }
   
   if(k==2)
   {
     log.lambda.sf[,s]=increase*log.lambda.maxvar
     logit.psi.sf[,s]=top*logit.occu.maxvar/2
   }

   if(k==3 && k!=genera.num.species[g])
    {
     log.lambda.sf[,s]=decrease*log.lambda.maxvar
     logit.psi.sf[,s]=bottom*logit.occu.maxvar/2
   }

   if(k==genera.num.species[g])
   {
     log.lambda.sf[,s]=bottom*log.lambda.maxvar/2
     logit.psi.sf[,s]=decrease*logit.occu.maxvar
   }
  
  s=s+1
 } 
}

lambda=exp(rep(1,Nf)%*%t(log.lambda.sp)+log.lambda.sf)
psi=ilogit(rep(1,Nf)%*%t(logit.psi.sp)+logit.psi.sf)

# Superspecies:
lambda[,N.sp]=lambda.sp[N.sp]
psi[,N.sp]=psi.sp[N.sp]


# Common detection dynamics:

log.lambda.common=top*0.5


lambda.obs=lambda
for(f in 1:Nf)
  lambda.obs[f,]=exp(log(lambda[f,])+log.lambda.common[f])
  
  





psilambda=psi*lambda




# Regional occupancy:
# Leave two formations unoccupied for each species (except the Superspecies)::
occu.reg=array(1,c(Nf,N.sp))

s=1
for(g in 1:N.genera)
{
 for(k in 1:genera.num.species[g])
 {
  for(f in 1:Nf)
  {
    if((f+k)%%floor(Nf/2)==0)
      occu.reg[f,s]=0
  }
  s=s+1
 }
}



rel.abundance.real=psilambda*occu.reg
for(i in 1:Nf)
  rel.abundance.real[i,]=rel.abundance.real[i,]/sum(rel.abundance.real[i,])



# Unidentifiabilty
# Only the case for each second genera:
# Let the first genus with unidentified colonies
# have decending identifiability rate, the second
# ascending identifiability rate and the rest
# alternating between high and low identifiability
# rate

unid.gen=rep(c(FALSE,TRUE),N.genera/2)
gam.real=array(1, c(Nf,N.genera))
for(g in 2*(1:(N.genera/2)))
{
 if(g/2==1)
   gam.real[,g]=sort(ilogit(seq(logit(0.3),logit(0.9),length.out=Nf)),decreasing=FALSE)
 if(g/2==2)
   gam.real[,g]=sort(ilogit(seq(logit(0.3),logit(0.9),length.out=Nf)),decreasing=TRUE)
 if(g/2>=3)
   gam.real[,g]=rep(c(0.3,0.8),Nf/2+1)[1:Nf]
}
gam=gam.real



# Occupancy state of each formation/site/species combination
occu.state=array(0.0, c(Nf,N.site,N.sp))
for(f in 1:Nf)
 for(s in 1:N.sp)
  occu.state[f,,s]=occu.reg[f,s]*sample(1:0,N.site, replace=TRUE,
                          prob=c(psi[f,s],1-psi[f,s]))

