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

occu.state=array(0.0, c(Nf,N.site,N.sp))
for(f in 1:Nf)
 for(s in 1:N.sp)
  occu.state[f,,s]=sample(1:0,N.site, replace=TRUE,
                          prob=c(psi[f,s],1-psi[f,s]))


# Variation on species+site level, shell level and
# species+shell level.

s.species.site=log(3)/1.96
s.shell=log(2.2)
s.species.shell=0.1*log(2)*sqrt(N.shell)


# Errors:
p.err=c(0,0.0001,0.0003,0.001,0.003,0.01,0.03,0.1,0.3)
N.err=length(p.err)




###############
# Sampling
###############

v.species.site=array(0.0, c(Nf, N.site, N.sp))
for(site in 1:N.site)
 for(s in 1:N.sp)
  v.species.site[,site,s]=rnorm(Nf,0,s.species.site)

v.shell=array(0.0, c(Nf,N.site,N.shell))
for(site in 1:N.site)
 for(shell in 1:N.shell)
  v.shell[,site,shell]=rnorm(Nf,0,s.shell)

v.species.shell=array(0.0, c(Nf,N.site,N.shell,N.sp))
for(site in 1:N.site)
 for(shell in 1:N.shell)
  for(s in 1:N.sp)
   v.species.shell[,site,shell,s]=rnorm(Nf,0,s.shell)

# Sample on the shell level
shell.abundance=array(0.0,c(Nf,N.site,N.shell,N.sp))
shell.samples=array(0,c(Nf,N.site,N.shell,N.sp,N.err))
shell.transfers=array(0,c(Nf,N.site,N.shell,N.sp,N.err))

for(f in 1:Nf)
 for(site in 1:N.site)
  for(shell in 1:N.shell)
   for(s in 1:N.sp)
   {
     #show(c(f,site,shell,s))
     shell.abundance[f,site,shell,s]=
      occu.state[f,site,s]*lambda[f,s]*
      exp(v.species.site[f,site,s]+v.shell[f,site,shell]+
          v.species.shell[f,site,shell,s])

     shell.samples[f,site,shell,s,1]=
      rpois(1,shell.abundance[f,site,shell,s])
   }

t=0
for(f in 1:Nf)
 for(site in 1:N.site)
  for(shell in 1:N.shell)
   for(s in 1:N.sp)
    for(k in 2:N.err)
   {
     shell.samples[f,site,shell,s,k]=shell.samples[f,site,shell,s,1]

     # Undercount:
     shell.samples[f,site,shell,s,k]=shell.samples[f,site,shell,s,k]-
       rbinom(1,shell.samples[f,site,shell,s,k],p.err[k])

     # overcount:
     shell.samples[f,site,shell,s,k]=shell.samples[f,site,shell,s,k]+
       rbinom(1,shell.samples[f,site,shell,s,k],p.err[k])

     # misclassifcy:
     mis=rbinom(1,shell.samples[f,site,shell,s,k],p.err[k])
     if(mis>0)
      for(m in 1:mis)
       {
         p=apply(occu.state[f,,],2,mean)*lambda[f,s]
	 exclude=setdiff(1:N.sp,s)
	 p=p[exclude]
	 p=p/sum(p)

         s.new=sample(exclude,1,prob=p)

         shell.transfers[f,site,shell,s.new,k]=shell.transfers[f,site,shell,s.new,k]+1
	 t=t+1
       }
     shell.samples[f,site,shell,s,k]=shell.samples[f,site,shell,s,k]-mis
   }

shell.samples=shell.samples+shell.transfers


shell.total.samples=array(0,c(Nf,N.site,N.shell,N.err))
for(f in 1:Nf)
 for(site in 1:N.site)
  for(shell in 1:N.shell)
   for(k in 1:N.err)
   shell.total.samples[f,site,shell,k]=
     sum(shell.samples[f,site,shell,,k])


# Aggregate to site level, colony counts and 
# occupied shell counts:
site.colony.count=array(0,c(Nf,N.site,N.sp,N.err))
site.shell.count=array(0,c(Nf,N.site,N.sp,N.err))
for(f in 1:Nf)
 for(site in 1:N.site)
  for(s in 1:N.sp)
   for(k in 1:N.err)
  {
   site.colony.count[f,site,s,k]=sum(shell.samples[f,site,,s,k])
   site.shell.count[f,site,s,k]=sum(shell.samples[f,site,,s,k]>0)
  }

max(shell.samples) # Should be around 100
max(shell.total.samples) # should be around 117
max(site.shell.count) # should be around 92


for(k in 1:N.err)
{
colony.data=as.data.frame(array(0,c(Nf*N.site,5+N.sp)))
names(colony.data)=c("SAMPLE_ID","Formation_name","form.nr","Total","num.with.colonies",sprintf("Species_%02d", 1:(N.sp-1)),
"Species_Superspecies")
colony.data$SAMPLE_ID=as.character(colony.data$SAMPLE_ID)
colony.data$Formation_name=
  as.character(colony.data$Formation_name)
for(f in 1:Nf)
 for(site in 1:N.site)
 {
  i=site+(f-1)*N.site
  colony.data$SAMPLE_ID[i]=sprintf("Site_f%02d_s%02d", f,site)
  colony.data$Formation_name[i]=sprintf("F%02d",f)
  colony.data$form.nr[i]=f
  colony.data$Total[i]=N.shell
  colony.data$num.with.colonies[i]=
     sum(shell.total.samples[f,site,,k]>0)

  for(s in 1:N.sp)
  {
    colony.data[i,5+s]=site.colony.count[f,site,s,k]
  }
 }


shell.data=as.data.frame(array(0,c(Nf*N.site,5+N.sp)))
names(shell.data)=c("SAMPLE_ID","Formation_name","form.nr","Total","num.with.colonies",sprintf("Species_%02d", 1:(N.sp-1)),
"Species_Superspecies")
shell.data$SAMPLE_ID=as.character(shell.data$SAMPLE_ID)
shell.data$Formation_name=
  as.character(shell.data$Formation_name)
for(f in 1:Nf)
 for(site in 1:N.site)
 {
  i=site+(f-1)*N.site
  shell.data$SAMPLE_ID[i]=sprintf("Site_f%02d_s%02d", f,site)
  shell.data$Formation_name[i]=sprintf("F%02d",f)
  shell.data$form.nr[i]=f
  shell.data$Total[i]=N.shell
  shell.data$num.with.colonies[i]=
     sum(shell.total.samples[f,site,,k]>0)

  for(s in 1:N.sp)
  {
    shell.data[i,5+s]=site.shell.count[f,site,s,k]
  }
 }


if(!exists("save.nr"))
{
  show("save.nr not given!")
}


if(exists("save.nr"))
{
  write.table(colony.data, 
    sprintf("sim_colony_s%03d_e%02d.csv",save.nr,k), 
    col.names=T, row.names=F, sep=";")
  write.table(shell.data, 
    sprintf("sim_shell_s%03d_e%02d.csv",save.nr,k), 
    col.names=T, row.names=F, sep=";")
}
}







