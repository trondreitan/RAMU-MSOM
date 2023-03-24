
# Requires data and betabin definitions to have been read already.

source("negbinom.R")


hyper.this=list(logsd.mu=0, logsd.sd=1.52847, # 95% cred within (1/20,20) 
                log.kappa.mu=0, log.kappa.sd=5,
                logit.mu=0, logit.sd=2.5, 
                log.lambda.mu=-3, log.lambda.sd=2.5)
hyper=hyper.this
hyper$N=N
hyper$N.sp=N.sp
hyper$N.un=N.un
hyper$N.ge=N.ge
hyper$Nf=N.forms
hyper$N.forms=N.forms
hyper$forms=forms
hyper$K1=K1
hyper$K2=K2
hyper$span=K2-K1
hyper$K=K
hyper$sp=sp
hyper$genera=genera
hyper$genus.nr=genus.nr
hyper$genusnames=genusnames
hyper$unid.genus.nr=unid.genus.nr
hyper$inter=inter
hyper$no.presence=no.presence
hyper$n.no.presence=n.no.presence





init.params.flat=function(hyp)
{
  N.sp=hyp$N.sp
  N=hyp$N
  N.forms=hyp$N.forms
  Nf=hyp$Nf
  N.un=hyp$N.un

  if(N.un>0)
  {
  ret=list(
   log.kappa=rnorm(1, hyp$log.kappa.mu, hyp$log.kappa.sd/4),
   
   mu.const.occu=rnorm(1,hyp$logsd.mu,hyp$logsd.sd/4),
   lsd.const.occu=rnorm(1,hyp$logsd.mu,hyp$logsd.sd/4),
   mu.const.lambda=rnorm(1,hyp$logsd.mu,hyp$logsd.sd/4),
   lsd.const.lambda=rnorm(1,hyp$logsd.mu,hyp$logsd.sd/4),

   mu.lsd.sf.occu=rnorm(1, hyp$logsd.mu, hyp$logsd.sd/4),
   lsd.lsd.sf.occu=rnorm(1,hyp$logsd.mu,hyp$logsd.sd/4),
   mu.lsd.sf.lambda=rnorm(1, hyp$logsd.mu, hyp$logsd.sd/4),
   lsd.lsd.sf.lambda=rnorm(1,hyp$logsd.mu,hyp$logsd.sd/4),
   
   log.const.lambda=rnorm(N.sp,hyp$log.lambda.mu,hyp$log.lambda.sd),
   logit.const.occu=rnorm(N.sp-1,hyp$logit.mu,hyp$logit.sd),
   logit.gamma.mu=rnorm(1, hyp$logsd.mu, hyp$logsd.sd/4),
   logit.gamma.lsd=rnorm(1,hyp$logsd.mu,hyp$logsd.sd/4),
   logit.gamma=rnorm(N.un*Nf,hyp$logit.mu,hyp$logit.sd),
   log.sd.f.occu=rnorm(1,hyper$logsd.mu, hyper$logsd.sd/4),
   log.sd.f.lambda=rnorm(1,hyper$logsd.mu, hyper$logsd.sd/4),
   log.sd.sf.occu=rnorm(N.sp-1,hyper$logsd.mu-2, hyper$logsd.sd/4),
   log.sd.sf.lambda=rnorm(N.sp,hyper$logsd.mu-2, hyper$logsd.sd/4),

   f.occu=rep(0,Nf),
   f.lambda=rep(0,Nf),
   sf.occu=matrix(0, N.forms, N.sp-1),
   sf.lambda=matrix(0, N.forms, N.sp)

   )
 }
 if(N.un==0)
 {
   ret=list(
   log.kappa=rnorm(1, hyp$log.kappa.mu, hyp$log.kappa.sd/4),
   
   mu.const.occu=rnorm(1,hyp$logsd.mu,hyp$logsd.sd/4),
   lsd.const.occu=rnorm(1,hyp$logsd.mu,hyp$logsd.sd/4),
   mu.const.lambda=rnorm(1,hyp$logsd.mu,hyp$logsd.sd/4),
   lsd.const.lambda=rnorm(1,hyp$logsd.mu,hyp$logsd.sd/4),

   mu.lsd.sf.occu=rnorm(1, hyp$logsd.mu, hyp$logsd.sd/4),
   lsd.lsd.sf.occu=rnorm(1,hyp$logsd.mu,hyp$logsd.sd/4),
   mu.lsd.sf.lambda=rnorm(1, hyp$logsd.mu, hyp$logsd.sd/4),
   lsd.lsd.sf.lambda=rnorm(1,hyp$logsd.mu,hyp$logsd.sd/4),
   
   log.const.lambda=rnorm(N.sp,hyp$log.lambda.mu,hyp$log.lambda.sd),
   logit.const.occu=rnorm(N.sp-1,hyp$logit.mu,hyp$logit.sd),
   log.sd.f.occu=rnorm(1,hyper$logsd.mu, hyper$logsd.sd/4),
   log.sd.f.lambda=rnorm(1,hyper$logsd.mu, hyper$logsd.sd/4),
   log.sd.sf.occu=rnorm(N.sp-1,hyper$logsd.mu-2, hyper$logsd.sd/4),
   log.sd.sf.lambda=rnorm(N.sp,hyper$logsd.mu-2, hyper$logsd.sd/4),

   f.occu=rep(0,Nf),
   f.lambda=rep(0,Nf),
   sf.occu=matrix(0, N.forms, N.sp-1),
   sf.lambda=matrix(0, N.forms, N.sp) )
 }
   
 ret$f.occu=rnorm(Nf, 0, exp(ret$log.sd.f.occu))
 ret$f.lambda=rnorm(Nf, 0, exp(ret$log.sd.f.lambda))
 
 for(i in 1:(N.sp-1))
  ret$sf.occu[,i]=rnorm(Nf, 0, exp(ret$log.sd.sf.occu[i])/4)
    
 for(i in 1:N.sp)
  ret$sf.lambda[,i]=rnorm(Nf, 0, exp(ret$log.sd.sf.lambda[i])/4)
 
   if(sum(hyp$inter==0)>0)
   {
    ret$logit.reg.occu=rnorm(1,hyp$logit.mu+log(2),hyp$logit.sd/4)
    ret$reg.occu.state=rep(0.5, sum(hyp$inter==0))
   }
 
 return(unlist(ret))
}



init.params.flat.prevmodel=function(hyp,prevmodel.mu,prevmodel.sd)
{
  N.sp=hyp$N.sp
  N=hyp$N
  N.forms=hyp$N.forms
  Nf=hyp$Nf
  N.un=hyp$N.un

  if(N.un>0)
  {
  ret=list(
   log.kappa=rnorm(1, prevmodel.mu[1], prevmodel.sd[1]), 

   mu.const.occu=rnorm(1,hyp$logsd.mu,hyp$logsd.sd/4),
   lsd.const.occu=rnorm(1,hyp$logsd.mu,hyp$logsd.sd/4),
   mu.const.lambda=rnorm(1,hyp$logsd.mu,hyp$logsd.sd/4),
   lsd.const.lambda=rnorm(1,hyp$logsd.mu,hyp$logsd.sd/4),

   mu.lsd.sf.occu=rnorm(1, hyp$logsd.mu, hyp$logsd.sd/4),
   lsd.lsd.sf.occu=rnorm(1,hyp$logsd.mu,hyp$logsd.sd/4),
   mu.lsd.sf.lambda=rnorm(1, hyp$logsd.mu, hyp$logsd.sd/4),
   lsd.lsd.sf.lambda=rnorm(1,hyp$logsd.mu,hyp$logsd.sd/4),
   
   log.const.lambda=rnorm(N.sp,prevmodel.mu[1+1:N.sp], prevmodel.sd[1+1:N.sp]),
   logit.const.occu=rnorm(N.sp-1,prevmodel.mu[1+N.sp+1:(N.sp-1)],
                          prevmodel.sd[1+N.sp+1:(N.sp-1)]),
   logit.gamma.mu=rnorm(1, hyp$logsd.mu, hyp$logsd.sd/4),
   logit.gamma.lsd=rnorm(1,hyp$logsd.mu,hyp$logsd.sd/4),
   logit.gamma=rnorm(N.un*Nf,prevmodel.mu[2*N.sp+1:(N.un*Nf)],
                          prevmodel.sd[2*N.sp+1:(N.un*Nf)]),
			  
   log.sd.f.occu=rnorm(1,prevmodel.mu[2*N.sp+N.un*Nf+1], prevmodel.sd[2*N.sp+1]),
   log.sd.f.lambda=rnorm(1,prevmodel.mu[2*N.sp+N.un*Nf+2], prevmodel.sd[2*N.sp+2]),

   log.sd.sf.occu=rnorm(N.sp-1,hyper$logsd.mu-2, hyper$logsd.sd/4),
   log.sd.sf.lambda=rnorm(N.sp,hyper$logsd.mu-2, hyper$logsd.sd/4),

   f.occu=rnorm(Nf,prevmodel.mu[2*N.sp+N.un*Nf+2+1:Nf],prevmodel.sd[2*N.sp+2+1:Nf]),
   f.lambda=rnorm(Nf,prevmodel.mu[Nf+2*N.sp+N.un*Nf+2+1:Nf],prevmodel.sd[Nf+2*N.sp+2+1:Nf]),

   sf.occu=matrix(0, N.forms, N.sp-1),
   sf.lambda=matrix(0, N.forms, N.sp)

   )
 }
 if(N.un==0)
 {
  ret=list(
   log.kappa=rnorm(1, prevmodel.mu[1], prevmodel.sd[1]), 

   mu.const.occu=rnorm(1,hyp$logsd.mu,hyp$logsd.sd/4),
   lsd.const.occu=rnorm(1,hyp$logsd.mu,hyp$logsd.sd/4),
   mu.const.lambda=rnorm(1,hyp$logsd.mu,hyp$logsd.sd/4),
   lsd.const.lambda=rnorm(1,hyp$logsd.mu,hyp$logsd.sd/4),

   mu.lsd.sf.occu=rnorm(1, hyp$logsd.mu, hyp$logsd.sd/4),
   lsd.lsd.sf.occu=rnorm(1,hyp$logsd.mu,hyp$logsd.sd/4),
   mu.lsd.sf.lambda=rnorm(1, hyp$logsd.mu, hyp$logsd.sd/4),
   lsd.lsd.sf.lambda=rnorm(1,hyp$logsd.mu,hyp$logsd.sd/4),
   
   log.const.lambda=rnorm(N.sp,prevmodel.mu[1+1:N.sp], prevmodel.sd[1+1:N.sp]),
   logit.const.occu=rnorm(N.sp-1,prevmodel.mu[1+N.sp+1:(N.sp-1)],
                          prevmodel.sd[1+N.sp+1:(N.sp-1)]),
			  
   log.sd.f.occu=rnorm(1,prevmodel.mu[2*N.sp+N.un*Nf+1], prevmodel.sd[2*N.sp+1]),
   log.sd.f.lambda=rnorm(1,prevmodel.mu[2*N.sp+N.un*Nf+2], prevmodel.sd[2*N.sp+2]),

   log.sd.sf.occu=rnorm(N.sp-1,hyper$logsd.mu-2, hyper$logsd.sd/4),
   log.sd.sf.lambda=rnorm(N.sp,hyper$logsd.mu-2, hyper$logsd.sd/4),

   f.occu=rnorm(Nf,prevmodel.mu[2*N.sp+N.un*Nf+2+1:Nf],prevmodel.sd[2*N.sp+2+1:Nf]),
   f.lambda=rnorm(Nf,prevmodel.mu[Nf+2*N.sp+N.un*Nf+2+1:Nf],prevmodel.sd[Nf+2*N.sp+2+1:Nf]),

   sf.occu=matrix(0, N.forms, N.sp-1),
   sf.lambda=matrix(0, N.forms, N.sp) )
 }


 ret$f.occu=rnorm(Nf, 0, exp(ret$log.sd.f.occu)/4)
 ret$f.lambda=rnorm(Nf, 0, exp(ret$log.sd.f.lambda)/4)

 for(i in 1:(N.sp-1))
  ret$sf.occu[,i]=rnorm(Nf, 0, exp(ret$log.sd.sf.occu[i])/4)
    
 for(i in 1:N.sp)
  ret$sf.lambda[,i]=rnorm(Nf, 0, exp(ret$log.sd.sf.lambda[i])/4)
 
   
   if(sum(hyp$inter==0)>0)
   {
     ret$logit.reg.occu=rnorm(1,hyp$logit.mu+log(2),hyp$logit.sd/4)
     ret$reg.occu.state=rep(0.5, sum(hyp$inter==0))
   }
   

 return(unlist(ret))
}

logprior.flat=function(par, hyp)
{
  N.sp=hyp$N.sp
  N=hyp$N
  N.forms=hyp$N.forms
  Nf=hyp$Nf
  N.un=hyp$N.un

  i=0
  
  # Get top parameters:
  log.kappa=par[(i+1)]; i=i+1

  mu.const.occu=par[(i+1)]; i=i+1
  lsd.const.occu=par[(i+1)]; i=i+1
  mu.const.lambda=par[(i+1)]; i=i+1
  lsd.const.lambda=par[(i+1)]; i=i+1

  mu.lsd.sf.occu=par[(i+1)]; i=i+1
  lsd.lsd.sf.occu=par[(i+1)]; i=i+1
  mu.lsd.sf.lambda=par[(i+1)]; i=i+1
  lsd.lsd.sf.lambda=par[(i+1)]; i=i+1
  
  log.const.lambda=par[(i+1):(i+N.sp)]; i=i+N.sp
  logit.const.occu=par[(i+1):(i+N.sp-1)]; i=i+N.sp-1
  if(N.un>0)
  {
   logit.gamma.mu=par[(i+1)]; i=i+1
   logit.gamma.lsd=par[(i+1)]; i=i+1
   logit.gamma=par[(i+1):(i+Nf*N.un)]; i=i+Nf*N.un
  }
  log.sd.f.occu=par[i+1]; i=i+1
  log.sd.f.lambda=par[i+1]; i=i+1
  log.sd.sf.occu=par[(i+1):(i+(N.sp-1))]; i=i+N.sp-1
  log.sd.sf.lambda=par[(i+1):(i+N.sp)]; i=i+N.sp
  
  
  # Get random variables ...
  f.occu=as.numeric(par[(i+1):(i+N.forms)]); i=i+N.forms
  f.lambda=as.numeric(par[(i+1):(i+N.forms)]); i=i+N.forms
  
  sf.occu=matrix(par[(i+1):(i+Nf*(N.sp-1))], ncol=N.sp-1); i=i+Nf*(N.sp-1)
  sf.lambda=matrix(par[(i+1):(i+Nf*N.sp)], ncol=N.sp); i=i+Nf*N.sp
  
  ret=0
  if(sum(hyp$inter==0)>0)
  {
    logit.reg.occu=par[(i+1)]; i=i+1
    ret=sum(dnorm(logit.reg.occu, hyp$logit.mu,hyp$logit.sd,
                  log=T))
  }

  ret=ret+sum(dnorm(log.kappa, hyper$log.kappa.mu, hyper$log.kappa.sd,
   log=T)) +
   
  sum(dnorm(mu.const.occu, 
    hyper$logit.mu, hyper$logit.sd, log=T)) +
  sum(dnorm(lsd.const.occu, 
    hyper$logsd.mu, hyper$logsd.sd, log=T)) +
  sum(dnorm(mu.const.lambda, 
    hyper$log.lambda.mu, hyper$log.lambda.sd, log=T)) +
  sum(dnorm(lsd.const.lambda, 
    hyper$logsd.mu, hyper$logsd.sd, log=T)) +

  sum(dnorm(mu.lsd.sf.occu, 
    hyper$logsd.mu, hyper$logsd.sd, log=T)) +
  sum(dnorm(lsd.lsd.sf.occu, 
    hyper$logsd.mu, hyper$logsd.sd, log=T)) +
  sum(dnorm(mu.lsd.sf.lambda, 
    hyper$logsd.mu, hyper$logsd.sd, log=T)) +
  sum(dnorm(lsd.lsd.sf.lambda, 
    hyper$logsd.mu, hyper$logsd.sd, log=T)) +

  # PS: Random effects:
  sum(dnorm(log.const.lambda, mu.const.lambda, exp(lsd.const.lambda), log=T)) +
  sum(dnorm(logit.const.occu, mu.const.occu, exp(lsd.const.occu), log=T)) +

  # Back to parameters:
  sum(dnorm(log.sd.f.occu, hyper$logsd.mu, hyper$logsd.sd, log=T)) +
  sum(dnorm(log.sd.f.lambda, hyper$logsd.mu, hyper$logsd.sd, log=T)) +
  # Another set of random effects:
  sum(dnorm(log.sd.sf.occu, mu.lsd.sf.occu, exp(lsd.lsd.sf.occu), log=T))  +
  sum(dnorm(log.sd.sf.lambda[1:(N.sp-1)],
      mu.lsd.sf.lambda, exp(lsd.lsd.sf.lambda), log=T))+
   # Superspecies has it's own prior for the variation in log-lambda
  dnorm(log.sd.sf.lambda[N.sp],hyper$logsd.mu, hyper$logsd.sd, log=TRUE)

  # Genera-wise identifiabiliy probabilities:
  if(N.un>0)
  {
    ret=ret+
      sum(dnorm(logit.gamma.mu, hyper$logit.mu, hyper$logit.sd, log=T)) +
      sum(dnorm(logit.gamma.lsd, hyper$logsd.mu, hyper$logsd.sd, log=T)) +
      sum(dnorm(logit.gamma, logit.gamma.mu, exp(logit.gamma.lsd), log=T)) 
  }

  return(ret)
}

loglik.flat=function(Y,X,par)
# PS: we're going to take data Y and X right from new.data anyhow, no
# need to send it, really...
{
  N.sp=X$N.sp
  N=X$N
  N.forms=X$N.forms
  Nf=X$Nf
  N.un=X$N.un
  inter=X$inter  
  no.presence=X$no.presence
  n.no.presence=X$n.no.presence
  genus.nr=X$genus.nr
  unid.genus.nr=X$unid.genus.nr
  
  counts=as.matrix(Y[,which(substr(names(Y),1,7)=="Species")])
  
  i=0
  
  # Get top parameters:
  log.kappa=par[(i+1)]; i=i+1
  kappa=exp(log.kappa)
  
  mu.const.occu=par[(i+1)]; i=i+1
  lsd.const.occu=par[(i+1)]; i=i+1
  mu.const.lambda=par[(i+1)]; i=i+1
  lsd.const.lambda=par[(i+1)]; i=i+1
  mu.lsd.sf.occu=par[(i+1)]; i=i+1
  lsd.lsd.sf.occu=par[(i+1)]; i=i+1
  mu.lsd.sf.lambda=par[(i+1)]; i=i+1
  lsd.lsd.sf.lambda=par[(i+1)]; i=i+1
  
  log.const.lambda=par[(i+1):(i+N.sp)]; i=i+N.sp
  logit.const.occu=par[(i+1):(i+N.sp-1)]; i=i+N.sp-1
  logit.gamma=NA
  if(N.un>0)
  { 
    logit.gamma.mu=par[(i+1)]; i=i+1
    logit.gamma.lsd=par[(i+1)]; i=i+1
    logit.gamma=matrix(par[(i+1):(i+Nf*N.un)], ncol=N.un); i=i+Nf*N.un
  }
  log.sd.f.occu=par[i+1]; i=i+1
  sd.f.occu=exp(log.sd.f.occu)
  log.sd.f.lambda=par[i+1]; i=i+1
  sd.f.lambda=exp(log.sd.f.lambda)
  log.sd.sf.occu=par[(i+1):(i+N.sp-1)]; i=i+N.sp-1
  sd.sf.occu=exp(log.sd.sf.occu)
  log.sd.sf.lambda=par[(i+1):(i+N.sp)]; i=i+N.sp
  sd.sf.lambda=exp(log.sd.sf.lambda)
  
  # Get random variables ...
  f.occu=as.numeric(par[(i+1):(i+N.forms)]); i=i+N.forms
  f.lambda=as.numeric(par[(i+1):(i+N.forms)]); i=i+N.forms
  
  sf.occu=matrix(par[(i+1):(i+Nf*(N.sp-1))], ncol=N.sp-1); i=i+Nf*(N.sp-1)
  sf.lambda=matrix(par[(i+1):(i+Nf*N.sp)], ncol=N.sp); i=i+Nf*N.sp
  
  reg.occu=0
  reg.occu.state=rep(0,0)
  if(sum(inter==0)>0)
  {
    logit.reg.occu=par[(i+1)]; i=i+1
    reg.occu=ilogit(logit.reg.occu)
    reg.occu.state=par[(i+1):(i+sum(inter==0))]; i=i+sum(inter==0)
  }
  
  ll=0
  kappa.mat=rep(1,N)%*%t(rep(kappa,N.sp))
  
  ll=ll+sum(dnorm(reg.occu.state, mean=qnorm(reg.occu),sd=1,log=T))
  
  inter.active=which(inter[,1:(N.sp-1)]==1,arr.ind=TRUE)
  inter.inactive=which(inter[,1:(N.sp-1)]==0,arr.ind=TRUE)
  
  reg.occu.state2=array(1, c(Nf,N.sp-1))
  reg.occu.state2[inter.inactive]=reg.occu.state
  
  # Random effect contributions
  ll=ll+sum(dnorm(f.occu, 0, sd.f.occu, log=T))
  f.occu.use=f.occu%*%t(rep(1,N.sp-1))
  
  ll=ll+sum(dnorm(f.lambda, 0, sd.f.lambda, log=T))
  f.lambda.use=f.lambda%*%t(rep(1,N.sp))
  
  sd.sf.occu.use=rep(1,Nf)%*%t(sd.sf.occu)
  ll=ll+sum(dnorm(sf.occu, 0, sd.sf.occu.use,log=T))
  
  sd.sf.lambda.use=rep(1,Nf)%*%t(sd.sf.lambda)
  ll=ll+sum(dnorm(sf.lambda, 0, sd.sf.lambda.use,log=T))
  
  const.occu.mat=rep(1,N)%*%t(logit.const.occu)
  occu.mat=array(1, c(N, N.sp))
  occu.mat[,(1:(N.sp-1))]=ilogit(const.occu.mat+sf.occu[Y$form.nr,] +
                                 f.occu.use[Y$form.nr,] )
  occu.mat[,(1:(N.sp-1))][reg.occu.state2[Y$form.nr,]<0]=0

  log.const.lambda.mat=rep(1,N)%*%t(log.const.lambda)
  lambda.mat=exp(log.const.lambda.mat + f.lambda.use[Y$form.nr,] + sf.lambda[Y$form.nr,])
    
  lambda.obs=lambda.mat
  if(N.un>0)
  {
   for(s in 1:(N.sp-1))
   {
    uid=which(genus.nr[s]==unid.genus.nr)
    if(length(uid)==1)
    {
      for(f in 1:Nf)
        lambda.obs[Y$form.nr==f,s]=lambda.obs[Y$form.nr==f,s]*ilogit(logit.gamma[f,uid])
     }	
   }
  }

  if(sum(lambda.obs==0 & counts>0)>0)
    return(-1e+200)

  
  nrep=Y$Total%*%t(rep(1,N.sp))
  #ll=ll+sum(dlbetabin.zero(counts, nrep, 1-exp(-lambda.obs), betabin.s.use, occu.mat))
  ll=ll+sum(dlnegbinom.zero(counts, nrep*lambda.obs, kappa.mat, occu.mat))

  if(N.un>0)
  {
  ng=as.matrix(Y[,which(substr(names(Y),1,12)==
         "Unidentified")])
  tg=as.matrix(Y[,which(substr(names(Y),1,5)==
                   "Genus")[X$unid.genus.nr]]  )
  ig=tg-ng
  gamma.mat=ilogit(matrix(logit.gamma[y$form.nr,], ncol=dim(logit.gamma)[2]))
  ll=ll+sum(dnegbinom1(ng,ig+1,1-gamma.mat,uselog=TRUE))
  }

  return(ll)
}

# NB: method to extract regional occupancy state for 
# species s, formation f:
# par[?]





