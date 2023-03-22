source("simoccu_settings.R")





# Variation on species+site level, shell level and
# species+shell level.

s.species.site=log(3)/1.96
s.shell=log(2.2)
s.species.shell=0.1*log(2)*sqrt(N.shell)

v.species.site=array(0.0, c(Nf, N.site, N.sp))
for(site in 1:N.site)
 for(s in 1:N.sp)
  v.species.site[,site,s]=rnorm(Nf,-0.5*s.species.site^2,s.species.site)

v.shell=array(0.0, c(Nf,N.site,N.shell))
for(site in 1:N.site)
 for(shell in 1:N.shell)
  v.shell[,site,shell]=rnorm(Nf,-0.5*s.shell^2,s.shell)

v.species.shell=array(0.0, c(Nf,N.site,N.shell,N.sp))
for(site in 1:N.site)
 for(shell in 1:N.shell)
  for(s in 1:N.sp)
   v.species.shell[,site,shell,s]=rnorm(Nf,-0.5*s.species.shell^2,s.species.shell)

# Sample on the shell level
shell.abundance=array(0.0,c(Nf,N.site,N.shell,N.sp))
shell.samples=array(0,c(Nf,N.site,N.shell,N.sp))

for(f in 1:Nf)
 for(site in 1:N.site)
  for(shell in 1:N.shell)
   for(s in 1:N.sp)
   {
     #show(c(f,site,shell,s))
     shell.abundance[f,site,shell,s]=
      occu.state[f,site,s]*lambda.obs[f,s]*
      exp(v.species.site[f,site,s]+v.shell[f,site,shell]+
          v.species.shell[f,site,shell,s])

     shell.samples[f,site,shell,s]=
      rpois(1,shell.abundance[f,site,shell,s])
   }

shell.total.samples=array(0,c(Nf,N.site,N.shell))
for(f in 1:Nf)
 for(site in 1:N.site)
  for(shell in 1:N.shell)
   shell.total.samples[f,site,shell]=
     sum(shell.samples[f,site,shell,])


# Aggregate to site level, colony counts:
site.colony.count=array(0,c(Nf,N.site,N.sp))
for(f in 1:Nf)
 for(site in 1:N.site)
  for(s in 1:N.sp)
  {
   site.colony.count[f,site,s]=sum(shell.samples[f,site,,s])
  }

for(s in 1:N.sp)
{
 if(max(site.colony.count[,,s])==0)
  print.srcref(sprintf("%s: min=%d max=%d mean=%f", sp[s],min(site.colony.count[,,s]),
     max(site.colony.count[,,s]),mean(site.colony.count[,,s])))
}





#max(shell.samples) # Should be around 100
#max(shell.total.samples) # should be around 117

colony.data=as.data.frame(array(0,c(Nf*N.site,5+N.sp)))
names(colony.data)=c("SAMPLE_ID","Formation_name","form.nr","Total","num.with.colonies",sprintf("Species_%s", sp))
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
     sum(shell.total.samples[f,site,]>0)

  for(s in 1:N.sp)
  {
    colony.data[i,5+s]=site.colony.count[f,site,s]
  }
 }



colony.data.unid=colony.data
nc=dim(colony.data)[2]
colony.data.unid[,(nc+1):(nc+sum(unid.gen)+N.genera)]=0
names(colony.data.unid) [(nc+1):(nc+sum(unid.gen))]=
  sprintf("Unidentified_%s", genera.names[unid.gen==TRUE])
names(colony.data.unid) [(nc+sum(unid.gen)+1):(nc+sum(unid.gen)+N.genera)]=
  sprintf("Genus_%s", genera.names)

for(f in 1:Nf)
 for(site in 1:N.site)
{
 i=site+(f-1)*N.site
 s=1
 unid.g=0
 for(g in 1:N.genera)
 {
  if(unid.gen[g])
    unid.g=unid.g+1
     
  unid=0
  gen=0
  for(k in 1:genera.num.species[g])
  {
   if(unid.gen[g])
   {
     #show(c(f,site,g,k,s))
     id=rbinom(1, size=colony.data[i,5+s], prob=gam[f,g])
     unid=unid+colony.data[i,5+s]-id
     colony.data.unid[i,5+s]=id
   }
   gen=gen+colony.data[i,5+s]
   
   s=s+1
  }
  
  if(unid.gen[g])
    colony.data.unid[i,nc+unid.g]=unid
  colony.data.unid[i,nc+sum(unid.gen)+g]=gen

 }
}

# Check consistency between species sum, genus sum and unidentified colonies sum
for(i in 1:(Nf*N.site))
{
 unid.g=0
 for(g in which(unid.gen==TRUE))
 {
   unid.g=unid.g+1
   id=sum(colony.data.unid[i,5+which(species.genera==g)])
   unid=colony.data.unid[i,nc+unid.g]
   total=colony.data.unid[i,nc+sum(unid.gen)+g]

   if(id+unid!=total)
   {
     show(c(i,g,id,unid,total))
   }
 }
 
 for(g in which(unid.gen==FALSE))
 {
   id=sum(colony.data.unid[i,5+which(species.genera==g)])
   total=colony.data.unid[i,nc+sum(unid.gen)+g]

   if(id!=total)
   {
     show(c(i,g,id,total))
   }
 }
}
# No output= everything good


if(!exists("save.nr"))
{
  show("save.nr not given!")
}


if(exists("save.nr"))
{
  write.table(colony.data.unid, 
    sprintf("sim_colony_%03d.csv",save.nr), 
    col.names=T, row.names=F, sep=";")
}




# make and save the extra source matrix:
inter=array(1,c(Nf,N.sp))
for(f in 1:Nf)
 for(s in 1:(N.sp-1))
   inter[f,s]=occu.reg[f,s]*((f+s)%%2==1)
inter[f,N.sp]=1

write.table(inter, sprintf("inter_%03d.csv",save.nr), 
    col.names=T, row.names=F, sep=";")







# make and save the extra source matrix:
inter=array(1,c(Nf,N.sp))
for(f in 1:Nf)
 for(s in 1:(N.sp-1))
   inter[f,s]=((f+s)%%2==1)
inter[f,N.sp]=1
