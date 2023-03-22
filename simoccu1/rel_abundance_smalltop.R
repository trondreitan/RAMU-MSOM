source("rel_abundance_smalltop_read_data.R")



c(mean(rmse.mod),mean(rmse.mod.uncorr),mean(rmse.raw), mean(rmse.raw.uncorr))
# 0.01573273 0.02076907 0.02440976 0.02418065 # 100 datasets
# 0.01570361 0.02072453 0.02428980 0.02412719 # 72 datasets

# old:
# 0.01357421 0.01888175 0.02323400 0.02320144



c(mean(log.rmse.mod),mean(log.rmse.mod.uncorr),mean(log.rmse.raw),mean(log.rmse.raw.uncorr))
# 0.3978234 0.5403042 0.5535490 0.7008994
# 0.3970074 0.5402071 0.5535101 0.6989714

# old:
# 0.3698906 0.5189491 0.5392649 0.6972314
 
 


sim.rel.abundance=0*occu.reg
mse.tot=0
for(s in 1:N.sp)
{
  sim.rel.abundance[,s]=apply(rel.abundance.colony.mod[,,s],2,mean)

  png(sprintf("rel_abundance_species_%02d.png",s),
      height=1200,width=1600)
  par(cex=2.2)
  plot(1:Nf, rel.abundance.real[,s],type="b",
     ylim=c(min(rel.abundance.colony.mod[,,s]),
            max(rel.abundance.colony.mod[,,s])),
     xlab="time",ylab="Rel.abundance",
     main=sprintf("%s",sp[s]))
  lines( (1:Nf)[occu.reg[,s]==0],
      rep(1,sum(occu.reg[,s]==0)),type="h",lwd=50,col="lightgrey")
  lines(1:Nf, rel.abundance.real[,s],type="b",lwd=3)
  lines(1:Nf, sim.rel.abundance[,s], pch=2, lty=2, lwd=3, col="red")
  for(f in 1:Nf)
    points(rep(f, N.sim)+runif(N.sim,-0.1,0.1), 
      rel.abundance.colony.mod[,f,s])
  for(i in 1:N.sim)
    mse.tot=mse.tot+sum((rel.abundance.colony.mod[i,,s]-rel.abundance.real[,s])^2)
  dev.off()
}
mse.tot=mse.tot/N.sp/Nf/N.sim
rmse.tot=sqrt(mse.tot)
rmse.tot
#  0.01584088 # Pretty comparable...



# Try to fit it in one big figure
sim.rel.abundance=0*occu.reg
png("big_rel_abundance.png",width=4000,4000)
par(mfrow=c(5,5))
maxes.raw=rep(0,N.sp)
for(s in 1:N.sp)
  maxes.raw[s]=max(rel.abundance.colony.mod[,,s])
maxes=ceiling(maxes.raw*10)/10
maxes[maxes>=0.7]=1
maxes[maxes.raw<0.05]=0.05
maxes[maxes.raw<0.02]=0.02
maxes[maxes.raw<0.01]=0.01
seqs=list()
for(s in 1:N.sp)
{
  seqs[[s]]=NA
  if(maxes[s]==0.4)
    seqs[[s]]=seq(0,0.4,0.1)
  if(maxes[s]==0.3)
    seqs[[s]]=seq(0,0.3,0.1)
  if(maxes[s]==0.2)
    seqs[[s]]=seq(0,0.2,0.1)
  if(maxes[s]==0.1)
    seqs[[s]]=seq(0,0.1,0.05)
  if(maxes[s]==0.05)
    seqs[[s]]=seq(0,0.05,0.025)
  if(maxes[s]==0.02)
    seqs[[s]]=seq(0,0.02,0.01)
  if(maxes[s]==0.01)
    seqs[[s]]=seq(0,0.01,0.005)
  if(maxes[s]==1)
    seqs[[s]]=seq(0,1,0.5)
}
par(cex.axis=9.1,cex.lab=9.1,cex.main=8.1,cex.sub=9.1,
  mar=c(17,18.2,8,5),tcl=-2,font.main=3,omi=c(0.4,2,0.6,1))
for(s in 1:N.sp)
{
  show(s)
  sim.rel.abundance[,s]=apply(rel.abundance.colony.mod[,,s],2,mean)
  plot(1:Nf, rel.abundance.real[,s],type="b",
     ylim=c(0,maxes[s]),xlim=c(1,Nf),
     xlab="",ylab="", axes=F,
     main=sprintf("%s",sp[s]))
  axis(1, at=c(1,5,10), label=as.character(c(1,5,10)), padj=+1.5)
  axis(2, at=seqs[[s]], label=seqs[[s]], padj=-0.5)
  box(lwd=5)
  lines( (1:Nf)[occu.reg[,s]==0],
      rep(1,sum(occu.reg[,s]==0)),type="h",lwd=50,col="lightgrey")
  lines(1:Nf, sim.rel.abundance[,s], lwd=6, col="red")
  lines(1:Nf, rel.abundance.real[,s],type="l",lwd=3)
  for(f in 1:Nf)
    points(rep(f, N.sim)+runif(N.sim,-0.1,0.1), 
      rel.abundance.colony.mod[,f,s])
}
mtext(side=1, "Time interval", outer=T, cex=8, padj=+0.2)
mtext(side=2, "Relative abundance", outer=T, cex=8, padj=-0.5)
dev.off()







 

plot(as.vector(rel.abundance.colony.mod),as.vector(rel.abundance.from.means))




save(rel.abundance.colony.mod, file="rel_abundance_colony_mod.Rdata")
save(rel.abundance.colony.mod.uncorr,
     file="rel_abundance_colony_mod_uncorr.Rdata")
save(rel.abundance.colony.raw, file="rel_abundance_colony_raw.Rdata")
save(rel.abundance.colony.raw.uncorr, file="rel_abundance_colony_uncorr.Rdata")

save(rel.abundance.colony.upper, file="rel_abundance_colony_upper.Rdata")
save(rel.abundance.colony.lower, file="rel_abundance_colony_lower.Rdata")




for(f in 1:Nf)
 for(s in 1:N.sp)
  if(inter[f,s]==1)
   show(min(reg.occu.mean[,f,s]))
# 1 - Good


ocs=sum(inter==0 & occu.reg==1)
tab=as.data.frame(list(f=rep(0,ocs),s=rep(0,ocs),reg.occu.mean=rep(0,ocs), lambda=rep(0,ocs), psi=rep(0,ocs)))
i=1
for(f in 1:Nf)
 for(s in 1:N.sp)
  if(inter[f,s]==0 & occu.reg[f,s]==1)
  {
   tab[i,]=c(f,s,mean(reg.occu.mean[,f,s]),lambda[f,s],psi[f,s])
   show(c(f,s,mean(reg.occu.mean[,f,s]),lambda[f,s],psi[f,s]))
   i=i+1
  }

png("regional_abundance_where_regionally_occupied.png",
  height=1200,width=1600)
par(cex=2.4)  
i=1
mr=rep(0,ocs)
for(f in 1:Nf)
 for(s in 1:N.sp)
  if(inter[f,s]==0 & occu.reg[f,s]==1)
  {
   if(i==1)
    plot(rep(i,N.sim),reg.occu.mean[,f,s],
      xlim=c(1,ocs),ylim=c(0,1))
   if(i>1)
    points(rep(i,N.sim),reg.occu.mean[,f,s])
   mr[i]=mean(reg.occu.mean[,f,s])
   i=i+1
   if(sum(reg.occu.mean[,f,s]<0.5)>0)
     show(c(i,f,s,max(reg.occu.mean[reg.occu.mean[,f,s]<0.5,f,s]), min(reg.occu.mean[,f,s])))
   }
lines(1:ocs, mr)
dev.off()

lowreg=rep(NA,0)
i=1
for(f in 1:Nf)
 for(s in 1:N.sp)
  if(inter[f,s]==0 & occu.reg[f,s]==1)
  {
    np=which(no.presence.all[,f,s]>0)
    #if(sum(reg.occu.mean[,f,s]<0.99)>0)
    if(length(np)>0)
      lowreg=c(lowreg,reg.occu.mean[np,f,s])
  }
  
png("regpresence_nositepresence.png",width=2000,height=1600)
par(cex=5,cex.main=0.9)
hist(lowreg,xlab="Pr(regional occupancy)",freq=F,ylab="Frequency",
  main="Histogram of estimated probability of regional occupancy\nfor simulation+species+formation combinations\nwith regional presence but without detected\npresence in any of the sites")
dev.off()




plot(tab$lambda,tab$reg.occu.mean)
plot(tab$psi,tab$reg.occu.mean)
plot(tab$psi*tab$lambda,tab$reg.occu.mean)


png("regional_abundance_where_regionally_absent.png",
  height=1200,width=1600)
par(cex=2.4)  
ocs=sum(inter==0 & occu.reg==0)
tab=as.data.frame(list(f=rep(0,ocs),s=rep(0,ocs),reg.occu.mean=rep(0,ocs), lambda=rep(0,ocs), psi=rep(0,ocs)))
i=1
for(f in 1:Nf)
 for(s in 1:N.sp)
  if(inter[f,s]==0 & occu.reg[f,s]==0)
  {
   if(i==1)
    plot(rep(i,N.sim),reg.occu.mean[,f,s],
      xlim=c(1,ocs),ylim=c(0,1))
   if(i>1)
    points(rep(i,N.sim),reg.occu.mean[,f,s])

 tab[i,]=c(f,s,mean(reg.occu.mean[,f,s]),lambda[f,s],psi[f,s])
   show(c(f,s,mean(reg.occu.mean[,f,s]),lambda[f,s],psi[f,s]))
   i=i+1
  }
dev.off()


lowreg2=rep(NA,0)
i=1
for(f in 1:Nf)
 for(s in 1:N.sp)
  if(inter[f,s]==0 & occu.reg[f,s]==0)
  {
    np=which(no.presence.all[,f,s]>0)
    #if(sum(reg.occu.mean[,f,s]<0.99)>0)
    if(length(np)>0)
      lowreg2=c(lowreg2,reg.occu.mean[np,f,s])
  }
  
png("regabsence_nositepresence.png",width=2000,height=1600)
par(cex=5,cex.main=0.9)
hist(lowreg2,xlab="Pr(regional occupancy)",freq=F,ylab="Frequency",
  main="Histogram of estimated probability of regional occupancy\nfor simultion+species+formation combinations\nwith regional absence.")
dev.off()


plot(tab$lambda,tab$reg.occu.mean)
plot(tab$psi,tab$reg.occu.mean)
plot(tab$psi*tab$lambda,tab$reg.occu.mean)






# Look at reconstruction of per site occupancy, for
# species,formation combinations with regional occupancy.






sim.occu.mean=0*occu.reg
png("big_occupancy.png",width=4000,4000)
par(mfrow=c(5,5))
par(cex.axis=9.1,cex.lab=9.1,cex.main=8.1,cex.sub=9.1,
  mar=c(17,18.2,8,5),tcl=-2,font.main=3,omi=c(0.4,2,0.6,1))
for(s in 1:(N.sp-1))
{
  sim.occu.mean[,s]=apply(reg.occu.mean[,,s]*psi.mean[,,s],2,mean)
  plot(1:Nf, occu.reg[,s]*psi[,s],type="b",
     ylim=c(0,1),xlim=c(1,Nf),
     xlab="",ylab="", axes=F,
     main=sprintf("%s",sp[s]))
  axis(1, at=c(1,5,10), label=c(1,5,10), padj=+1.5)
  axis(2, at=seq(0,1,0.5), label=seq(0,1,0.5), padj=-0.5)
  box(lwd=5)
  lines( (1:Nf)[occu.reg[,s]==0],
      rep(2,sum(occu.reg[,s]==0)),type="h",lwd=50,col="lightgrey")
  lines(1:Nf, sim.occu.mean[,s], lwd=6, col="red")
  lines(1:Nf, occu.reg[,s]*psi[,s],type="l",lwd=3)
  for(f in 1:Nf)
    points(rep(f, N.sim)+runif(N.sim,-0.1,0.1), 
      reg.occu.mean[,f,s]*psi.mean[,f,s])
}
mtext(side=1, "Time interval", outer=T, cex=8, padj=+0.2)
mtext(side=2, "Occupancy probability", outer=T, cex=8, padj=-0.5)
dev.off()




sim.lambda=0*occu.reg
png("big_lambda.png",width=4000,4000)
par(mfrow=c(5,5))
maxes.raw=rep(0,N.sp)
for(s in 1:N.sp)
  maxes.raw[s]=max(lambda.mean[,,s])
maxes=ceiling(maxes.raw*10)/10
maxes[maxes>=0.7]=1
maxes[maxes.raw<0.05]=0.05
maxes[maxes.raw<0.02]=0.02
maxes[maxes.raw<0.01]=0.01
seqs=list()
for(s in 1:N.sp)
{
  seqs[[s]]=NA
  if(maxes[s]==0.6)
    seqs[[s]]=seq(0,0.6,0.3)
  if(maxes[s]==0.5)
    seqs[[s]]=seq(0,0.5,0.25)
  if(maxes[s]==0.4)
    seqs[[s]]=seq(0,0.4,0.2)
  if(maxes[s]==0.3)
    seqs[[s]]=seq(0,0.3,0.1)
  if(maxes[s]==0.2)
    seqs[[s]]=seq(0,0.2,0.1)
  if(maxes[s]==0.1)
    seqs[[s]]=seq(0,0.1,0.05)
  if(maxes[s]==0.05)
    seqs[[s]]=seq(0,0.05,0.025)
  if(maxes[s]==0.02)
    seqs[[s]]=seq(0,0.02,0.01)
  if(maxes[s]==0.01)
    seqs[[s]]=seq(0,0.01,0.005)
  if(maxes[s]==1)
    seqs[[s]]=seq(0,1,0.5)
}
par(cex.axis=9.1,cex.lab=9.1,cex.main=8.1,cex.sub=9.1,
  mar=c(17,18.2,8,5),tcl=-2,font.main=3,omi=c(0.4,2,0.6,1))
for(s in 1:N.sp)
{
  sim.lambda[,s]=apply(lambda.mean[,,s],2,mean)
  plot(1:Nf, lambda[,s]*occu.reg[,s],type="b",
     ylim=c(0,maxes[s]),xlim=c(1,Nf),
     xlab="",ylab="", axes=F,
     main=sprintf("%s",sp[s]))
  axis(1, at=c(1,5,10), label=as.character(c(1,5,10)), padj=+1.5)
  axis(2, at=seqs[[s]], label=seqs[[s]], padj=-0.5)
  box(lwd=5)
  lines( (1:Nf)[occu.reg[,s]==0],
      rep(1,sum(occu.reg[,s]==0)),type="h",lwd=50,col="lightgrey")
  lines(1:Nf, sim.lambda[,s], lwd=6, col="red")
  lines(1:Nf, lambda[,s]*occu.reg[,s],type="l",lwd=3)
  for(f in 1:Nf)
    points(rep(f, N.sim)+runif(N.sim,-0.1,0.1), 
      lambda.mean[,f,s])
}
mtext(side=1, "Time interval", outer=T, cex=8, padj=+0.2)
mtext(side=2, "Abundance given occupancy", outer=T, cex=8, padj=-0.5)
dev.off()




# Identification probability

n.unid=sum(unid.gen==T)
gen.unid=which(unid.gen==T)
sim.id.mean=0*gam.real
png("big_identification.png",width=5000,height=2000)
par(mfrow=c(1,n.unid))
par(cex.axis=9.1,cex.lab=9.1,cex.main=8.1,cex.sub=9.1,
  mar=c(17,18.2,8,5),tcl=-2,font.main=3,omi=c(0.4,2,0.6,1))
for(g in 1:n.unid)
{
  sim.id.mean[,g]=apply(gam.mean[,,g],2,mean)
  plot(1:Nf, gam.real[,gen.unid[g]],type="b",
     ylim=c(0,1),xlim=c(1,Nf),
     xlab="",ylab="", axes=F,
     main=sprintf("Genus %s",genera.names[gen.unid[g]]))
  axis(1, at=c(1,5,10), label=c(1,5,10), padj=+1.5)
  axis(2, at=seq(0,1,0.5), label=seq(0,1,0.5), padj=-0.5)
  box(lwd=5)
  lines(1:Nf, sim.id.mean[,g], lwd=6, col="red")
  lines(1:Nf, gam.real[,gen.unid[g]],type="l",lwd=3)
  for(f in 1:Nf)
    points(rep(f, N.sim)+runif(N.sim,-0.1,0.1), 
      gam.mean[,f,g])
}
mtext(side=1, "Time interval", outer=T, cex=8, padj=+0.2)
mtext(side=2, "Occupancy probability", outer=T, cex=8, padj=-0.5)
dev.off()


lowreg=rep(NA,0)
i=1
for(f in 1:Nf)
 for(s in 1:N.sp)
  if(inter[f,s]==0 & occu.reg[f,s]==1)
  {
    np=which(no.presence.all[,f,s]>0)
    #if(sum(reg.occu.mean[,f,s]<0.99)>0)
    if(length(np)>0)
      lowreg=c(lowreg,reg.occu.mean[np,f,s])
  }
  
lowreg2=rep(NA,0)
i=1
for(f in 1:Nf)
 for(s in 1:N.sp)
  if(inter[f,s]==0 & occu.reg[f,s]==0)
  {
    np=which(no.presence.all[,f,s]>0)
    #if(sum(reg.occu.mean[,f,s]<0.99)>0)
    if(length(np)>0)
      lowreg2=c(lowreg2,reg.occu.mean[np,f,s])
  }
  
png("regoccu_hists.png",width=3000,height=1600)
par(mfrow=c(1,2))
par(cex=5,cex.main=0.9)
hist(lowreg,xlab="Pr(regional occupancy)",freq=F,ylab="Frequency",
  main="(a)")

hist(lowreg2,xlab="Pr(regional occupancy)",freq=F,ylab="Frequency",
  main="(b)")
dev.off()
