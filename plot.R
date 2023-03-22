model="smalltop"

source("read_data.R")

sp2=gsub("_"," ",sp)

load(sprintf("mean_ra_%s.Rdata",model));
load(sprintf("median_ra_%s.Rdata",model));
load(sprintf("lower_ra_%s.Rdata",model));
load(sprintf("upper_ra_%s.Rdata",model));
load(file="mean_abundance_given_occupancy.Rdata");
load(file="median_abundance_given_occupancy.Rdata");
load(file="lower_abundance_given_occupancy.Rdata");
load(file="upper_abundance_given_occupancy.Rdata");
load(file="mean_occupancy.Rdata");
load(file="median_occupancy.Rdata");
load(file="lower_occupancy.Rdata");
load(file="upper_occupancy.Rdata");
load(file="mean_site_occupancy.Rdata");
load(file="median_site_occupancy.Rdata");
load(file="lower_site_occupancy.Rdata");
load(file="upper_site_occupancy.Rdata");
load(file="mean_reg_occu.Rdata")

form=read.csv("formation_ages.csv",sep=";",header=T)


# RELATIVE ABUNDANCE

# Try to fit it in one big figure
maxes.raw=rep(0,N.sp)
for(s in 1:N.sp)
  maxes.raw[s]=max(upper.rel.abundance[,s])
maxes=ceiling(maxes.raw*100)/100
maxes[maxes.raw>0.1]=ceiling(maxes.raw[maxes.raw>0.1]*10)/10
maxes[maxes.raw>0.4]=ceiling(maxes.raw[maxes.raw>0.4]*2)/2
seqs=list()
for(s in 1:N.sp)
{
  seqs[[s]]=NA
  if(maxes[s]<=0.01)
    seqs[[s]]=seq(0,0.01,0.0025)
  if(maxes[s]>0.01)
    seqs[[s]]=seq(0,0.02,0.01)
  if(maxes[s]>=0.02)
    seqs[[s]]=seq(0,0.05,0.01)
  if(maxes[s]>=0.05)
    seqs[[s]]=seq(0,0.1,0.05)
  if(maxes[s]>=0.1)
    seqs[[s]]=seq(0,0.2,0.05)
  if(maxes[s]>=0.2)
    seqs[[s]]=seq(0,0.3,0.1)
  if(maxes[s]>=0.3)
    seqs[[s]]=seq(0,0.5,0.25)
  if(maxes[s]>0.5)
    seqs[[s]]=seq(0,1,0.25)
}
png("emp_big_rel_abundance.png",width=4000,4000)
par(mfrow=c(5,5))
par(cex.axis=9.1,cex.lab=9.1,cex.main=7.1,cex.sub=9.1,
  mar=c(17,18.2,8,5),tcl=-2,font.main=3,omi=c(1,2,0.6,1))
for(s in 1:N.sp)
{
  show(s)
  # PS: Want to replace formation numbers with reverse age
  plot(-form$mean.age[1:Nf], mean.rel.abundance[,s],type="b",
     ylim=c(0,maxes[s]),
     xlab="",ylab="", axes=F, lwd=5,
     main=sprintf("%s",sp2[s]))
  axis(1, at=c(-2,-1.5,-1,-0.5), label=c("2","","1",""), padj=+1.5)
  axis(2, at=seqs[[s]], label=seqs[[s]], padj=-0.5)

  #reg.abs.pos=which(inter[,s]==0 & no.presence[,s]==1)
  #if(length(reg.abs.pos)>0)
  #  lines( (1:Nf)[reg.abs.pos],
  #      rep(1,length(reg.abs.pos)),type="h",lwd=50,col="lightgrey")
  #lines(1:Nf, mean.rel.abundance[,s],type="b",lwd=5)
  
  box(lwd=5)
  lines(-form$mean.age[1:Nf], lower.rel.abundance[,s],type="b",col="darkgrey", lwd=5)
  lines(-form$mean.age[1:Nf], upper.rel.abundance[,s],type="b",col="darkgrey", lwd=5)
  lines(-form$mean.age[1:Nf], mean.rel.abundance[,s],type="b",col="black", lwd=5)
}
mtext(side=1, "Age (millions of years)", outer=T, cex=8, padj=+0.2)
mtext(side=2, "Relative abundance", outer=T, cex=8, padj=-0.5)
dev.off()



# OCCUPANCY

# Try to fit it in one big figuremaxes.raw=rep(0,N.sp)
png("emp_big_occupancy.png",width=4000,4000)
par(mfrow=c(5,5))
par(cex.axis=9.1,cex.lab=9.1,cex.main=7.1,cex.sub=9.1,
  mar=c(17,18.2,8,5),tcl=-2,font.main=3,omi=c(1,2,0.6,1))
for(s in 1:(N.sp-1))
{
  show(s)
  # PS: Want to replace formation numbers with reverse age
  plot(-form$mean.age[1:Nf], mean.occupancy[,s],type="b",
     ylim=c(0,1),
     xlab="",ylab="", axes=F, lwd=5,
     main=sprintf("%s",sp2[s]))
  #if(s<=5)
  #{
    axis(1, at=c(-2,-1.5,-1,-0.5), label=c("2","","1",""), padj=+1.5)
  #}
  #else
  #{
  #  axis(1, at=c(1,5,10), label=FALSE, padj=+1.5)
  #}
  if(s%%5==1)
  {
    axis(2, at=seq(0,1,0.5), label=seq(0,1,0.5), padj=-0.5)
  }
  else
  {
    axis(1, at=c(1,5,10), label=FALSE, padj=-0.5)
  }

  reg.abs.pos=which(inter[,s]==0 & no.presence[,s]==1)
  if(length(reg.abs.pos)>0)
    lines( (1:Nf)[reg.abs.pos],
        rep(1,length(reg.abs.pos)),type="h",lwd=50,col="#eeeeee")
 
  box(lwd=5)
  lines(-form$mean.age[1:Nf], lower.occupancy[,s],type="b",col="darkgrey", lwd=5)
  lines(-form$mean.age[1:Nf], upper.occupancy[,s],type="b",col="darkgrey", lwd=5)
  lines(-form$mean.age[1:Nf], mean.occupancy[,s],type="b",lwd=5)
}
mtext(side=1, "Age (millions of years)", outer=T, cex=8, padj=+0.2)
mtext(side=2, "Occupancy probability", outer=T, cex=8, padj=-0.5)
dev.off()


# ABUNDANCE GIVEN OCCUPANCY

# Try to fit it in one big figure
maxes.raw=rep(0,N.sp)
for(s in 1:N.sp)
maxes.raw[s]=max(c(upper.abundance.given.occupancy[,s],mean.abundance.given.occupancy[,s]))
maxes=ceiling(maxes.raw*100)/100
maxes[maxes.raw>0.1]=ceiling(maxes.raw[maxes.raw>0.1]*10)/10
maxes[maxes.raw>0.4]=ceiling(maxes.raw[maxes.raw>0.4]*2)/2
seqs=list()
for(s in 1:N.sp)
{
  seqs[[s]]=NA
  if(maxes[s]<=0.01)
    seqs[[s]]=seq(0,0.01,0.0025)
  if(maxes[s]>0.01)
    seqs[[s]]=seq(0,0.02,0.01)
  if(maxes[s]>0.02)
    seqs[[s]]=seq(0,0.05,0.01)
  if(maxes[s]>0.05)
    seqs[[s]]=seq(0,0.1,0.025)
  if(maxes[s]>0.1)
    seqs[[s]]=seq(0,0.2,0.05)
  if(maxes[s]>0.2)
    seqs[[s]]=seq(0,0.4,0.2)
  if(maxes[s]>0.4)
    seqs[[s]]=seq(0,0.5,0.25)
  if(maxes[s]>0.5)
    seqs[[s]]=seq(0,1,0.5)
  if(maxes[s]>1)
    seqs[[s]]=seq(0,3,1)
  if(maxes[s]>3)
    seqs[[s]]=seq(0,5,2.5)
  if(maxes[s]>5)
    seqs[[s]]=seq(0,10,5)
    
}
png("emp_big_abundance_given_occupancy.png",width=4000,4000)
par(mfrow=c(5,5))
par(cex.axis=9.1,cex.lab=9.1,cex.main=7.1,cex.sub=9.1,
  mar=c(17,18.2,8,5),tcl=-2,font.main=3,omi=c(1,2,0.6,1))
for(s in 1:N.sp)
{
  show(s)
  # PS: Want to replace formation numbers with reverse age
  plot(-form$mean.age[1:Nf], mean.abundance.given.occupancy[,s],type="b",
     ylim=c(0,maxes[s]),
     xlab="",ylab="", axes=F, lwd=5,
     main=sprintf("%s",sp2[s]))
  axis(1, at=c(-2,-1.5,-1,-0.5), label=c("2","","1",""), padj=+1.5)
  axis(2, at=seqs[[s]], label=seqs[[s]], padj=-0.5)
    
  reg.abs.pos=which(inter[,s]==0 & no.presence[,s]==1)
  if(length(reg.abs.pos)>0)
    lines( -form$mean.age[reg.abs.pos],
        rep(10,length(reg.abs.pos)),type="h",lwd=15,col="#eeeeee")
  
  box(lwd=5)
  lines(-form$mean.age[1:Nf], lower.abundance.given.occupancy[,s],type="b",col="darkgrey", lwd=5)
  lines(-form$mean.age[1:Nf], upper.abundance.given.occupancy[,s],type="b",col="darkgrey", lwd=5)
  lines(-form$mean.age[1:Nf], mean.abundance.given.occupancy[,s],type="b",col="black", lwd=5)
}
mtext(side=1, "Age (millions of years)", outer=T, cex=8, padj=+0.2)
mtext(side=2, "Abundance given occupancy", outer=T, cex=8, padj=-0.5)
dev.off()



# Cases where mean>upper
index=which(mean.abundance.given.occupancy>upper.abundance.given.occupancy,arr.ind=T)

for(i in 1:dim(index)[1])
{
  s=index[i,2]
  f=index[i,1]
  print.srcref(sprintf("%s, formation %d: mean=%f > upper=%f",
       sp2[s],f,mean.abundance.given.occupancy[f,s],upper.abundance.given.occupancy[f,s]))
}
#Exochella armata, formation 2: mean=3.706593 > upper=0.256636
#Exochella armata, formation 6: mean=7.259068 > upper=0.279801











# Minimal regional occupancy probability:
i=which(mean.reg.occu==min(mean.reg.occu),arr.ind=T)
f=i[1]
s=i[2]
print.srcref(sprintf("%s, formation %d: Pr(reg. occu.)=%7.2f%%",
  sp2[s],f,100*mean.reg.occu[f,s]))
# Calloporina angustipora, formation 3: Pr(reg. occu.)=  26.15%
  




