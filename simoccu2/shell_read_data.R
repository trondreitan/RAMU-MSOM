
##########################
# Fetch and process data:
##########################

y=read.csv(sprintf("sim_shell_%03d.csv",sim.nr),sep=";",header=T)

N.sp=sum(substr(names(y),1,7)=="Species")
N.ge=sum(substr(names(y),1,5)=="Genus")
N.un=sum(substr(names(y),1,12)=="Unidentified")

N.sites=dim(y)[1]
N.forms=length(unique(y$Formation_name))

logit=function(x) log(x/(1-x))
ilogit=function(x) 1/(1+exp(-x))

#####################################################

#Shortform for essential data sizes:
# Number of species:
S=N.sp
# Number of sites:
N=N.sites
# Number of formations:
Nf=N.forms



forms=sort(unique(y$form.nr))
Nf2=length(forms) # check that this is the same as Nf
Nf3=max(forms) # check that this is also the same as Nf

K1=1:Nf-0.5
K2=1:Nf+0.5
if(exists("time.start",where=y) & exists("time.end",where=y))
 for(i in 1:Nf)
 {
  K1[i]=y$time.start[y$form.nr==i][1]
  K2[i]=y$time.end[y$form.nr==i][1]
 }

K.mid=(K1+K2)/2

K=K.mid
Ks=K1
Ke=K2



sp=names(y)[substr(names(y),1,8)=="Species_"]
sp=substr(sp,9,nchar(sp))


genus=sp
sl=strsplit(sp,"_")
for(i in 1:length(sp))
  genus[i]=sl[[i]][1]
genera=unique(genus)

genus.nr=rep(1,length(sp))
for(i in 1:length(sp))
  genus.nr[i]=which(genus[i]==genera)

genusnames=names(y)[substr(names(y),1,5)=="Genus"]
genusnames=substr(genusnames,7,nchar(genusnames))

G=max(genus.nr)


unid.genus.nr=rep(1,N.un)
if(N.un>0)
 for(i in 1:N.un)
{
  genusname1=names(y)[substr(names(y),1,12)=="Unidentified"][i]
  genusname1=substr(genusname1,14,nchar(genusname1))

  unid.genus.nr[i]=which(genusname1==genusnames)
}



##### Set interaction file #####
inter=as.matrix(array(1, c(Nf, N.sp)))
 

