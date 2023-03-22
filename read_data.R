
##########################
# Fetch and process data:
##########################

y=read.csv("expanded_infile.csv",sep=";",header=T)
y$std.weight=(y$std.weight/y$Total)/mean(y$std.weight/y$Total)*y$Total
y$std.area=(y$std.area/y$Total)/mean(y$std.area/y$Total)*y$Total


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

forms=rep("",Nf)
for(i in 1:Nf)
  forms[i]=as.character(y$Formation_name[min(which(y$form.nr==i))])
  

K1=rep(NA,Nf)
K2=rep(NA,Nf)
for(i in 1:Nf)
{
  K1[i]=y$K1[y$form.nr==i][1]
  K2[i]=y$K2[y$form.nr==i][1]
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
for(i in 1:N.un)
{
  genusname1=names(y)[substr(names(y),1,12)=="Unidentified"][i]
  genusname1=substr(genusname1,14,nchar(genusname1))

  unid.genus.nr[i]=which(genusname1==genusnames)
}



##### Read interaction file #####

inter.r=read.table("interaction_indicators.csv",sep=";",
  header=TRUE)
inter.r$Superspecies=1
o.rows=1:Nf
for(i in 1:Nf)
  o.rows[i]=which(forms[i]==row.names(inter.r))
inter.r=inter.r[o.rows,]

inter=array(0,c(Nf,N.sp))
for(i in 1:N.sp)
{
  j=which(sp[i]==names(inter.r))
  if(length(j)==1)
    inter[,i]=inter.r[,j]
}
inter=as.matrix(inter)



no.presence=array(0, c(Nf,N.sp))
for(s in 1:(N.sp-1))
 for(f in 1:Nf)
 {
   j=which(names(y)==sprintf("Species_%s",sp[s]))
   if(sum(y[y$form.nr==f,j])==0)
     no.presence[f,s]=1
 }

n.no.presence=sum(inter==0 & no.presence==1)






