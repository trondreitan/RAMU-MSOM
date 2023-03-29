

library(plyr)
library(gdata)


##################
# Read shells:
##################


shell=read.xls("Sample_Shell_Colonies_25.07.2022.xlsx",
  sheet="Shells",header=T)

shell_id=sprintf("%s:%s", shell$SAMPLE_ID, shell$Shell_ID)
shell$Shell_ID=shell_id

shell$area=shell$weight^(2/3)

shell$std.weight=shell$weight/mean(shell$weight)
shell$std.area=shell$area/mean(shell$area)



#######################
# Read colonies:
#######################

colonies=read.xls("Sample_Shell_Colonies_25.07.2022.xlsx",
  sheet="Colonies",header=T)


shell_id=sprintf("%s:%s", colonies$SAMPLE_ID, colonies$Shell_ID)
colonies$Shell_ID=shell_id

species=sprintf("%s_%s", colonies$Genus, colonies$Species)
species[species=="Microporella_aff. speculum"]="Microporella_aff_speculum"
colonies$Species=species


# Check if there's any sample id's in colonies not found in shells or vice versa
col.sample=sort(unique(colonies$SAMPLE_ID))
sh.sample=sort(unique(shell$SAMPLE_ID))

setdiff(col.sample, sh.sample)
# none
setdiff(sh.sample, col.sample)
# "203" "26"
# Ask Emanuela and LH

# Check if there's any shell id's in colonies not found in shells:
col.shells=sort(unique(colonies$Shell_ID))
sh.shells=sort(unique(shell$Shell_ID))

setdiff(col.shells, sh.shells)
# none
setdiff(sh.shells, col.shells)
# quite a few (2034)

# fractional number of colonies?
sum(floor(shell$Colonies)!=ceiling(shell$Colonies))
#[1] 0
#No, luckily!




###################################
# Add colony info to the shells:
###################################

shell2=shell
species=sort(unique(colonies$Species))

shell2[,dim(shell)[2]+1:length(species)]=array(0,c(dim(shell)[1],length(species)))

names(shell2)[dim(shell)[2]+1:length(species)]=sprintf("Species_%s",species)

for(i in 1:dim(shell2)[1])
 for(j in 1:length(species))
  shell2[i, dim(shell)[2]+j]=sum(colonies$Shell_ID==shell2$Shell_ID[i] &
                                 colonies$Species==species[j])

for(j in 1:length(species))
  show(c(sum(shell2[,dim(shell)[2]+j]),
         sum(colonies$Species==species[j])))
#[1] 720 720
#[1] ...
# Look good!

# Add Superspecies:
shell2$Species_Superspecies=rep(0,dim(shell2)[1])
for(i in 1:dim(shell2)[1])
  shell2$Species_Superspecies[i]=shell2$Colonies[i]-
      shell2$Cyclostomes[i]-
      sum(shell2[i,dim(shell)[2]+1:length(species)])

# Any negatives?
shell2$Shell_ID[shell2$Species_Superspecies<0]
# none!
# 22.07.2022: [1] "16:6"   "44:17"  "208:47", yes! :(



# Code for removing negatives, not necessary now:
shell2=shell2[shell2$Species_Superspecies>=0,]



# Add genera info:
genuses=sort(unique(colonies$Genus))

shell2[,dim(shell2)[2]+1:length(genuses)]=array(0,c(dim(shell2)[1],length(genuses)))

names(shell2)[dim(shell2)[2]+1:length(genuses)-length(genuses)]=sprintf("Genus_%s",genuses)

for(i in 1:dim(shell2)[1])
 for(j in 1:length(genuses))
  shell2[i, dim(shell2)[2]+j-length(genuses)]=
    sum(colonies$Shell_ID==shell2$Shell_ID[i] &
        colonies$Genus==genuses[j])




write.table(shell2, "shell_extra.csv",sep=";",
  col.names=TRUE, row.names=FALSE) 


#######################################
# Now make sample (site) level info:
#######################################

# Check if positition's are unique to each site
sites=sort(unique(shell2$SAMPLE_ID))
loc=sprintf("%s:%s",as.character(shell2$Lat),as.character(shell2$Long))
species=sort(unique(colonies$Species))


for(i in 1:length(sites))
{
  locs=unique(loc[shell2$SAMPLE_ID==sites[i]])
  if(length(locs)>2)
    show(c(sites[i],locs))
}
# No non-unique sites. Good!


# Check if formation name is unique
for(i in 1:length(sites))
{
  forms=unique(shell2$Formation_name[shell2$SAMPLE_ID==sites[i]])
  if(length(forms)>2)
    show(c(sites[i],forms))
}

sample=read.xls("Sample_Shell_Colonies_25.07.2022.xlsx",
  sheet="Samples",header=T)


# Add species info:
orig.col=dim(sample)[2]
sample[,orig.col+1:(length(species)+5)]=
  array(0,c(dim(sample)[1],length(species)+5))
names(sample)[orig.col+1:(length(species)+5)]=
  c("Total","Colonies","std.weight","std.area",names(shell2)[min(which(substr(names(shell2),1,8)=="Species_"))-1+1:(length(species)+1)])

# Compare number of sites in shell info with sample info:
c(length(sites), dim(sample)[1])
# 144 144

setdiff(sites,sample$SAMPLE_ID)
# none

setdiff(sample$SAMPLE_ID,sites)
# None

# Code for removing sites from shell data:
extrasites=setdiff(sites,sample$SAMPLE_ID)
for(site in extrasites)
{
  shell2=shell2[shell2$SAMPLE_ID!=site,]
}



for(i in 1:dim(sample)[1])
{
  # Traverse the species colony counts
  for(j in which(substr(names(sample),1,8)=="Species_"))
  {
    show(c(i,j))
    sample[i,j]=sum(shell2[shell2$SAMPLE_ID==sample$SAMPLE_ID[i],
      which(names(shell2)==names(sample)[j])]>0)
  }
}

# Worry: Do we have all the shells, or just shells with something on them?
sum(shell2$Colonies==0 && shell2$Finichnus=="N")
# 0
# Seems like we only have shells with something on them!

# Fill Total number of shells
# Fill year 1 from old data file:
oldsamp=read.csv("allsamples_with_counts_and_metainfo.csv",sep=";",header=T)

site.1=sort(unique(shell$SAMPLE_ID[shell$year==1]))
site.old=sort(unique(oldsamp$SAMPLE_ID))

setdiff(site.1,site.old)
# None
setdiff(site.old,site.1)
#[1] "27" "33" "35" "36" "37" "39" "42" "51" "53" "62" "63" "64"
# These, Emanuela has not gone through yet.

# Fill std.weight and std.area:
for(i in 1:dim(sample)[1])
{
  j=which(shell2$SAMPLE_ID==sample$SAMPLE_ID[i])
  sample$std.weight[i]=sum(shell2$std.weight[j])
  sample$std.area[i]=sum(shell2$std.area[j])
}

# Fill Total for year 1
sample$Total=0
for(i in 1:length(site.1))
{
  j=which(sample$SAMPLE_ID==site.1[i])
  k=which(oldsamp$SAMPLE_ID==site.1[i])
  if(length(j)!=1 | length(k)!=1)
    show(i)
  if(length(j)==1 & length(k)==1)
  {
   sample$Total[j]=oldsamp$Total[k]
   numshells=sum(shell2$SAMPLE_ID==site.1[i])
   sample$std.weight[j]=sample$std.weight[j]*sample$Total[j]/numshells
   sample$std.area[j]=sample$std.area[j]*sample$Total[j]/numshells
  }
}

# Fill Total for year>1
shell23=shell2[shell2$year>1,]
site.23=sort(unique(shell23$SAMPLE_ID))
for(i in 1:length(site.23))
{
  j=which(sample$SAMPLE_ID==site.23[i])
  if(length(j)!=1)
    show(i)
  if(sample$Total[j]!=0)
    show(c(siter.23[i],j,sample$Total[j]))
  sample$Total[j]=sum(shell23$SAMPLE_ID==site.23[i])
}

# Check total vs number of shells in shell info:
for(i in 1:dim(sample)[1])
  if(sample$Total[i]<sum(shell2$SAMPLE_ID==sample$SAMPLE_ID[i]))
    show(i)




names(sample)[names(sample)=="Colonies"]="num.with.colonies"

sample=sample[,c(sort(which(names(sample)!="Species_Microporella_sp.")),
   which(names(sample)=="Species_Microporella_sp."))]
names(sample)[names(sample)=="Species_Microporella_sp."]="Unidentified_Microporella"

# Add genus information:
genuses=sort(unique(colonies$Genus))

sample[,dim(sample)[2]+1:(length(genuses))]=
  array(0,c(dim(sample)[1],length(genuses)))
names(sample)[(dim(sample)[2]-length(genuses)+1):dim(sample)[2]]=
  sprintf("Genus_%s",genuses)

for(i in 1:dim(sample)[1])
{
  # Traverse the genus colony counts
  for(j in (dim(sample)[2]-length(genuses)+1):dim(sample)[2])
    sample[i,j]=sum(shell2[shell2$SAMPLE_ID==sample$SAMPLE_ID[i],
      which(names(shell2)==names(sample)[j])]>0)
}


# Check if species/genus count is higher than total count
for(j in which(substr(names(sample),1,5)=="Genus" | 
  substr(names(sample),1,5)=="Species"))
  if(sum(sample[,j]>sample$Total)>0)
    c(names(sample)[j], which(sample[,j]>sample$Total))
# nothing




# Add formation numbers
excel.file=file.path("Ecological_samples_masterfile_NEW_27.09.2021.xlsx")
form=read.xls(excel.file, sheet="Formations", header=T, dec=",")


sample$form.nr=rep(0, dim(sample)[1])
sample$orig.form.nr=rep(0, dim(sample)[1])
sample$K1=rep(0, dim(sample)[1])
sample$K2=rep(0, dim(sample)[1])
sample$K.mid=rep(0, dim(sample)[1])
for(i in 1:dim(sample)[1])
{
  j=which(form$Formation_name==sample$Formation_name[i])
  if(length(j)!=1)
    show(i)
  if(length(j)==1)
  {
    sample$orig.form.nr[i]=form$No.[j]
    sample$K1[i]=form$K1[j]
    sample$K2[i]=form$K2.[j]
    sample$K.mid[i]=mean(c(form$K1[j],form$K2.[j]))
  }
}
uform=sort(unique(sample$orig.form.nr),decreasing=TRUE)
for(i in 1:length(uform))
  sample$form.nr[sample$orig.form.nr==uform[i]]=i

sample=sample[,c(1:which(names(sample)=="Formation_name"),
   which(names(sample)=="orig.form.nr"),
   which(names(sample)=="form.nr"),
   which(names(sample)=="K1"),
   which(names(sample)=="K2"),
   which(names(sample)=="K.mid"),
   which(names(sample)=="Lat"):max(which(substr(names(sample),1,5)=="Genus" | substr(names(sample),1,5)=="Species" | names(sample)=="Microporella_sp.")))]




write.table(sample, file="analysis_infile.csv", sep=";",
  col.names=TRUE, row.names=FALSE)












# Aggregate number of colony data (instead of number of shells with colony):
sample2=sample
for(j in which(substr(names(sample2),1,8)=="Species_"))
  sample2[,j]=NA
for(j in which(substr(names(sample2),1,6)=="Genus_"))
  sample2[,j]=NA
  

for(i in 1:dim(sample2)[1])
{
  # Traverse the species colony counts
  for(j in which(substr(names(sample2),1,8)=="Species_"))
  {
    show(c(i,j))
    sample2[i,j]=sum(shell2[shell2$SAMPLE_ID==sample2$SAMPLE_ID[i],
      which(names(shell2)==names(sample2)[j])])
  }
}

for(i in 1:dim(sample2)[1])
{
  # Traverse the species colony counts
  for(j in which(substr(names(sample2),1,6)=="Genus_"))
    sample2[i,j]=sum(shell2[shell2$SAMPLE_ID==sample2$SAMPLE_ID[i],
      which(names(shell2)==names(sample2)[j])])
}


for(i in 1:dim(sample2)[1])
  sample2$Unidentified_Microporella[i]=sum(shell2$Species_Microporella_sp.[shell2$SAMPLE_ID==sample2$SAMPLE_ID[i]])
  


write.table(sample2, file="expanded_infile.csv", sep=";",
  col.names=TRUE, row.names=FALSE)










# Sample tests:
which(sample$Species_Microporella_discors>0)
i=5
sample$SAMPLE_ID[i]
# "102"

sample$Species_Microporella_discors[i]
# 1

colonies[colonies$SAMPLE_ID=="102" & colonies$Species=="Microporella_discors",]
#    SAMPLE_ID Shell_ID Colony_ID        Genus              Species
#3624       102   102:15        20 Microporella Microporella_discors

#ok


i=1
sample$SAMPLE_ID[i]
# "1"

sample$Species_Microporella_discors[i]
# 0

colonies[colonies$SAMPLE_ID=="1" & colonies$Species=="Microporella_discors",]
# <0 rows>

# OK


i=7
sample$SAMPLE_ID[i]
# "104"

sample$Species_Microporella_discors[i]
# 3

colonies[colonies$SAMPLE_ID=="104" & colonies$Species=="Microporella_discors",]
#    SAMPLE_ID Shell_ID Colony_ID        Genus              Species
#3659       104    104:2         2 Microporella Microporella_discors
#3670       104    104:5         4 Microporella Microporella_discors
#3671       104    104:6         1 Microporella Microporella_discors



i=136
sample$SAMPLE_ID[i]
# "56"

c(sample$Species_Antarctothoa_tongima[i], sample2$Species_Antarctothoa_tongima[i])
# 3 4


colonies[colonies$SAMPLE_ID=="56" & colonies$Species=="Antarctothoa_tongima",]
#    SAMPLE_ID Shell_ID Colony_ID        Genus              Species
#4650        56     56:2         1 Antarctothoa Antarctothoa_tongima
#4651        56     56:2         2 Antarctothoa Antarctothoa_tongima
#4662        56    56:16         1 Antarctothoa Antarctothoa_tongima
#4685        56    56:27         1 Antarctothoa Antarctothoa_tongima

# OK






















