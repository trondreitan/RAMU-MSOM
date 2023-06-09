This folder contains the original file used for making the input file (expanded_infile.csv 
in the main folde) to RAMU-MSOM. There is one Excel file (Ecological_samples_masterfile_NEW_27.09.2021.xlsx) 
with sheets for the site, subsample and species data and another for the formation data. Also contains a 
script to read the old input file data from Reitan et al. 2022, in order to fill in some information for 
the material collected during the first field season (2014) that ended up being the data material in this 
study. 
 

(PS: A new file formation age estimates was used for the plotting, included in the main folder). The 
script that takes these different sources, checks and corrects them, collates the information to the 
site level and saves to the occupancy analysis input file ("epxanded_infile.csv") which is found 
in the main folder.

R script:

* check_and_write4.R - check_and_write4.R - Reads the input file sheets for subsamples (shells), 
individuals (colony) and sites (samples). Checks these sheets separately, and confirms that nothing is wrong. 
Makes a subsample data structure with individual information added and makes the superspecies category 
and genera level data. Then aggregates this information to the site level. Reads the old input file from the
previous occupancy study (Reitan et al. 2022) in order to fetch the total number of subsamples in each
site for the 2014 field season that is our data material (as this information was not present in our 
the current subsample data), disregarding the sites that were in the old but not the new dataset. 
Fills in the toal number of subsamples for the other years directly from the subsample data sheet. 
Fills in standardized weights and areas from the subsample information (shell). Adds genus information, 
including total number of individuals of each genus and the number of unidentified individuals belonging 
to each genus. Adds formation information from a data sheet from an older Excel file 
(Ecological_samples_masterfile_NEW_27.09.2021.xlsx). Writes the resulting data structure 
to a file that will be used as the input file for the occupancy analysis, "expanded_infile.csv".



Input data:

* Sample_Shell_Colonies_25.07.2022.xlsx - Contains 3 datasheets called "Shells" (subsamples),
"Colonies" (individuals) and "Samples" (sites). This is the main data source for our occupancy
data file. 

* Ecological_samples_masterfile_NEW_27.09.2021.xlsx - Contains an extra data sheet called "Formations" 
from where we fetch the formation information.

* allsamples_with_counts_and_metainfo.csv - Contains number of subsamples from field season of 2014 
(which was missing from our current subsample data sheet called "Shells"). 





