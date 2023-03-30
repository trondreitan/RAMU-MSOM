This folder contains the model script, the run script and the MCMC summary script for the weighted analysis of number of subsamples versus standardized volume and area derived from weight presented in the SI of the ms only. Other analyses scripts from the main folder are need to run this model. In addition to showing how this study was performed, this also gives an example of an alternative model.

R scripts:
* model_weight_exposure.R - A model that makes a weighted average of number of subsamples, standardized volume and standardized area, with the weighted average parameters determining how much each of these 3 factors scales with the expected value.

* run_weight_expsoure.R - An example run script for the model. The scripts "model_weight_exposure.R", "3stage_cpu.R" (found in the main folder) and the scripts that "3stage_cpu.R" relies on are needed in order to run this.

* lookat_weight_exposure.R - Fetches MCMC samples and analyzing them, looking specifically at the parameters representing the weighted average of umber of subsamples, standardized volume and standardized area.
