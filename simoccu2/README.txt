This folder contains the necessary scripts to run simulation 2 in our study, which is concerned with estimating relatve abundance using
both individual count and subsample count data, as well as raw and model estimates. The complications of unidentified individuals and
the possibility of regional absence has been removed, in order to make the comparison between individual counts and subsample counts fair
(as the subsample count model does not allow for these things). The summary scripts compares RMSEs for model and raw (ratio) estimates
for both subsample count and individual count data. Note that the model for individual counts is called "smalltop" while files concerned with
individual counts themselves are lebeled "colony" as those are the individuals in our case. Scripts concerned with subsample counts are labeled
with "shell", since those are our subsamples.

In a subfolder, "simdata", we have put the simulation files we used. In another, we will put the sub-study concerning observational errors. In another called "obs_error", we put the substudy concerning observational errors.

Most of the structure and scripts here are the same as for the empirical study, but we give complete files so that the scripts are 
self-contained in this folder.

R Scripts:

Script for making the simulations:
    
* simoccu_save.R - Defines the simulation settings. Requires a variable 'save.nr" representing the number of the simulation to be predefined (as it is in the make_sim.R" script. Then defines the variation components, calculates the matrix of site-occupancy, abundance given occupancy, regional occupancy, additional information, and identifiability probabilities, before sampling the number of individuals on the subsample level for each species (using "rpois"). This is then aggregated to site-level, boths as individual count data and as subsample count data. Then the script creates the output structure, and saves the two datasets on two files that can be used as input files for the analysis scripts. Lastly, also makes a file for additional information concerning regional occupancy.


Post-analyses scripts:
    
* rel_abundance_smalltop.R - Reads and collates analyses results for the individual count data. Calculates RMSEs for model and raw (ratio) estimates. Compares these to previously calculated RMSEs from subsample count data. Also saves the relative abundance estimates, if further analyses or plotting is wanted.

* rel_abundance_shell.R - Reads and collates analyses results for the subsample count data and calculates RMSEs for model and raw (ratio) estimates.


Files necessary for performing MCMC analyses:

* run_smalltop.R - An example file for how the individual count data analysis can be run. Basically only defines the model ("smalltop"), the simulation number, a run identifiactor and some MCMC run variables. (Number of MCMC samples, thinning and number of CPUs used. Note that the LaplacesDemon used for MCMC sampling creates huge temporary files if number of CPUs is bigger than 1. I prefer just starting many processes with the use of one CPU per process, now). it then call the routine "smalltop_3stage_cpu.R", which performs the MCMC sampling. with help from some other sub-scripts. The model is called "smalltop" since the last change to out model was to introduce species random effects and so reduce the number of top parameters.

* run_shell.R - An example file for how the subsample count data analysis can be run. Basically only defines the model ("shell"), the simulation number, a run identifiactor and some MCMC run variables. (Number of MCMC samples, thinning and number of CPUs used. Note that the LaplacesDemon used for MCMC sampling creates huge temporary files if number of CPUs is bigger than 1. I prefer just starting many processes with the use of one CPU per process, now). it then call the routine "shell_3stage_cpu.R", which performs the MCMC sampling. with help from some other sub-scripts. The model is called "shell" since those were our subsample units..

* smalltop_3stage_cpu.R - MCMC script targeted for individual count data. Uses the LaplacesDemon package (https://github.com/LaplacesDemonR/LaplacesDemon) to perform MCMC sampling with the help of some sub-scripts. Does this in 3 stages and allows for running on multiple CPUs at a time (hence the name). Two of these stages are burn-in sampling using two different MCMC methods, before the real sampling starts using a third method. These different methods were used in order to manually optimize the chance of getting convergence in the MCMC samples. This script uses the "colony_read_data.R" script to read the dataset, the "model_smalltop.R" script to get the definition of the hyperparameters, the prior distribution and the likelihood, the "init_reruns.R" script in order to have routines for trying to find the best starting parameters, the "make_laplace_wrapper.R" script for creating and interface between the prior and likelihood functions and the LaplacesDemon package and the "find_best_par.R" script for finding and returning the parameter set with the highest prior*likelihood from an MCMC sample from the burn-in phases. The best parameters from "init_reruns.R" is used in as the initial parameter set in the first phase, the best parameer set from the first phase is used as the initial parameer set in the second phase, and the best parameters from the second phase is used as the intial parameter set in the third stage.

shell_3stage_cpu.R - A variant of "smalltop_3stage_cpu.R" MCMC script targeted for subsample count data instead of individual count data. Uses the LaplacesDemon package (https://github.com/LaplacesDemonR/LaplacesDemon) to perform MCMC sampling with the help of some sub-scripts. Does this in 3 stages and allows for running on multiple CPUs at a time (hence the name). Two of these stages are burn-in sampling using two different MCMC methods, before the real sampling starts using a third method. These different methods were used in order to manually optimize the chance of getting convergence in the MCMC samples. This script uses the "shell_read_data.R" script to read the dataset, the "model_shell.R" script to get the definition of the hyperparameters, the prior distribution and the likelihood, the "init_reruns.R" script in order to have routines for trying to find the best starting parameters, the "make_laplace_wrapper.R" script for creating and interface between the prior and likelihood functions and the LaplacesDemon package and the "find_best_par.R" script for finding and returning the parameter set with the highest prior*likelihood from an MCMC sample from the burn-in phases. The best parameters from "init_reruns.R" is used in as the initial parameter set in the first phase, the best parameer set from the first phase is used as the initial parameer set in the second phase, and the best parameters from the second phase is used as the intial parameter set in the third stage.

    model_smalltop.R - This script defines the model for individual count data. See otherwise the description of the script with the same name in the main folder. Defines hyperparameters, initialization functions and prior and likelihood on the log-scale. 
    
    model_shell.R - This script defines the model for subsample count data,a nd is the same as the "abundance-focused model" in Reitan, Ergon & Liow (2022). Defines hyperparameters, initialization functions and prior and likelihood on the log-scale. 

* init_rerun.R - Called by the "...stage_cpu.R" scripts to sample initial parameters many times, and return the "best one" in terms of highest prior*likelihood. Contains two functions. "init.params.flat2" uses the "init.params.flat" function defined in "model_smalltop.R" to sample without any previously defined model. "init.params.flat2.prevmodel" instead uses "init.params.flat.prevmodel" (also defined in "model_smalltop.R") instead for sampling from a previous run.

* negbinom.R - Called by the likelihood calculation in "model_smalltop.R". Contains functions defining the negative binomial distribution for the original parametrization and our parametrization ("dnegbinom1 and "dnegbinom2" respectively) and also the zero-inflated negative binomial distribution with our parametrization, on the original scale ("dnegbinom.zero") and the log-scale ("dlnegbinom.zero").

* betabin.R - Called by the likelihood calculation in "model_shell.R". Contains functions defining the beta-binomial distribution and also the zero-inflated beta-binomial distribution, on the original scale ("dbetabin.zero") and the log-scale ("dlbetabin.zero").

* colony_read_data.R - Reads the individual count data from the file "sim_colony_<sim.nr>.R", where "sim.nr" is a predfined variable representing the simulation number. Then sets variables representing the number of species, genera, genera with unidentified individuals, formations and sites. Also defines the logit and inverse logit funciton ("logit" and "ilogit" respectively). Creates a species name vector ("sp") from the columns in the data file, and similarly with the genus names ("genusnames"). Defines a structgure called "inter" with indicator values for whether additional information suggest regional occupancy for eahc species or not. In this case, where regional occupancy is not a thing, we specify that by giving a "1" for all species+formation combinations, indicating the information for regional presence exists.

* find_best_par.R - Contains code for running through MCMC samples returned from the LaplacesDemon package and find the sampled parameter set with the highest prior*likelihood. This function comes in two versions, "find.best.par.hpc" for LaplacesDemon runs with multiple cpus (really jsut cores) run in parallel, and "find.best.par", which does the same for single CPU LaplacesDemon runs.

* make_laplace_wrapper.R - Contains method that presents the data from "read_data.R" and the prior/likelihood functions defined in "model_smalltop.R" in a format that the LaplacesDemon package understands.


shell scripts:
* make_runs - A shell script that starts analyses on all the simulated indivdual count and subsample count data files for a machine running the SLURM system. PS: only run the SLURM script for the "shell" (subsample" scripts, the "colony" (individual count) scripts has been superseded by the "smalltop" model.

* make_runs_smalltop -  A shell script that starts analyses on all the simulated individual count data files for a machine running the SLURM system.

