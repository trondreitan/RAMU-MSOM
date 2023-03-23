# RAMU-MSOM - sumlation 2 files
Relative Abundance estimation with Unidentified individuals using Multi-Species Occupancy Modelling

This folder contains the necessary scripts and data files to
run our analyses, as well as running your own analysis (by replacing 
our files with your own).


Data files:

* expanded_infile.csv - This is the individuaol count 
    occupancy analysis data file. Contains one row per site.
    There is one column for each species, where the cells
    has the number of individuals for each site for each species. 
    Each row also has a formation number, a formation name
    (which can be matched with other data files), the total
    number of subsamples (shells), the sum weight and standardized
    area and volume of these shells and the number of unidentified
    individuals belonging to the Microporella genus. See the
    "empirical dataset" folder for how this file was created.

    The following columns are (so far) manditory:
    * form.nr - The number of the formation/time period. Should
      run from 1 to number of formations. 
    * Formation_name - The name of the formation/time interval.
      Make sure that the formation names and the formation numbers 
      are consistent! Formation names are used for creating the 
      additional information file in the "collate_interactions" 
      folder, though not in the analysiss cripts themselves. 
      The formation names are however defined in the "read_data.R"
      script, so the analysis will throw an error as the code is now.
      This could however be remedied in future versions, in which
      case you can create the additional information file without
      formation names, this column can be omitted. 
    * K1/K2/K.mid - Should contain the start, end and middle of the 
      formations/time intervals. This information is collated in the 
      "read_data.R" script so is manditory so far. However, it is 
      not used in the analysis scripts, so could be omitted in a 
      future version.
    * Total - The number of subsamples at the site.
    * std.weight/std.area - Used for the alternative scaling
      analysis in the SI, but not in the analysis presented here.
      However, the "read_data.R" script re-standardizes these
      values, so the script will crash without these columns. Could
      be omitted in future versions.    
    * Species_<species_name> - This is a whole set of columns, one 
      for each species. <species_name> is here the name of each species,
      written as <genus>_<species>, where <genus> is the genus name
      and <species> is the particular species name within the genus.
      So "Species_Microporella_speculum" is the column for
      the Microporella speculum species. Each cell should contain the
      the number of individuals found of that particular species
      for that particular site.
    * Genus_<genus> - This is also a whole set of columns, one for 
      each genus, where the name of the genus is <genus>. So
      "Genus_Microporella" represents the Microprella species.
      Should contain the number of individuals found in each genus
      for each site. Used for calculating the number of identified
      individuals in a genus. Thus only necessary if there are genera
      with unidentified individuals.  
    * Unidentified_<genus> - One column for each genera with
      unidentified individuals, containing the number of unidentified 
      individuals belonging to the genus for each site. <genus> is 
      the genus name, so for instance "Unidentified_Microporella" 
      contains the number of unidentified individuals belonging to the
      genus Microporella. Thus only necessary if there are genera
      with unidentified individuals. 

* interaction_indicators.csv - Contains one row per formation
    and one column per species that there is additional information for.
    Each cell indicates whether there exists additional information
    suggesting regional presence for a species in a formation (1) or 
    not (0). See the "collate_interactions" folder for how this
    data file was created.
  
R Scripts:

1) Files necessary for performing MCMC analyses:

* run_smalltop.R - An example file for how the analysis can
    be run. Basically only defines the model ("smalltop")
    a run identifiactor and some MCMC run variables. (Number of 
    MCMC samples, thinning and number of CPUs used. Note that 
    the LaplacesDemon used for MCMC sampling creates huge temporary 
    files if number of CPUs is bigger than 1. I prefer just 
    starting many processes with the use of one CPU per process, 
    now). it then call the routine "3stage_cpu.R", which performs
    the MCMC sampling. with help from some other sub-scripts.
    The model is called "smalltop" since the last change to out 
    model was to introduce species random effects and so reduce the
    number of top parameters.

* 3stage_cpu.R - Uses the LaplacesDemon package 
    (https://github.com/LaplacesDemonR/LaplacesDemon) to
    perform MCMC sampling with the help of some sub-scripts.
    Does this in 3 stages and allows for running on multiple 
    CPUs at a time (hence the name). TTwo of these stages are
    burn-in sampling using two different MCMC methods, before
    the real sampling starts using a third method. These different
    methods were used in order to manually optimize the 
    chance of getting convergence in the MCMC samples. 
    This script uses the "read_data.R" script to read the dataset,
    the "model_smalltop.R" script to get the definition of
    the hyperparameters, the prior distribution and the likelihood,
    the "init_reruns.R" script in order to have routines for
    trying to find the best starting parameters, the
    "make_laplace_wrapper.R" script for creating and interface
    between the prior and likelihood functions and the 
    LaplacesDemon package and the "find_best_par.R" script for
    finding and returning the parameter set with the highest 
    prior*likelihood from an MCMC sample from the burn-in phases.
    The best parameters from "init_reruns.R" is used in as the 
    initial parameter set in the first phase, the best parameer 
    set from the first phase is used as the initial parameer set 
    in the second phase, and the best parameters from the second 
    phase is used as the intial parameter set in the third stage.
 
* model_smalltop.R - This script defines the model. (There are many
    scripts called "model_<modelname>.R" in our working directory,
    representing different stages in our development of the model
    we ended up with. One can copy and modify this script and call it
    something starting with "model_<...>" where "<...> stands for the 
    new model name. and the other scripts can then
    call that model instead of the one given here. This setup made
    for rapid testing and running of new models, while still keeping
    the old ones. The model we ended up with is called "smalltop", 
    since the last thing we did was to make species constants into 
    species random effects and thus create a small set of top 
    parameters. The script contains 4 features:  
    1) The definition of the hyperparameters (i.e. the variables 
    that define the prior distribution of the top parameters). This
    is a list called "hyper" that contains the hyperparameters, but
    also some additional information that is sent to the prior
    and the likelihood (such as the number of species, number of 
    genera, number of genera with unidentified individuals, number 
    of formations, a set of indicators of whether a species has no 
    detected presence in each formation ("no.presence") and another 
    set of indicators for each species of whether there is additional 
    information suggesting regional presence or not for each formation 
    ("inter"). 
    2) A set of initial value functions called "init.params.flat"
    and "init.params.flat.previousmodel". "init.params.flat" 
    is called in the "init_reruns.R" script if no previous model is 
    specified, while "init.params.flat.previousmodel" is called if a
    previous model is specified. The object of these two functions 
    are to sample an intitial value for the parameter set.
    "init.params.flat" does this completely at random (within the 
    bounds of hyperparameters. "init.params.flat.previousmodel" sets
    the parameters/random factors that the model has in common with the
    previous model from the mean and standard deviation of MCMC
    runs for the previous model. The parameters that are not common
    are set at random. For our current run, we did not use
    parameter samples from the previous model, as the model seemed
    to converge anyhow. However, in the past we have found the use
    of previous (simpler) models useful.
    3) The prior distribution function, called "logprior.flat". 
    reads a parameter vector and the hyperparameter list and assigns 
    the different values to different top parameters. Then calculates 
    the logarithmic of the probability density of the prior 
    distribution for the specified values of the top parameter.
    4) The likelihood function, called "loglik.flat", reads a 
    parameter vector, a data file and a set of additional information 
    (from the "hyper" list) and outputs the logartihmic likelihood
    value. It first reads the parameter and random effect values 
    from the parameter vector, then calculates the likelihood
    contribution from the random effect values (this could instead
    be in the prior distribution function, but we opted for putting 
    it here). It then derives from the random effect values the 
    occupancy and abundance given occupancy matrices (size: number of 
    sites x number of species). It then derives the regional occupancy 
    state from the regional occupancy random effects plus the 
    additional data indicating regional presence and modifies the 
    occupancy matrix accordingly. It also creates a modified abundance 
    given occupancy matrix which takes into account identification 
    probabilities of genera with unidentified individuals. It then 
    calculates the log-likelihood value for the individual species 
    counts based on teh occupancy and modified abundance given 
    occupancy matrices. Lastly, if there are genera with unidentified 
    individuals, it calculates the likelihood contribution for the
    unidentified individuals given the identified ones.
    
* init_rerun.R - Called by "stage_cpu.R" to sample initial parameters
    many times, and return the "best one" in terms of 
    highest prior*likelihood. Contains two functions.
    "init.params.flat2" uses the "init.params.flat" function defined 
    in "model_smalltop.R" to sample without any previously defined 
    model. "init.params.flat2.prevmodel" instead uses 
    "init.params.flat.prevmodel" (also defined in "model_smalltop.R")
    instead for sampling from a previous run. 

* negbinom.R - Called by "model_smalltop.R". Contains functions
    defining the negative binomial distribution for the original 
    parametrization and our parametrization ("dnegbinom1 and
    "dnegbinom2" respectively) and also the zero-inflated 
    negative binomial distribution with our parametrization, on 
    the original scale ("dnegbinom.zero") and the log-scale
    ("dlnegbinom.zero").

* read_data.R - Reads the individual count data from the file 
    "expanded_infile.R". Then standardizes the subsample weights
    and areas (if applicable) and sets variables representing
    the number of species, genera, genera with unidentified 
    individuals, formations and sites. Also defines the logit
    and inverse logit funciton ("logit" and "ilgoit" respectively). 
    Calculates the midpoint of each formation age. Creates 
    a species name vector ("sp") from the columns in the data file,
    and similarly with the genus names ("genusnames"). Creates 
    a lookup table for the genus numbers of genera with
    unidentified individuals ("unid.genus.nr"). Then reads
    a file called "interaction_indicators.csv" with indicator 
    values for whether additional information suggest regional 
    occupancy for eahc species or not. Fills in with zeros the 
    species that appear in the individual count data but not
    in the additional information file. Lastly, it creates
    a number of formation x number of species matrix called 
    "no.presence" which indicates whether no regional presence
    was found in the individual count data for each formation+species
    combination. "1" means no presence was found, while "0" means
    presence was found. The variable "n.no.presence" is the seet to
    the sum of this matrix, so represents the number of 
    species+formation combinations with no deteciton in the 
    individual count data.

* find_best_par.R - Contains code for running through MCMC samples
    returned from the LaplacesDemon package and find the
    sampled parameter set with the highest prior*likelihood.
    This function comes in two versions, "find.best.par.hpc"
    for LaplacesDemon runs with multiple cpus (really jsut cores) 
    run in parallel, and "find.best.par", which does the same
    for single CPU LaplacesDemon runs. 

* make_laplace_wrapper.R - Contains method that presents the data 
    from "read_data.R" and the prior/likelihood functions defined 
    in "model_smalltop.R" in a format that the LaplacesDemon 
    package understands.
    
Scripts that collates the output from the MCMC sampling:
    
* lookat.R - Reads in the MCMC samples from multiple runs. The output
    from 3stage_cpu.R are store on files called <model>_<run-number>.Rdata, 
    so for instance the first output file from the "smalltop" model (the 
    model presented in our study) is called "smalltop_01.Rdata". In the
    current script, 100 MCMC chains are fetched, though it is easy enough 
    to change that, just set the "N.runs" variable to the number of output 
    files you have. The script looks at the Gelman index, and plots the
    4 "worst" parameters according to that index and calculates the 
    autocorrelation length of those. Store parameter means and 95% credibility 
    intervals, in case of later use in other scripts. Plots MCMC samples
    with prior and also compares those analytically for all top parameters.
    Outputs some statistics about regional occupancy estimates.

* rel_abundance_smalltop.R - Reads in the MCMC samples from multiple runs.
    (See "lookat.R" for how this is done). Then collates all the MCMC chains
    into one, and fetches random variable values in order to piece together
    posterior samples of site-occupancy, identification and regional occupancy 
    probabilities as well as abundance given occupancy for each species+formation
    combinations. This again is used for getting the posterior (MCMC) samples of 
    relative abundance for each species+formation. These samples are then further 
    collated by finding the mean (estimates) and 2.5% and 97.5 quantiles (limits
    of the 95% credibility bands) for each of these quantities for each 
    species+formation combination. These are then saved on files, for later 
    plotting or analyses. A preliminary plot of relative abundance (one plot for
    each species) is done at the end of the script, to show how these things can
    be plotted.
    
    
    
