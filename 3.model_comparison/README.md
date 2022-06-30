# Overview
The *model_comparison.R* code is designed such that it can be run either locally or on the JASMIN HPC. The installation instructions given here are for running the code on JASMIN. Note: that preprocessing steps must be completed prior to using this code.

Comparison metrics for evaluating the basic and advanced models will be produced, as well as various plots to help the user visualise the characteristics of the data for their specific application. 

The *model_comparison.R* script is the main code for analysis and should not need to be altered by the user (unless the underlying functionality has to be changed). The *config.R* script contains variables that the user should alter depending on the experiments they wish to carry out. 


# Required Data
All of the required data can be obtained by running the *preprocessing.R* code in advance. The required data can be found in the *output data* folder associated with preprocessing. 


# JASMIN Setup 
These steps assume that the user has already completed the steps required for the *preprocessing.R* code.
 
## Package Installation
The first step is to ensure the necessary R packages are installed on JASMIN: 
 1. Load the jaspy module: `module load jasr`
 2. Start the R command line: `R`
 3. Install the required packages available on CRAN (choose to install all packages and select a server location close to yours): `install.packages(c("pdp","mgcv","fitdistrplus","fields"))`
 4. Install the dggridR package from Github: `devtools::install_github("r-barnes/dggridR", branch="install_github")`


## Directory Setup
The correct directory structure must be setup for correct functionality. 
 1. Make a project folder for the *model_comparison.R*  and  *config.R*  files to live in
 2. Within the project folder, create a data folder (specified in config)
 3. Within the data folder, create a folder named  *input data*  (make sure to use this name). Place the required files (for each species) contained within the preprocessing *output data* folder here

In summary, the directory structure should be as follows: 
    
    project folder
        -> data folder
            -> input data 
                -> gis-data.gpkg (each species)
                -> ebd_processed_output.csv (each species)
                -> covariate_checklists.csv (each species)
                -> covariate_predictions.csv (each species) 
                -> prediction-surface.tif (each species)


# Expected Output
After a run has been completed (see below), the following outputs will be generated and placed in the *output data* subfolder: 
1. Spearman's rank correlation scores: `spearman_ranks.csv`
2. MAD scores: `mad_ranks.csv`
3. Deviance explained scores: `deviance_ranks.csv`

In addition, various plots are generated to give the user insight into the characteristics of the data. These plots are placed in the *analytics* folder. The following plots are generated: 
1. Species count distributions (histograms)
2. Proxy predictor fitted smooth plot
3. Optimal time of day for observation plots (basic and advanced)
4. Relative abundance and uncertainty prediction maps (basic and advanced)

In summary, the final directory structure will be as follows: 

     project folder
            -> data folder
                -> input data 
                    -> gis-data.gpkg 
                    -> ebd_processed_output.csv 
                    -> covariate_checklists.csv 
                    -> covariate_predictions.csv 
                    -> prediction-surface.tif
                -> output data
                    -> spearman_ranks.csv
                    -> mad_ranks.csv
                    -> deviance_ranks.csv
                -> analytics 
                    -> Species count distributions
                    -> Proxy predictor fitted smooth plot
                    -> Optimal time of day for observation plots 
                    -> Relative abundance and uncertainty prediction maps



# Completing a Run
Once all of the above steps have been completed, the user is ready to start a run. JASMIN uses the SLURM system for scheduling jobs - once a job is complete the outputs described above should be available for transfer. An example workflow for submitting a job is given here:
1. Write a Slurm submission script (example given below)
2. Submit the job: `sbatch <slurm file>`
3. Check the status of the job: `sacct`
4. Check the statistics of the completed run: `seff <job ID>`
5. Check the *analytics* and *output data* folders for the relevant files and transfer them (either within JASMIN or to your local machine using `rsync`)

**Example Slurm job script**

    #!/bin/bash 
    #SBATCH --partition=test 
    #SBATCH -o %j.out 
    #SBATCH -e %j.err 
    #SBATCH --time=01:00:00 
    #SBATCH --mem=0 
    
    # Must add R to workspace 
    module add jasr 
    
    # Go to relevant folder 
    cd HuntingPressure/3.model_comparison 
    
    # Execute the job in question 
    Rscript model_comparison.R

