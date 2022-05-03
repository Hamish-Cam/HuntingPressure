# Get Abundance Map
The *get_abundance.R* code is written in the R language and is designed such that it can be run either locally or on the JASMIN HPC. The installation instructions given here are for running the code on JASMIN. 

Following the steps outlined here, after having completed the preprocessing steps, will result in a relative abundance map for the users species of choice. In addition, various plots are generated to help the user visualise the characteristics for their specific species of interest. 

The *get_abundance.R* script is the main code for preprocessing the data and should not need to be altered by the user (unless the model behaviour or characteristics of a given species are very different to most). The *config.R* script contains variables that the user should alter depending on the experiments they wish to carry out. 


# Required Data
All of the required data for obtaining a relative abundance map can be obtained by running the *preprocessing.R* code previously. The required data can be found in the *output data* folder associated with this process. 


# JASMIN Setup 
These steps assume that the user has already completed the steps required for the *preprocessing.R* code.
 
## Package Installation
The first step is to ensure the necessary R packages are installed on JASMIN: 
 1. Load the jaspy module: `module load jasr`
 2. Start the R command line: `R`
 3. Install the required packages available on CRAN (choose to install all packages and select a server location close to yours): `install.packages(c("pdp","mgcv","fitdistrplus","fields"))`
 4. Install the dggridR package from Github: `devtools::install_github("r-barnes/dggridR", branch="install_github")`


## Directory Setup
Within the project folder, create a folder for the abundance map data with the name *abun-SpeciesName*.  Make a subfolder within here named *input data* and copy the required data from the preprocessing step into this folder. 

In summary, the directory structure should be as follows: 
    
    project folder
        -> abun-<species name> folder
            -> input data folder
                -> gis-data.gpkg 
                -> ebd_processed_output.csv 
                -> landcover_and_elevation_checklists.csv 
                -> landcover_and_elevation_prediction.csv 
                -> prediction-surface.tif
                -> prediction-surface.tif.aux.xml


# Expected Output
After a run has been completed (see below), the following outputs will be generated and placed in the *output data* subfolder: 
1. TIFF file containing the values of relative abundance across the area of interest: `abundance_map.tif`
2. TIFF file containing the values of standard error across the area of interest: `abundance_map_uncertainty.tif`

In addition, various plots are generated to give the user insight into the characteristics of the data. These plots are placed in the *analytics* folder. The following plots are generated: 
1. Species count distributions (histograms)
2. Model covariate dependencies 
3. Optimal time of day for observation
4. Relative abundance and uncertainty maps

In summary, the final directory structure will be as follows: 
            
  

     project folder
            -> abun-<species name> folder
                -> input data folder
                    -> gis-data.gpkg 
                    -> ebd_processed_output.csv 
                    -> landcover_and_elevation_checklists.csv 
                    -> landcover_and_elevation_prediction.csv 
                    -> prediction-surface.tif
                    -> prediction-surface.tif.aux.xml
                -> output data
                    -> abundance_map.tif
                    -> abundance_map_uncertainty.tif
                -> analytics 
                    -> count_distributions.png
                    -> model_covariates_dependency.png
                    -> optimal_time_of_day.png
                    -> abundance_maps.png



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
    #SBATCH --time=04:00:00 
    #SBATCH --mem=0 
    
    # Must add R to workspace 
    module add jasr 
    
    # Go to relevant folder 
    cd HuntingPressure/data_analysis/eBird 
    
    # Execute the job in question 
    Rscript get_abundance.R

