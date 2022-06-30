# Preprocessing Overview
The preprocessing code is designed such that it can be run either locally or on the JASMIN HPC. The installation instructions given here are for running the code on JASMIN. 

Following the steps outlined here and running the preprocessing code will result in the download of all required datasets and processing as required for further analysis. In addition, various plots are generated to help the user visualise the characteristics of the data for their specific species of interest. 

The *preprocessing.R* script is the main code for preprocessing the data and should not need to be altered by the user unless functionality needs to be added/changed. The *config.R* script contains variables that the user should alter depending on the experiments they wish to carry out. 


# Data Download
For reproducing results, the user can download all required data from the [Zenodo project page](https://doi.org/10.5281/zenodo.%20%20%206761830).

For up-to-date data, the preprocessing code requires manual download of 7 datasets prior to running of the preprocessing code. Note that some of these datasets require an approved account before download is possible:
1. [eBird Basic Dataset (EBD) - Sampling event data](https://ebird.org/data/download): Required for zero-filling checklist data (account setup required)
2. [eBird Basic Dataset (EBD) - Species data](https://ebird.org/data/download) - Custom download of species data. The species of interest should be selected, but all other filters left blank (account setup required)
3. [BirdLife Range Polygons](http://datazone.birdlife.org/species/requestdis) - Range polygons for all monitored bird species (data access requested via email)
4. [NASA Landcover Data](https://lpdaac.usgs.gov/products/mcd12q1v006/) - Landcover classification data using Modis satellite imagery (1km resolution). A section of optional code in *preprocessing.R* allows download of the necessary tiles (Earthdata account setup required)
5. [EarthEnv Elevation Dataset](http://www.earthenv.org/topography) - Derived topographic continuous variables (elevation, median, GMTED2010, 1km)
6. [Harfoot et al. Threat Maps](https://www.nature.com/articles/s41559-021-01542-9) - Pressure maps of 6 major threats as defined by the IUCN Red List: agriculture, logging, hunting, invasives, pollution and climate change (whole world, 50km resolution). Email lead author to request data
7. [Accessibility Dataset](https://malariaatlas.org/research-project/accessibility-to-cities/) - Minimum travel time to urban area of >50,000 people (minutes)

# JASMIN Setup 
These steps assume that the user already: 

 - Has a JASMIN login account 
 - Can SSH into the JASMIN remote server (`ssh -A <username>@login1.jasmin.ac.uk`)
 - Knows how to transfer files onto/off of JASMIN (`rsync` recommended for these operations)
 
## Package Installation
The first step is to ensure the necessary R packages are installed on JASMIN: 
 1. Load the jaspy module: `module load jasr`
 2. Start the R command line: `R`
 3. Install the required packages available on CRAN (choose to install all packages and select a server location close to yours): `install.packages(c("sf","rnaturalearth","dplyr","raster","exactextractr","viridis","tidyverse","auk","lubridate","gridExtra","devtools","stars","ggplot2"))`
 4. (Optional Modis download) Install the MODIS package from Github: `devtools::install_github("MatMatt/MODIS")`
 5. (Optional Modis download) Test your Earthdata credentials (may need 24 hours after account creation):  `MODIS::EarthdataLogin(usr = "username", pwd = "password")
MODISoptions(check_earthdata_login = TRUE)`


## Directory Setup
The correct directory structure must be setup for correct preprocessing functionality: 
1. Place all data except eBird species checklist data into the permanent files folder (specified in config)
2. Make a project folder for the *preprocessing.R* and *config.R* files to live in
3. Within the project folder, create a data folder (specified in config)
4. Within the data folder, create a folder named *input data* (make sure to use this name) and place the eBird species checklist data inside

In summary, the directory structure should be as follows:

    perm_files directory
        -> eBird sampling data (specified in config)
        -> landcover data (.tif folder)
        -> elevation data
        -> range polygons (.gdb folder)
        -> threat map shape files
        -> accessibility data
        
    project folder
        -> data folder (specified in config)
            -> preprocessing.R
            -> config.R
            -> input data 
                -> eBird species checklist data

# Expected Output
Note that all output files will be prefaced with the species eBird short-code that the data corresponds to. After a run has been completed (see below), the following outputs will be generated and placed in the *output data* subfolder (within the data folder). 
1. GIS buffered range data: `gis-data.gpkg`
2. eBird processed output data: `ebd_processed_output.csv`
3. Landcover and elevation data at checklist locations: `covariate_checklists.csv`
4. Landcover and elevation data at all locations on prediction surface (evenly spaced points across buffered range): `covariate_predictions.csv`
5. Raster file defining all points on prediction surface geographically: `prediction-surface.tif`

In addition, various plots are generated to give the user insight into the characteristics of each of the datasets being processed. These plots are placed in the *analytics* subfolder (within the data folder). The following plots are generated: 
1. Buffered range map (each species) 
2. eBird checklist map (each species)
3. Effort covariate insights (x4 each species)
4. Urban landcover map (each species)
5. Elevation map (each species)
6. Harfoot whole world pressure map 
7. Harfoot pressure map (each species)
8. Accessibility map (each species) 


In summary, the final directory structure will be as follows: 

     perm_files directory
        -> eBird sampling data
        -> landcover data (.tif folder)
        -> elevation data
        -> range polygons (.gdb folder)
        -> threat map shape files
        -> accessibility data
            
    project folder
        -> data folder
            -> preprocessing.R
            -> config.R
            -> input data
                -> eBird checklist data
            -> output data  
                -> gis-data.gpkg
                -> ebd_processed_output.csv
                -> covariate_checklists.csv.csv
                -> covariate_predictions.csv
                -> prediction-surface.tif
            -> analytics folder 
                -> Buffered range map
                -> eBird checklist map
                -> effort covariate insights (x4)
                -> urban landcover map
                -> elevation map
                -> Harfoot whole world pressure map
                -> Harfoot pressure map
                -> accessibility map

# Completing a Run
Once all of the above steps have been completed, the user is ready to start a run. JASMIN users the SLURM system for scheduling jobs - once a job is complete the outputs described above should be available for transfer. An example workflow for submitting a job is given here:
1. Write a Slurm submission script (example given below)
2. Submit the job: `sbatch <slurm file>`
3. Check the status of the job: `sacct`
4. Check the statistics of the completed run: `seff <job ID>`
5. Check the *analytics* and *output data* folders for the relevant files and transfer them (either within JASMIN or to your local machine using `rsync`)

## Example Slurm job script

    #!/bin/bash 
    #SBATCH --partition=short-serial 
    #SBATCH -o %j.out 
    #SBATCH -e %j.err 
    #SBATCH --time=15:00:00 
    #SBATCH --mem=0 
    
    # Must add R to workspace 
    module add jasr 
    
    # Go to relevant folder 
    cd HuntingPressure/1.preprocessing
    
    # Execute the job in question 
    Rscript preprocessing.R

