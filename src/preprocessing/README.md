

# Preprocessing Overview
The preprocessing code is written in the R language and is designed such that it can be run either locally or on the JASMIN HPC. The installation instructions given here are for running the code on JASMIN. 

Following the steps outlined here and running the preprocessing code will result in all required datasets being downloaded and processed as required for further analysis. In addition, various plots are generated to help the user visualise the characteristics for their specific species of interest. 

The *preprocessing.R* script is the main code for preprocessing the data and should not need to be altered by the user. The *config.R* script contains variables that the user should alter depending on the experiments they wish to carry out. 


# Data Download
The preprocessing code requires manual download of 3 datasets prior to running the program. Two of these datasets are from [eBird](https://ebird.org/home) and require an approved account to access the data. Once an account has been approved, follow the 3 steps outlined below:
1. [eBird Basic Dataset (EBD)](https://ebird.org/data/download) - Sampling event data (required for obtaining presence/absence data)
2. [eBird Basic Dataset (EBD)](https://ebird.org/data/download) - Checklist data (custom download). The species of interest should be selected, but all other filters left blank
3. [EarthEnv Elevation Dataset](http://www.earthenv.org/topography) - Derived topographic continuous variables (elevation, median, GMTED2010, 1km)

# JASMIN Setup 
These steps assume that the user already: 

 - Has a JASMIN login account 
 - Can SSH into the JASMIN remote server (`ssh -A <username>@login1.jasmin.ac.uk`)
 - Knows how to transfer files onto/off of JASMIN (`rsync` recommended for these operations)
 
## Package Installation
The first step is to ensure the necessary R packages are installed on JASMIN: 
 1. Load the jaspy module: `module load jasr`
 2. Start the R command line: `R`
 3. Install the required packages available on CRAN (choose to install all packages and select a server location close to yours): `install.packages(c("sf","rnaturalearth","dplyr","raster","exactextractr","viridis","tidyverse","auk","lubridate","gridExtra","devtools"))`
 4. Install the MODIS package from Github: `devtools::install_github("MatMatt/MODIS")`
 5. Make a [NASA Earthdata](https://urs.earthdata.nasa.gov/home) account for accessing the MODIS data 
 6. Test your Earthdata credentials (may need 24 hours after account creation):  `MODIS::EarthdataLogin(usr = "username", pwd = "password")
MODISoptions(check_earthdata_login = TRUE)`


## Directory Setup
Once the relevant datasets have been downloaded (info above), along with the *preprocessing.R*/*config.R* scripts, and transferred onto the JASMIN servers, the correct directory structure must be setup (naming convention must match those stated below): 
1. Place the eBird sampling data and elevation data into your home directory (use `cd ~` to find where this is)
2. Make a project folder where for the *preprocessing.R* and *config.R* files to live in
3. Within the project folder, create a *data-SpeciesName* folder
4. Within the data folder, create a folder named *input data* and place the eBird checklist data within this folder

In summary, the directory structure should be as follows:

    home directory
        -> eBird sampling data
        -> elevation data
        
    project folder
        -> data-<species name> folder
            -> input data folder
                -> eBird checklist data

# Expected Output
After a run has been completed (see below), the following outputs will be generated and placed in the *output data* subfolder: 
1. GIS area of interest (AOI data: `gis-data.gpkg`
2. eBird processed output data: `ebd_processed_output.csv`
3. Landcover and elevation data at checklist locations: `landcover_and_elevation_checklists.csv`
4. Landcover and elevation data at all locations on prediction surface (evenly spaced points across AOI): `landcover_and_elevation_prediction.csv`
5. Raster file defining all points on prediction surface geographically: `prediction-surface.tif`

In addition, various plots are generated to give the user insight into the characteristics of each of the datasets being processed. These plots are placed in the *analytics* folder. The following plots are generated: 
1. AOI map 
2. eBird checklist map
3. Effort covariate insights (x4)
4. Urban landcover map
5. Elevation map  

In summary, the final directory structure will be as follows: 

     home directory
            -> eBird sampling data
            -> elevation data
            
    project folder
        -> data-<species name> folder
            -> input data folder
                -> eBird checklist data
            -> output data folder 
                -> gis-data.gpkg
                -> ebd_processed_output.csv
                -> landcover_and_elevation_checklists.csv
                -> landcover_and_elevation_prediction.csv
                -> prediction-surface.tif
            -> analytics folder 
                -> AOI map
                -> eBird checklist map
                -> effort covariate insights (x4)
                -> urban landcover map
                -> elevation map




# Completing a Run
Once all of the above steps have been completed, the user is ready to start a run. JASMIN users the SLURM system for scheduling jobs - once a job is complete the outputs described above should be available for transfer. An example workflow for submitting a job is given here:
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
    #SBATCH --time=03:00:00 
    #SBATCH --mem=3G 
    
    # Must add R to workspace 
    module add jasr 
    
    # Go to relevant folder 
    cd HuntingPressure/data_analysis/eBird 
    
    # Execute the job in question 
    Rscript preprocessing.R

