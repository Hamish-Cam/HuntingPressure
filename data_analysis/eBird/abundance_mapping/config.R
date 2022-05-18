# Config script for 'get_abundance.R'


################ Parameters ######################

# Common name of species 
species_name <- "Great Hornbill"  

# Name of the folder that contains input data folder
data_folder <- "abun-ghb"

# Selection of habitat covariates to be used for both models
# Selection is species dependent and can be informed using the IUCN RedList
habitat_covariates <- c("elevation_median",
                        "pland_02_evergreen_broadleaf", 
                        "pland_04_deciduous_broadleaf", 
                        "pland_05_mixed_forest",
                        "pland_12_cropland")

# Selection of non-habitat variables to include in the advanced model 
non_habitat_covariates <- c("hunting_mean",
                            "invasives_mean",
                            "pollution_mean",
                            "climate_mean")   

# Training proportion of train/test split
train_prop <- 0.8

# GAM degrees of freedom 
k <- 5

# GAM degrees of freedom for cyclic time of day smooth
k_time <- 7 

# Month and Day to use for abundance map prediction (day of year variable)
time_of_year <- "06-15"

# Define a threshold below which abundance values are treated as 0
zero_threshold <- 0.05

################ end ######################







