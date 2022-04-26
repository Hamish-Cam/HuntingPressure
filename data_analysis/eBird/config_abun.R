# Config script for 'get_abundance.R'


################ Parameters ######################

# Common name and scientific name of species 
species_name <- "Red Grouse"            

# Relevant habitat covariates in the form: pland_<landcover code>_<landcover name> 
# Selection is species dependent and can be informed using the IUCN RedList
habitat_covariates <- c("pland_02_evergreen_broadleaf", 
                        "pland_04_deciduous_broadleaf", 
                        "pland_05_mixed_forest",
                        "pland_12_cropland")

# Name of the folder that contains input data folder
data_folder <- "abun-ghb"

# Train proportion used for modelling 
train_prop <- 0.8

#GAM
# degrees of freedom 
k <- 5

# degrees of freedom for cyclic time of day smooth
k_time <- 7 

# Month and Day to use for abundance prediction (day of year)
time_of_year <- "06-15"

# Define a threshold below which abundance values are treated as 0
zero_threshold <- 0.05

################ end ######################







