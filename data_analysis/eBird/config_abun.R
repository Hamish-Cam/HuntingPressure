# Config script for 'get_abundance.R'


################ Parameters ######################

# Common name of species 
species_name <- "Red Grouse"  

# Name of the folder that contains input data folder
data_folder <- "abun-grouse"

# Relevant habitat covariates in the form: pland_<landcover code>_<landcover name> 
# Selection is species dependent and can be informed using the IUCN RedList
habitat_covariates <- c("pland_05_mixed_forest",
                        "pland_07_open_shrubland", 
                        "pland_10_grassland")

# GAM degrees of freedom 
k <- 5

# GAM degrees of freedom for cyclic time of day smooth
k_time <- 7 

# Month and Day to use for abundance map prediction (day of year variable)
time_of_year <- "06-15"

# Define a threshold below which abundance values are treated as 0
zero_threshold <- 0.05

################ end ######################







