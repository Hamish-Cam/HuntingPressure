# Config script for 'abundance_model_comparison.R'


################ Parameters ######################


# Specify which species to study - Entries must be in the form:
#     c(common name, scientific name used by IUCN)
selected_species <- data.frame(c("Great Hornbill", "Buceros bicornis"),
                               c("Large Cuckooshrike", "Coracina javensis"))

# Selection of non-habitat variables to include in the advanced model 
non_habitat_covariates <- c("hunting_mean")   

# Name of the folder that contains input data folder
data_folder <- "data"

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



