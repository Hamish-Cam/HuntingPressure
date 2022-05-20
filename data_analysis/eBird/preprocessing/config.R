# Config script for 'preprocessing.R'


################ Parameters ######################


# Specify which species to study - Entries must be in the form:
#     c(common name, scientific name used by IUCN)
selected_species <- data.frame(c("Great Hornbill", "Buceros bicornis"),
                               c("Large Cuckooshrike", "Coracina javensis")) 

# Continent of study 
continent <- "Asia"

# Earthdata login credentials
username = "hrac2"
password = "mQ5taFmenGs9oHS"

# Name of the folder that contains input data folder for species specific data
data_folder <- "data"

# Path to permanent location containing general data (sampling, elevation, range etc.)
perm_files_location <- "~/hunting_data"

# Name of downloaded sampling data file 
sampling_data <- "ebd_sampling_relApr-2022.txt"


################ end ######################
