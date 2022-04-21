# Config script for 'preprocessing.R'


################ Parameters ######################

# Common name and scientific name of species 
species_name <- "Red Grouse"            
species_name_scientific <- "Lagopus lagopus"

# Choice of country that you want to analyse using ISO 2-letter country code: 
# https://en.wikipedia.org/wiki/ISO_3166-1_alpha-2 
# NOTE: To specify the whole world, just pass the empty character
country_choice <- "GB"

# Name of the folder that contains input data folder
data_folder <- "data-grouse"

# Name of downloaded species data file
species_data <- "ebd_wilpta1_relMar-2022.txt"

# Path to permanent location containing sampling and elevation data
perm_files_location = "~"

# Name of downloaded sampling data file 
sampling_data <- "ebd_sampling_relFeb-2022.txt"

################ end ######################