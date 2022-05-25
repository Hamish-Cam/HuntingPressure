# Config script for 'preprocessing.R'


################ Parameters ######################


# Specify which species to study - Entries must be in the form:
#     c(common name (according to eBird), scientific name (according to IUCN))
selected_species <- data.frame(c("Great Hornbill", "Buceros bicornis"),
                               c("Oriental Pied-Hornbill", "Anthracoceros albirostris"),
                               c("Wreathed Hornbill", "Rhyticeros undulatus"),
                               c("Scaly-breasted Partridge", "Tropicoperdix chloropus"),
                               c("Green Imperial-Pigeon", "Ducula aenea"),
                               c("Ashy-headed Green-Pigeon", "Treron phayrei"),
                               c("Alexandrine Parakeet", "Palaeornis eupatria"),
                               c("Blossom-headed Parakeet", "Himalayapsitta roseata"),
                               c("White-rumped Shama",	"Kittacincla malabarica"),
                               c("Indian White-eye", "Zosterops palpebrosus"))

# Country of study (within SE Asia)
country <- "Thailand"

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
