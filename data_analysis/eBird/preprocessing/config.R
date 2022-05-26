# Config script for 'preprocessing.R'


################ Species Sets ######################

# Specify the sets of hunted and non-hunted species used for analysis
# Entries must be in the form:
#   c("common name" (according to eBird), "scientific name" (according to IUCN))
hunted_species <- data.frame(c("Great Hornbill", "Buceros bicornis"),
                            c("Oriental Pied-Hornbill", "Anthracoceros albirostris"),
                            c("Wreathed Hornbill", "Rhyticeros undulatus"),
                            c("Scaly-breasted Partridge", "Tropicoperdix chloropus"),
                            c("Green Imperial-Pigeon", "Ducula aenea"),
                            #c("Ashy-headed Green-Pigeon", "Treron phayrei"),
                            c("Alexandrine Parakeet", "Palaeornis eupatria"),
                            c("Blossom-headed Parakeet", "Himalayapsitta roseata"),
                            c("White-rumped Shama",	"Kittacincla malabarica"),
                            c("Indian White-eye", "Zosterops palpebrosus"))

non_hunted_species <- data.frame(c("Banded Bay Cuckoo", "Cacomantis sonneratii"),
                                 c("Large Cuckooshrike", "Coracina javensis"),
                                 c("Greater Racket-tailed Drongo", "Dicrurus paradiseus"),
                                 c("Blyth's Frogmouth", "Batrachostomus javensis"),
                                 c("Large-tailed Nightjar",	"Caprimulgus macrurus"),
                                 c("Rufous Woodpecker",	"Micropternus brachyurus"),
                                 c("Plain Flowerpecker", "Dicaeum minullum"),
                                 c("Plain-backed Sparrow", "Passer flaveolus"),
                                 c("Sooty-headed Bulbul",	"Pycnonotus aurigaster"),
                                 c("Common Tailorbird",	"Orthotomus sutorius"))

################ end ######################


################ Parameters ######################

# Specify whether to use the hunted or non-hunted species sets
selected_species <- hunted_species

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
