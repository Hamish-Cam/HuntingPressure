# Hamish Campbell

# Config script for 'abundance_model_comparison.R'


#### Study Species ####

# Specify the sets of hunted and non-hunted species used for analysis
# Entries must be in the form:
#   c("common name" (according to eBird), "scientific name" (according to IUCN))
hunted_species <- data.frame(c("Great Hornbill", "Buceros bicornis"),
                            c("Oriental Pied-Hornbill", "Anthracoceros albirostris"),
                            c("Scaly-breasted Partridge", "Tropicoperdix chloropus"),
                            c("White-rumped Shama",	"Kittacincla malabarica"),
                            c("Indian White-eye", "Zosterops palpebrosus"))

non_hunted_species <- data.frame(c("Banded Bay Cuckoo", "Cacomantis sonneratii"),
                                 c("Greater Racket-tailed Drongo", "Dicrurus paradiseus"),
                                 c("Large-tailed Nightjar",	"Caprimulgus macrurus"),
                                 c("Plain Flowerpecker", "Dicaeum minullum"),
                                 c("Plain-backed Sparrow", "Passer flaveolus"),
                                 c("Sooty-headed Bulbul",	"Pycnonotus aurigaster"),
                                 c("Common Tailorbird",	"Orthotomus sutorius"))

#### end ####


#### Parameters ####

# Specify whether to use the hunted or non-hunted species sets
selected_species <- non_hunted_species

# Selection of non-habitat variables to include in the advanced model 
# Pressure maps: "hunting_mean"
# Accessibility: "access_mean"
non_habitat_covariates <- c("access_mean")   

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

# Define a threshold below which abundance values are treated as 0 for prediction maps
zero_threshold <- 0.05

#### end ####



