# Config script for 'abundance_model_comparison.R'


################ Species Sets ######################

# Specify the sets of hunted and non-hunted species used for analysis
# Entries must be in the form:
#   c("common name" (according to eBird), "scientific name" (according to IUCN))
hunted_species <- data.frame(c("Great Hornbill", "Buceros bicornis"),
                            c("Oriental Pied-Hornbill", "Anthracoceros albirostris"),
                            c("Wreathed Hornbill", "Rhyticeros undulatus"),
                            c("Scaly-breasted Partridge", "Tropicoperdix chloropus"),
                            c("Green Imperial-Pigeon", "Ducula aenea"),
                            c("Ashy-headed Green-Pigeon", "Treron phayrei"),
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



