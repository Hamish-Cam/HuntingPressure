# Hamish Campbell

# Calculation of presence proportions for selected provisional study species and 
# Thailand landcover types

library(lubridate)
library(sf)
library(raster)
library(dggridR)
library(pdp)
library(mgcv)
library(fitdistrplus)
library(viridis)
library(fields)
library(tidyverse)
library(auk)

# Select which packages we want to take duplicated function names from
select <- dplyr::select
map <- purrr::map
projection <- raster::projection

# Load the variable values from the config file
source("../1.preprocessing/config.R")


#### Species Proportions ####

# Transform species selection from config into a dataframe with correct form
selected_species <- data.frame(t(selected_species))
colnames(selected_species) <- c("common_name", "scientific_name_IUCN")

# Get eBird taxonomy for the selected species 
species_list <- get_ebird_taxonomy()
requested_species <- species_list[is.element(species_list$common_name, selected_species$common_name),]

# Extract the eBird scientific name (may differ to IUCN) and species code for requested species
species_eBird_data <- select(requested_species, common_name, scientific_name, species_code) %>%
  rename(scientific_name_eBird = scientific_name,
         species_code_eBird = species_code)

# Merge all of the taxonomic data we need for the selected species
species_data <- merge(selected_species, species_eBird_data)

# Make an empty table for containing positive sighting proportions
sighting_props <- data.frame(species=c(), pos_sighting_prop=c())

# Complete steps for each requested species 
for (row in 1:nrow(species_data)){
  current_species <- species_data[row, ]
  
  # Get the short code for the current species 
  short_code <- current_species$species_code
  
  # Load eBird data 
  ebird_data <- read_csv(file.path("../1.preprocessing", data_folder, "output data", 
                                      sprintf("%s_ebd_processed_output.csv", short_code))) %>% 
    # Make the 'protocol type' column categorical with 2 options 
    mutate(protocol_type = factor(protocol_type, 
                                  levels = c("Stationary" , "Traveling"))) %>%
    # remove observations with no count
    filter(!is.na(observation_count))
  
  # Calculate the proportion of positive sightings seen 
  total_checklists <- nrow(ebird_data)
  pos_sightings <- nrow(ebird_data[ebird_data$observation_count > 0,])
  prop <- (pos_sightings/total_checklists) * 100
  
  # Add the proportion to the data table 
  sighting_props <- rbind(sighting_props, c(current_species$common_name, prop))
}

# Rename the column titles
colnames(sighting_props) <- c("Species", "Proportion of Positive Sightings (%)")

# Save the resulting table after ordering from highest to lowest
sighting_props <- sighting_props[order(sighting_props$"Proportion of Positive Sightings (%)", decreasing = TRUE),]
write.csv(sighting_props, "sighting_props.csv", row.names = FALSE)

#### end ####


#### Landcover Proportions ####

# Load landcover data 
landcover_data <- read_csv(file.path("../1.preprocessing", data_folder, "output data", 
                                      sprintf("%s_covariate_predictions.csv", short_code)))

# Make a list of all possible landcover types
landcover_types <- c("pland_00_water", 
                     "pland_01_evergreen_needleleaf", 
                     "pland_02_evergreen_broadleaf", 
                     "pland_03_deciduous_needleleaf", 
                     "pland_04_deciduous_broadleaf", 
                     "pland_05_mixed_forest",
                     "pland_06_closed_shrubland", 
                     "pland_07_open_shrubland", 
                     "pland_08_woody_savanna", 
                     "pland_09_savanna", 
                     "pland_10_grassland", 
                     "pland_11_wetland", 
                     "pland_12_cropland", 
                     "pland_13_urban", 
                     "pland_14_mosiac", 
                     "pland_15_barren")

# Make an empty table for containing occurence proportions
occurence_props <- data.frame(type=c(), prop=c())

# Find proportion for each type
for (type in landcover_types){
  if(type %in% colnames(landcover_data)){
    # Find the total proportion within Thailand
    total_prop <- sum(landcover_data[,type])
    
    # Max proportion would be 1 for every 5km neighbourhood
    max_prop <- nrow(landcover_data)
    
    # Get the occurence proportion
    occurence_prop <- (total_prop/max_prop) * 100
  }
  # If the landcover type doesn't exist, that means the proportion is zero
  else{
    occurence_prop <- 0.0
  }
  
  # Add the proportion to the table
  occurence_props <- rbind(occurence_props, c(type, occurence_prop))
}

# Rename the column titles
colnames(occurence_props) <- c("Landcover Type", "Occurence Proportion (%)")

# Save the resulting table after ordering from highest to lowest
occurence_props$"Occurence Proportion (%)" <- as.numeric(as.character(occurence_props$"Occurence Proportion (%)"))
occurence_props <- occurence_props[order(occurence_props$"Occurence Proportion (%)", decreasing = TRUE),]
write.csv(occurence_props, "landcover_props.csv", row.names = FALSE)

#### end ####

