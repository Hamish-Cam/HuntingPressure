# Hamish Campbell
# Various useful statistics for hunted/non-hunted data 


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

# Select which packages want to take duplicated function names from
select <- dplyr::select
map <- purrr::map
projection <- raster::projection

# Load the variable values from the config file
source("config.R")


#### Taxonomic Data ####

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

#### end ####


# Make an empty table for containing positive sighting proportions
sighting_props <- data.frame(species=c(), pos_sighting_prop=c())

# Complete steps for each requested species 
for (row in 1:nrow(species_data)){
  current_species <- species_data[row, ]
  
  # Get the short code for the current species 
  short_code <- current_species$species_code
  
  
  #### Data Prep/Loading #### 
  
  # Load eBird data 
  ebird_data <- read_csv(file.path(data_folder, "input data", sprintf("%s_ebd_processed_output.csv", short_code))) %>% 
    # Make the 'protocol type' column categorical with 2 options 
    mutate(protocol_type = factor(protocol_type, 
                                  levels = c("Stationary" , "Traveling"))) %>%
    # remove observations with no count
    filter(!is.na(observation_count))
  
  #### end ####
  
  
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






