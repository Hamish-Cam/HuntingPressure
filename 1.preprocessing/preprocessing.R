# Hamish Campbell 

#> Code based heavily on:
#> https://cornelllabofornithology.github.io/ebird-best-practices/index.html

#> NOTE: the readme should be consulted before running this code.

library(sf)
library(rnaturalearth)
library(dplyr)
library(raster)
library(MODIS) 
library(exactextractr)
library(viridis)
library(tidyverse)
library(auk)
library(lubridate)
library(devtools)
library(gridExtra)
library(stars)
library(ggplot2)

# resolve namespace conflicts
select <- dplyr::select
map <- purrr::map
projection <- raster::projection

# Load the variable values from the config file
source("config.R")

# Check that the relevant output directories exist
if (!dir.exists(file.path(data_folder, "output data"))){
  dir.create(file.path(data_folder, "output data"))
}

if (!dir.exists(file.path(data_folder, "analytics"))){
  dir.create(file.path(data_folder, "analytics"))
}


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


# Complete preprocessing steps for each requested species 
for (row in 1:nrow(species_data)){
  current_species <- species_data[row, ]
  
  # Get the name of the species data file corresponding to the current species
  short_code <- current_species$species_code
  data_timestamp <- str_extract(sampling_data, "\\w{3}-\\d{4}")
  current_species_data <- sprintf("ebd_%s_rel%s.txt", short_code, data_timestamp)

  #### Species Range ####
  
  # File to save GIS data 
  f_ne <- file.path(data_folder, "output data", sprintf("%s_gis-data.gpkg", short_code))
  
  # Delete any existing GIS data
  unlink(f_ne)
  
  # Get species range polygon from BirdLife data for species requested in config
  range_query <- paste('SELECT * FROM "All_Species" WHERE sci_name = \'', 
                       current_species$scientific_name_IUCN, '\'', sep='')
  species_range <- st_read(dsn=file.path(perm_files_location, "BOTW.gdb"), query = range_query)
  
  # Get range polygon corresponding to extant, native and resident row
  species_range <- filter(species_range, presence==1, origin==1, seasonal==1) %>%
                      select(Shape)
  
  # Some of the range polygons have issues 
  if(!st_is_valid(species_range)){
    species_range <- st_make_valid(
      species_range,
      oriented = FALSE,
      s2_options = s2::s2_options(snap = s2::s2_snap_precision(1e+06))
    )
  }
  
  # Add a buffer of 100km around the inexact range polygon to allow for error
  species_range <- st_buffer(species_range, 100000)
  
  # Limit this extended range to that only within the country of study
  country_boundary <- ne_countries(country = country, scale = 50) %>%
            polygons(.) %>%
            st_as_sf(.)
  species_range <- st_intersection(species_range, country_boundary) 

  # Now get all Asian land and country boundaries for plotting purposes 
  # Land boundaries 
  ne_land <- ne_download(scale = 50, category = "cultural",
                         type = "admin_0_countries_lakes",
                         returnclass = "sf") %>%
    filter(CONTINENT == "Asia") %>%
    st_set_precision(1e6) %>%
    st_union()
  
  # Country boundaries (limit to continent by intersect command)
  ne_country_lines <- ne_download(scale = 50, category = "cultural",
                                  type = "admin_0_boundary_lines_land",
                                  returnclass = "sf") %>% 
    st_geometry()
  
  ne_country_lines <- st_intersects(ne_country_lines, ne_land, sparse = FALSE) %>%
    as.logical() %>%
    {ne_country_lines[.]}
  
  # Output relevant GIS data
  write_sf(ne_land, f_ne, "ne_land")
  write_sf(ne_country_lines, f_ne, "ne_country_lines")
  write_sf(species_range, f_ne, "range")
  
  #### end ####
  
  #### Species Range Analysis ####
  
  # Plot the land, country boundaries and range, using range bounds to limit axes
  pdf(file = file.path(data_folder, "analytics", sprintf("%s_range_map.pdf", short_code)))
  plot(ne_land, axes=TRUE, xlim = st_bbox(species_range)[c(1,3)], 
          ylim = st_bbox(species_range)[c(2,4)], main=paste(current_species$common_name, "Range in", country))
  plot(species_range, add=TRUE, col=alpha("orange",0.7), border='transparent')
  plot(ne_country_lines, add=TRUE)
  dev.off()
  
  #### end ####
  
  
  #### eBird Preprocess ####
  
  # Setup the eBird object for data manipulation
  ebd <- auk_ebd(file.path(getwd(), data_folder, "input data", current_species_data), 
                 file_sampling = file.path(perm_files_location, sampling_data))
  
  # Define/print the filters to use for the eBird data 
  # (https://cornelllabofornithology.github.io/auk/reference/index.html#section-filter)
  ebd_filters <- ebd %>% 
    # Restrict checklists to those within smallest bounding box containing range polygons
    auk_bbox(species_range) %>%
    # Restrict species to just that of interest
    auk_species(current_species$scientific_name_eBird) %>%
    # Restrict checklists to just the country that the user has specified 
    auk_country(country) %>%
    # restrict to the standard traveling and stationary count protocols
    auk_protocol(protocol = c("Stationary", "Traveling")) %>% 
    auk_complete()
  ebd_filters
  
  # Define output file names
  f_ebd <- file.path(data_folder, "ebd_filtered_output.txt")
  f_sampling <- file.path(data_folder, "ebd_samp_filtered_output.txt")
  
  # Apply filters to eBird data
  auk_filter(ebd_filters, file = f_ebd, file_sampling = f_sampling, overwrite = TRUE)
  
  # Read in the data and zero-fill using sampling info to produce full presence/absence data
  ebd_zf <- auk_zerofill(f_ebd, f_sampling, collapse = TRUE)
  
  # Needed to keep only the eBird checklists within the range polygon
  # Convert to sf object first
  ebd_zf_sf <- ebd_zf %>% 
    select(longitude, latitude) %>% 
    st_as_sf( coords = c("longitude", "latitude"), crs = 4326)
  
  # Ensure polygon and eBird data have the same crs
  poly_crop <- st_transform(species_range, crs = st_crs(ebd_zf_sf))
  
  # Identify checklists within range polygon
  in_range <- st_within(ebd_zf_sf, poly_crop, sparse = FALSE)
  
  # Subset the data to get only checklists in range
  ebd_zf <- ebd_zf[in_range[, 1], ]
  
  # Function to convert 'time observation' to hours (0-24)
  time_to_decimal <- function(x) {
    x <- hms(x, quiet = TRUE)
    hour(x) + minute(x) / 60 + second(x) / 3600
  }
  
  # Clean up data slightly
  ebd_zf <- ebd_zf %>% 
    mutate(
      # convert 'X' (species count not tracked)  to NA
      observation_count = if_else(observation_count == "X", 
                                  NA_character_, observation_count),
      observation_count = as.integer(observation_count),
      # convert effort_distance_km to 0 for non-travelling counts
      effort_distance_km = if_else(protocol_type != "Traveling", 
                                   0, effort_distance_km),
      # convert time to decimal hours using our function
      time_observations_started = time_to_decimal(time_observations_started),
      # split date into year and day of year
      year = year(observation_date),
      day_of_year = yday(observation_date)
    )
  
  # Additional filtering of checklists to remove some outliers with very high 
  # effort which may skew results 
  ebd_zf_filtered <- ebd_zf %>% 
    filter(
      # effort filters (< 5 hours and < 5km)
      duration_minutes <= 5 * 60,
      effort_distance_km <= 5,
      # last 10 years of data only
      year >= 2010,
      # 10 or fewer observers
      number_observers <= 10)
  
  # Remove some redundant columns not being used for subsequent analysis
  ebird <- ebd_zf_filtered %>% 
    select(checklist_id, observer_id, sampling_event_identifier,
           scientific_name,
           observation_count, species_observed, 
           state_code, locality_id, latitude, longitude,
           protocol_type, all_species_reported,
           observation_date, year, day_of_year,
           time_observations_started, 
           duration_minutes, effort_distance_km,
           number_observers)
  
  # Note: eBird analysis happens post covariate processing due to alterations
  
  #### end ####
  
  
  #### (Optional) MODIS Landcover Download ####
  
  # Source: https://lpdaac.usgs.gov/products/mcd12q1v006/
  # Run this section of code to get up to date landcover data or tiles for a region != Thailand.
  # Note: you will need an Earthdata login
  
  # # Read the GIS boundary chosen in the previous section
  # gis_land <- read_sf(f_ne, "range") %>%
  #   # project to the native modis projection
  #   st_transform(crs = paste("+proj=sinu +lon_0=0 +x_0=0 +y_0=0",
  #                            "+a=6371007.181 +b=6371007.181 +units=m +no_defs"))
  # 
  # # Get (and print) list of MODIS tiles required to cover this region
  # tiles <- getTile(gis_land)
  # tiles@tile
  # 
  # # Earliest/latest year of ebird data available
  # begin_year <- format(min(ebird$observation_date), "%Y.01.01")
  # end_year <- format(max(ebird$observation_date), "%Y.12.31")
  # 
  # # Store credentials for accessing MODIS data and test connection
  # MODIS::EarthdataLogin(usr = "username", pwd = "password")
  # MODISoptions(check_earthdata_login = TRUE)
  # 
  # # Download tiles for area of interest and combine into a single raster for each year
  # # then save into subfolder 'modis'
  # tifs <- runGdal(product = "MCD12Q1", collection = "006", SDSstring = "01",
  #                 extent = gis_land %>% st_buffer(dist = 10000),
  #                 begin = begin_year, end = end_year,
  #                 outDirPath = data_folder, job = "modis",
  #                 MODISserverOrder = "LPDAAC") %>%
  #   pluck("MCD12Q1.006") %>%
  #   unlist()
  # 
  # # Rename tifs to have more descriptive names
  # new_names <- format(as.Date(names(tifs)), "%Y") %>%
  #   sprintf("modis_mcd12q1_umd_%s.tif", .) %>%
  #   file.path(dirname(tifs), .)
  # file.rename(tifs, new_names)
  # 
  # # Now move modis tiles to permanent files location
  # dir.create(file.path(perm_files_location, "modis"))
  # file.copy(from=file.path(data_folder,"modis"), to=perm_files_location, overwrite = TRUE, recursive = TRUE)
  # unlink(file.path(data_folder,"modis"), recursive=TRUE)
  
  #### end ####
  
  #### MODIS Landcover Data ####
  
  # Read the GIS boundary chosen in the previous sections
  gis_land <- read_sf(f_ne, "range") %>% 
    # project to the native modis projection
    st_transform(crs = paste("+proj=sinu +lon_0=0 +x_0=0 +y_0=0",
                             "+a=6371007.181 +b=6371007.181 +units=m +no_defs"))
  
  # Load the landcover data into R as a single object, each year as a layer
  landcover <- list.files(file.path(perm_files_location, "modis"), "^modis_mcd12q1_umd", 
                          full.names = TRUE) %>% 
    stack()
  
  # Label the various layers with the corresponding year, print the object to inspect
  landcover <- names(landcover) %>% 
    str_extract("(?<=modis_mcd12q1_umd_)[0-9]{4}") %>% 
    paste0("y", .) %>% 
    setNames(landcover, .)
  landcover
  
  # Extract the most recent year of MODIS we have (will be less than eBird data)
  max_lc_year <- names(landcover) %>% 
    str_extract("[0-9]{4}") %>% 
    as.integer() %>% 
    max()
  
  # Checklist locations aren't exact and species interact with their habitat on 
  # a greater scale, thus, we want to summarise a species' habitat by using a 
  # neighborhood. Specifically, we use a 2.5km x 2.5km box around the checklist
  # location and take the proportions of landcover types found within this box
  # to describe the habitat in this area. We choose this size since we limited max
  # travel to 5km and size makes sense for scale of bird species
  
  # Get neighborhood radius using resolution of MODIS
  neighborhood_radius <- 5 * ceiling(max(res(landcover))) / 2
  
  # Get landcover neighborhood data for each checklist point
  ebird_buff <- ebird %>% 
    # Consider distinct checklists
    distinct(year = format(observation_date, "%Y"),
             locality_id, latitude, longitude) %>% 
    # For years not available, use max we have for landcover data. Add lc year being used
    mutate(year_lc = if_else(as.integer(year) > max_lc_year, 
                             as.character(max_lc_year), year),
           year_lc = paste0("y", year_lc)) %>% 
    # Convert location to spatial features
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
    # Transform to modis projection
    st_transform(crs = projection(landcover)) %>% 
    # Buffer to create neighborhood around each point
    st_buffer(dist = neighborhood_radius) %>% 
    # Nest by year to match of with lc data
    nest(data = c(year, locality_id, geometry))
  
  # Function to summarize landcover data for all checklists in a given year
  calculate_pland <- function(yr, regions, lc) {
    locs <- st_set_geometry(regions, NULL)
    exact_extract(lc[[yr]], regions, progress = FALSE) %>% 
      map(~ count(., landcover = value)) %>% 
      tibble(locs, data = .) %>% 
      unnest(data)
  }
  
  # Iterate over all years, extracting landcover summaries for each checklists
  lc_extract <- ebird_buff %>% 
    mutate(pland = map2(year_lc, data, calculate_pland, lc = landcover)) %>% 
    select(pland) %>% 
    unnest(cols = pland)
  
  # Convert landcover summaries into proportions
  pland <- lc_extract %>% 
    group_by(locality_id, year) %>% 
    mutate(pland = n / sum(n)) %>% 
    ungroup() %>% 
    select(-n) %>% 
    # remove NAs after tallying so pland is relative to total number of cells
    filter(!is.na(landcover))
  
  # Convert landcover names to be more descriptive
  lc_names <- tibble(landcover = 0:15,
                     lc_name = c("pland_00_water", 
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
                                 "pland_15_barren"))
  pland <- pland %>% 
    inner_join(lc_names, by = "landcover") %>% 
    arrange(landcover) %>% 
    select(-landcover)
  
  # Make each checklist a row, filling in implicit missing values with 0s 
  pland <- pland %>% 
    pivot_wider(names_from = lc_name, 
                values_from = pland, 
                values_fill = list(pland = 0))
  
  
  # Now we need to make the prediction surface using the landcover data. For this 
  # we need a regular grid with neighbourhoods the same size as we used before. For
  # each cell in this grid we also want to saves its PLAND values so we have the same 
  # metrics available for prediction as we had when the GAM model was fit 
  
  # Make an altered raster grid using the landcover data as a basis
  agg_factor <- round(2 * neighborhood_radius / res(landcover))
  r <- raster(landcover) %>% 
    aggregate(agg_factor) 
  
  # Assign any values inside GIS border a value of 1 and make those outside 0 
  r <- gis_land %>% 
    st_transform(crs = projection(r)) %>% 
    rasterize(r, field = 1) %>% 
    # remove any empty cells at edges
    trim()
  
  # Save the plain prediction surface raster
  r <- writeRaster(r, filename = file.path(data_folder, "output data", 
                      sprintf("%s_prediction-surface.tif", short_code)), overwrite = TRUE)
  
  # Get cell centers and create neighborhoods
  r_centers <- rasterToPoints(r, spatial = TRUE) %>% 
    st_as_sf() %>% 
    transmute(id = row_number())
  r_cells <- st_buffer(r_centers, dist = neighborhood_radius)
  
  # Extract landcover summaries within neighborhoods, only needed most recent year
  # (uses same process as used previously above)
  lc_extract_pred <- landcover[[paste0("y", max_lc_year)]] %>% 
    exact_extract(r_cells, progress = FALSE) %>% 
    map(~ count(., landcover = value)) %>% 
    tibble(id = r_cells$id, data = .) %>% 
    unnest(data)
  
  # Calculate the proportions for each landcover class
  pland_pred <- lc_extract_pred %>% 
    count(id, landcover) %>% 
    group_by(id) %>% 
    mutate(pland = n / sum(n)) %>% 
    ungroup() %>% 
    select(-n) %>% 
    filter(!is.na(landcover))
  
  # Convert names to be more descriptive (same names as before)
  pland_pred <- pland_pred %>% 
    inner_join(lc_names, by = "landcover") %>% 
    arrange(landcover) %>% 
    select(-landcover)
  
  # Transform to wide format, filling in implicit missing values with 0s
  pland_pred <- pland_pred %>% 
    pivot_wider(names_from = lc_name, 
                values_from = pland, 
                values_fill = list(pland = 0)) %>% 
    mutate(year = max_lc_year) %>% 
    select(id, year, everything())
  
  # Add lat/lon to landcover data for each checklist
  pland_coords <- st_transform(r_centers, crs = 4326) %>% 
    st_coordinates() %>% 
    as.data.frame() %>% 
    cbind(id = r_centers$id, .) %>% 
    rename(longitude = X, latitude = Y) %>% 
    inner_join(pland_pred, by = "id")
  
  #### end ####
  
  #### MODIS Landcover Analysis ####
  
  # Plot urban coverage as a landcover map test
  urban_cover <- pland_coords %>% 
    # convert to spatial features and project onto prediction surface
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
    st_transform(crs = projection(r)) %>% 
    # rasterize points
    rasterize(r, field = "pland_13_urban") %>%
    # project to correct crs
    projectRaster(crs = st_crs(4326)$proj4string, method = "ngb") %>%
    # trim off empty edges of raster
    trim()
  
  # make a map
  pdf(file = file.path(data_folder, "analytics", sprintf("%s_urban_coverage_map.pdf", short_code)))
  par(mar = c(3,3,2,5))
  t <- str_glue("Urban Coverage",", MODIS {max_lc_year}")
  plot(ne_land, axes=TRUE, xlim = st_bbox(species_range)[c(1,3)], 
          ylim = st_bbox(species_range)[c(2,4)], main=t)
  plot(species_range, add=TRUE, border='black')
  plot(urban_cover, add = TRUE, col = viridis(10))
  plot(ne_country_lines, add=TRUE)
  dev.off()
  
  #### end ####
  
  
  #### Elevation Data ####
  
  # Data download from: http://www.earthenv.org/topography
  elev <- raster(file.path(perm_files_location, "elevation_1KMmd_GMTEDmd.tif"))
  
  # crop and buffer prediction region by 10 km to provide a little wiggle room
  elev <- gis_land %>% 
    st_buffer(dist = 10000) %>% 
    st_transform(crs = projection(elev)) %>% 
    crop(elev, .) %>% 
    projectRaster(crs = projection(landcover))
  
  # Create neighbourhood area around each checklist location
  ebird_buff_noyear <- ebird %>% 
    distinct(locality_id, latitude, longitude) %>% 
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
    st_transform(crs = projection(elev)) %>% 
    st_buffer(dist = neighborhood_radius)
  
  # Extract elevation values within neighbourhoods to calculate median and standard 
  # deviation of elevation in that neighbourhood to summarize it 
  locs <- st_set_geometry(ebird_buff_noyear, NULL) %>% 
    mutate(id = row_number())
  elev_checklists <- exact_extract(elev, ebird_buff_noyear, progress = FALSE) %>% 
    map_dfr(~ tibble(elevation_median = mean(.$value, na.rm = TRUE),
                     elevation_sd = sd(.$value, na.rm = TRUE))) %>% 
    # join to lookup table to get locality_id
    bind_cols(locs, .)
  
  # Extract and calculate median and sd for prediction surface also
  elev_pred <- exact_extract(elev, r_cells, progress = FALSE) %>% 
    map_dfr(~ tibble(elevation_median = mean(.$value, na.rm = TRUE),
                     elevation_sd = sd(.$value, na.rm = TRUE))) %>% 
    # join to lookup table to get locality_id
    bind_cols(st_drop_geometry(r_cells), .)
  
  # Join the PLAND and elevation covariates for eBird data 
  pland_elev_checklist <- inner_join(pland, elev_checklists, by = "locality_id")
  pland_elev_pred <- inner_join(pland_coords, elev_pred, by = "id")
  
  #### end ####
  
  #### Elevation Analysis ####
  
  # Plot the elevation data as a map 
  med_elev <- pland_elev_pred %>% 
    # convert to spatial features and project onto prediction surface
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
    st_transform(crs = projection(r)) %>% 
    # rasterize points
    rasterize(r, field = "elevation_median") %>% 
    # project to correct crs
    projectRaster(crs = st_crs(4326)$proj4string, method = "ngb") %>%
    # trim off empty edges of raster
    trim()
  
  # make a map
  pdf(file = file.path(data_folder, "analytics", sprintf("%s_elevation_map.pdf", short_code)))
  par(mar = c(3,3,2,5))
  t <- str_glue("Median Elevation Values")
  plot(ne_land, axes=TRUE, xlim = st_bbox(species_range)[c(1,3)], 
          ylim = st_bbox(species_range)[c(2,4)], main=t)
  plot(species_range, add=TRUE, border='black')
  plot(med_elev, add = TRUE, col = plasma(10), main = t)
  plot(ne_country_lines, add=TRUE)
  dev.off()
  
  #### end ####
  
  
  #### Pressure Maps ####
  
  # Maps created by: https://www.nature.com/articles/s41559-021-01542-9
  # We are only interested in the non-habit covariates, so can disregard logging 
  # and agriculture. Note: 'threats' and 'pressure' used interchangeably
  
  # Load the threat data and keep only non-habitat threats for birds
  threats <- st_read(file.path(perm_files_location, "Grid_mammal_amph_threat_predictions.shp"))
  threats_birds <- select(threats, c("B5_1","B8_1","B9","B11","geometry"))
  
  # Convert the data to a raster so can treat the same way as the elevation data
  threats_birds_rast <- st_rasterize(threats_birds) %>%
                        as("Raster")
  
  # Buffer and crop data to prediction region (100km this time since 50km res)
  threats_birds_rast <- gis_land %>% 
    st_buffer(dist = 100000) %>% 
    st_transform(crs = projection(threats_birds_rast)) %>% 
    crop(threats_birds_rast, .) %>% 
    projectRaster(crs = projection(landcover))
  
  # Replace any NaN values with zero so we have valid data throughout region
  threats_birds_rast[is.na(threats_birds_rast[])] <- 0 
  
  # Create neighbourhood area around each checklist location
  ebird_buff_noyear <- ebird %>% 
    distinct(locality_id, latitude, longitude) %>% 
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
    st_transform(crs = projection(threats_birds_rast)) %>% 
    st_buffer(dist = neighborhood_radius)
  
  # Extract pressure values within neighbourhoods to calculate weighted mean value (based on proportions)
  # (this is different to the basic mean used for elevation data)
  locs <- st_set_geometry(ebird_buff_noyear, NULL) %>% 
    mutate(id = row_number())
  pressure_checklists <- exact_extract(threats_birds_rast, ebird_buff_noyear, progress = TRUE, 'mean') %>% 
    # rename the columns
    rename(
      hunting_mean = mean.B5_1,
      invasives_mean = mean.B8_1,
      pollution_mean = mean.B9,
      climate_mean = mean.B11
    ) %>%
    # join to lookup table to get locality_id
    bind_cols(locs, .) 
  
  # Extract and calculate weighted mean value for prediction surface also
  pressure_pred <- exact_extract(threats_birds_rast, r_cells, progress = TRUE, 'mean') %>% 
    # rename the columns 
    rename(
      hunting_mean = mean.B5_1,
      invasives_mean = mean.B8_1,
      pollution_mean = mean.B9,
      climate_mean = mean.B11
    ) %>% 
    # join to lookup table to get locality_id
    bind_cols(st_drop_geometry(r_cells), .)
  
  # Join the PLAND and elevation covariates, with the pressure data
  pland_elev_pressure_checklist <- inner_join(pland_elev_checklist, pressure_checklists, by = "locality_id")
  pland_elev_pressure_pred <- inner_join(pland_elev_pred, pressure_pred, by = "id")
  
  #### end ####
  
  #### Pressure Maps Analysis ####
  
  # Save a plot of hunting for whole world as illustration
  setEPS()
  postscript(file.path(data_folder, "analytics", "hunting_whole_world.eps"))
  plot(threats_birds["B5_1"], border = 'transparent', 
       pal=plasma, nbreaks=8, key.pos=4, main=NULL, key.length=lcm(8))
  dev.off()
  
  # Also save a plot of hunting threat within the range polygon
  mean_hunting_pressure <- pland_elev_pressure_pred %>% 
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
    st_transform(crs = projection(r)) %>% 
    rasterize(r, field = "hunting_mean") %>% 
    projectRaster(crs = st_crs(4326)$proj4string, method = "ngb") %>%
    trim()
  
  # Plot and save the map
  pdf(file = file.path(data_folder, "analytics", sprintf("%s_pressure_hunting.pdf", short_code)))
  par(mar = c(3,3,2,5))
  t <- str_glue("Mean Hunting Pressure")
  plot(ne_land, axes=TRUE, xlim = st_bbox(species_range)[c(1,3)], 
       ylim = st_bbox(species_range)[c(2,4)], main=t)
  plot(species_range, add=TRUE, border='black')
  plot(mean_hunting_pressure, add=TRUE, col = plasma(10))
  plot(ne_country_lines, add=TRUE)
  dev.off()
  
  #### end ####
  
  
  #### Accessibility Data ####
  
  # Data download from: https://malariaatlas.org/research-project/accessibility-to-cities/
  access <- raster(file.path(perm_files_location, "2015_accessibility_to_cities_v1.0.tif"))
  
  # crop and buffer prediction region by 10 km to provide a little wiggle room
  access <- gis_land %>% 
    st_buffer(dist = 10000) %>% 
    st_transform(crs = projection(access)) %>% 
    crop(access, .) %>% 
    projectRaster(crs = projection(landcover))
  
  # Create neighbourhood area around each checklist location
  ebird_buff_noyear <- ebird %>% 
    distinct(locality_id, latitude, longitude) %>% 
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
    st_transform(crs = projection(access)) %>% 
    st_buffer(dist = neighborhood_radius)
  
  # Extract mean time to city of more than 50,000 people within neighbourhood 
  locs <- st_set_geometry(ebird_buff_noyear, NULL) %>% 
    mutate(id = row_number())
  access_checklists <- exact_extract(access, ebird_buff_noyear, progress = FALSE, 'mean') %>% 
    # join to lookup table to get locality_id
    bind_cols(locs, .) %>%
    # Rename columns 
    rename(access_mean = ...3) 
  
  # Extract mean value for prediction surface also
  access_pred <- exact_extract(access, r_cells, progress = FALSE, 'mean') %>% 
    # join to lookup table to get locality_id
    bind_cols(st_drop_geometry(r_cells), .) %>%
    # Rename columns 
    rename(access_mean = ...2) 
  
  # Join the accessibility data with the rest of the covariate data
  covariate_checklists <- inner_join(pland_elev_pressure_checklist, access_checklists, by = "locality_id")
  covariate_preds <- inner_join(pland_elev_pressure_pred, access_pred, by = "id")
  
  ## Some points are calculated incorrectly due to poor match up of gis and accessibility data
  # For checklists - remove troublesome entries entirely.
  # Get problem checklists from covariate data
  problem_checklists <- covariate_checklists[is.na(covariate_checklists$access_mean), c("locality_id","access_mean")]
  
  # Check that we actually have some problem checklists
  if(nrow(problem_checklists) != 0){
    # For each problem checklist, remove any checklists at that location entirely
    for (row in 1:nrow(problem_checklists)){
      problem_location <- problem_checklists[row,]
      # Remove problem entries in covariate data
      covariate_checklists <- filter(covariate_checklists, locality_id != problem_location$locality_id)
      # Remove problem entries in ebird data
      ebird <- filter(ebird, locality_id != problem_location$locality_id)
    }
  }
  
  # For prediction points - interpolate missing data using nearest valid neighbour's value.
  # Get problem prediction points from covariate data and convert to an sf object
  problem_preds <- covariate_preds[is.na(covariate_preds$access_mean), c("id","latitude","longitude","access_mean")]
  problem_preds <- st_as_sf(problem_preds, coords = c("longitude", "latitude"), crs=4326)
  
  # Get the valid prediction points from covariate data and convert to an sf object 
  valid_preds <- covariate_preds[!is.na(covariate_preds$access_mean), c("latitude","longitude","access_mean")]
  valid_preds <- st_as_sf(valid_preds, coords = c("longitude", "latitude"), crs=4326)
  
  # Check that we actually have some problem prediction points
  if(nrow(problem_preds) != 0){
    # For each problem point, find its nearest valid neighbour and use its accessibility value
    for (row in 1:nrow(problem_preds)){
      problem_pred <- problem_preds[row,]
      # Find distances to each valid row
      distances <- st_distance(problem_pred[c("geometry","access_mean")], valid_preds)
      # Get the index of the closest point
      index_min <- which.min(distances)
      # Get accessibility at this closest point
      new_access <- st_drop_geometry(valid_preds[index_min,"access_mean"])
      # Save this value as the new value of the corresponding problem checklist
      covariate_preds[covariate_preds$id == problem_pred$id, "access_mean"] <- new_access
    }
  }
  
  # Save the updated covariate data
  write_csv(covariate_checklists, file.path(data_folder, "output data",
              sprintf("%s_covariate_checklists.csv", short_code)))
  write_csv(covariate_preds, file.path(data_folder, "output data",
              sprintf("%s_covariate_predictions.csv", short_code)))
  
  # Finally, save the updated eBird data also
  csv_name <- file.path(data_folder, "output data", sprintf("%s_ebd_processed_output.csv", short_code))
  write_csv(ebird, csv_name, na = "")
  
  #### end ####
  
  #### Accessibility Analysis ####
  
  # Plot the accessibility data as a map 
  access <- covariate_preds %>% 
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
    st_transform(crs = projection(r)) %>% 
    rasterize(r, field = "access_mean") %>% 
    projectRaster(crs = st_crs(4326)$proj4string, method = "ngb") %>%
    trim()

  # make a map
  pdf(file = file.path(data_folder, "analytics", sprintf("%s_accessibility_map.pdf", short_code)))
  par(mar = c(3,3,2,5))
  t <- str_glue("Mean Accessibility Values (mins)")
  plot(ne_land, axes=TRUE, xlim = c(95,107), 
       ylim = c(5,20), main=t)
  plot(species_range, add=TRUE, border='black')
  plot(access, add = TRUE, col = viridis(10), main = t)
  plot(ne_country_lines, add=TRUE)
  dev.off()
  
  #### end ####
  
  
  #### eBird Analysis ####
  
  # Plots only generated at this point due to updated eBird data from accessibility section
  
  # Start by producing a spatial map of all checklists recorded
  # Prepare ebird data for mapping
  ebird_sf <- ebird %>% 
    # convert to spatial points
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
    select(species_observed)
  
  # Set up a plot area and plot data/GIS
  pdf(file = file.path(data_folder, "analytics", sprintf("%s_eBird_checklists_map.pdf", short_code)))
  plot(ne_land, axes=TRUE, xlim = st_bbox(species_range)[c(1,3)], 
       ylim = st_bbox(species_range)[c(2,4)], main=paste(current_species$common_name, "eBird Observations"))
  plot(species_range, add=TRUE, col=alpha("orange",0.7), border='transparent')
  plot(ne_country_lines, add=TRUE)
  
  # eBird observations - colour black as not-observed and green as observed
  plot(st_geometry(ebird_sf),
       pch = 19, cex = 0.1, col = alpha("#555555", 0.25),
       add = TRUE)
  plot(filter(ebird_sf, species_observed) %>% st_geometry(),
       pch = 19, cex = 0.4, col = alpha("#4daf4a", 1),
       add = TRUE)
  
  # Add a legend and title
  legend("bottomright", bty = "n",
         col = c("#555555", "#4daf4a"),
         legend = c("eBird checklists", paste(current_species$common_name, "sightings")),
         pch = 19)
  dev.off()
  
  
  # Now lets explore how observations of the species are effected by each of the 
  # effort variables 
  
  # Time of day
  # summarize data by hourly bins
  breaks <- 0:24
  labels <- breaks[-length(breaks)] + diff(breaks) / 2
  ebird_tod <- ebird %>% 
    mutate(tod_bins = cut(time_observations_started, 
                          breaks = breaks, 
                          labels = labels,
                          include.lowest = TRUE),
           tod_bins = as.numeric(as.character(tod_bins))) %>% 
    group_by(tod_bins) %>% 
    summarise(n_checklists = n(),
              n_detected = sum(species_observed),
              det_freq = mean(species_observed))
  
  # histogram
  pdf(file = file.path(data_folder, "analytics", sprintf("%s_effort_time-of-day.pdf", short_code)))
  g_tod_hist <- ggplot(ebird_tod) +
    aes(x = tod_bins, y = n_checklists) +
    geom_segment(aes(xend = tod_bins, y = 0, yend = n_checklists),
                 color = "grey50") +
    geom_point() +
    scale_x_continuous(breaks = seq(0, 24, by = 3), limits = c(0, 24)) +
    scale_y_continuous(labels = scales::comma) +
    labs(x = "Hours since midnight",
         y = "# checklists",
         title = "Distribution of observation start times")
  
  # frequency of detection
  g_tod_freq <- ggplot(ebird_tod %>% filter(n_checklists > 100)) +
    aes(x = tod_bins, y = det_freq) +
    geom_line() +
    geom_point() +
    scale_x_continuous(breaks = seq(0, 24, by = 3), limits = c(0, 24)) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "Hours since midnight",
         y = "% checklists with detections",
         title = "Detection frequency")
  
  # combine
  grid.arrange(g_tod_hist, g_tod_freq)
  dev.off()
  
  
  # Checklist Duration
  # summarize data by 30 minute bins
  breaks <- seq(0, 5, by = 0.5)
  labels <- breaks[-length(breaks)] + diff(breaks) / 2
  ebird_dur <- ebird %>% 
    mutate(dur_bins = cut(duration_minutes / 60, 
                          breaks = breaks, 
                          labels = labels,
                          include.lowest = TRUE),
           dur_bins = as.numeric(as.character(dur_bins))) %>% 
    group_by(dur_bins) %>% 
    summarise(n_checklists = n(),
              n_detected = sum(species_observed),
              det_freq = mean(species_observed))
  
  # histogram
  pdf(file = file.path(data_folder, "analytics", sprintf("%s_effort_duration.pdf", short_code)))
  g_dur_hist <- ggplot(ebird_dur) +
    aes(x = dur_bins, y = n_checklists) +
    geom_segment(aes(xend = dur_bins, y = 0, yend = n_checklists),
                 color = "grey50") +
    geom_point() +
    scale_x_continuous(breaks = 0:5) +
    scale_y_continuous(labels = scales::comma) +
    labs(x = "Checklist duration (hours)",
         y = "# checklists",
         title = "Distribution of checklist durations")
  
  # frequency of detection
  g_dur_freq <- ggplot(ebird_dur %>% filter(n_checklists > 100)) +
    aes(x = dur_bins, y = det_freq) +
    geom_line() +
    geom_point() +
    scale_x_continuous(breaks = 0:5) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "Checklist duration (hours)",
         y = "% checklists with detections",
         title = "Detection frequency")
  
  # combine
  grid.arrange(g_dur_hist, g_dur_freq)
  dev.off()
  
  
  # Distance Traveled 
  # summarize data by 500m bins
  breaks <- seq(0, 5, by = 0.5)
  labels <- breaks[-length(breaks)] + diff(breaks) / 2
  ebird_dist <- ebird %>% 
    mutate(dist_bins = cut(effort_distance_km, 
                           breaks = breaks, 
                           labels = labels,
                           include.lowest = TRUE),
           dist_bins = as.numeric(as.character(dist_bins))) %>% 
    group_by(dist_bins) %>% 
    summarise(n_checklists = n(),
              n_detected = sum(species_observed),
              det_freq = mean(species_observed))
  
  # histogram
  pdf(file = file.path(data_folder, "analytics", sprintf("%s_effort_distance-travelled.pdf", short_code)))
  g_dist_hist <- ggplot(ebird_dist) +
    aes(x = dist_bins, y = n_checklists) +
    geom_segment(aes(xend = dist_bins, y = 0, yend = n_checklists),
                 color = "grey50") +
    geom_point() +
    scale_x_continuous(breaks = 0:5) +
    scale_y_continuous(labels = scales::comma) +
    labs(x = "Distance travelled (km)",
         y = "# checklists",
         title = "Distribution of distance travelled")
  
  # frequency of detection
  g_dist_freq <- ggplot(ebird_dist %>% filter(n_checklists > 100)) +
    aes(x = dist_bins, y = det_freq) +
    geom_line() +
    geom_point() +
    scale_x_continuous(breaks = 0:5) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "Distance travelled (km)",
         y = "% checklists with detections",
         title = "Detection frequency")
  
  # combine
  grid.arrange(g_dist_hist, g_dist_freq)
  dev.off()
  
  # Since most walks are short, landcover neighborhoods at location will describe 
  # habitat at sighting location well!
  
  
  # Number of Observers 
  # summarize data
  breaks <- 0:10
  labels <- 1:10
  ebird_obs <- ebird %>% 
    mutate(obs_bins = cut(number_observers, 
                          breaks = breaks, 
                          label = labels,
                          include.lowest = TRUE),
           obs_bins = as.numeric(as.character(obs_bins))) %>% 
    group_by(obs_bins) %>% 
    summarise(n_checklists = n(),
              n_detected = sum(species_observed),
              det_freq = mean(species_observed))
  
  # histogram
  pdf(file = file.path(data_folder, "analytics", sprintf("%s_effort_number-of-observers.pdf", short_code)))
  g_obs_hist <- ggplot(ebird_obs) +
    aes(x = obs_bins, y = n_checklists) +
    geom_segment(aes(xend = obs_bins, y = 0, yend = n_checklists),
                 color = "grey50") +
    geom_point() +
    scale_x_continuous(breaks = 1:10) +
    scale_y_continuous(labels = scales::comma) +
    labs(x = "# observers",
         y = "# checklists",
         title = "Distribution of the number of observers")
  
  # frequency of detection
  g_obs_freq <- ggplot(ebird_obs %>% filter(n_checklists > 100)) +
    aes(x = obs_bins, y = det_freq) +
    geom_line() +
    geom_point() +
    scale_x_continuous(breaks = 1:10) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "# observers",
         y = "% checklists with detections",
         title = "Detection frequency")
  
  # combine
  grid.arrange(g_obs_hist, g_obs_freq)
  dev.off()
  
  #### end ####
  
  
  # Clean up the directory 
  unlink(file.path(data_folder, "ebd_filtered_output.txt"))
  unlink(file.path(data_folder, "ebd_samp_filtered_output.txt"))
}



