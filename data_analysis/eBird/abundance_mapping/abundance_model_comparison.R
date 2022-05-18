# Hamish Campbell

#> Compares the accuracy of the "standard" model without any non-habitat covariates
#> with a model that includes such data i.e. pressure maps etc. 
#> 
#> Heavily based on code written and explained here: 
#> #> https://cornelllabofornithology.github.io/ebird-best-practices/abundance.html

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

# Select which packages want to take duplicated function names from
select <- dplyr::select
map <- purrr::map
projection <- raster::projection

# Set random number seed to ensure repeatable results
set.seed(1)

# Load the variable values from the config file
source("config.R")

# Setup output directories if don't exist
if (!dir.exists(file.path(data_folder, "output data"))){
  dir.create(file.path(data_folder, "output data"))
}
if (!dir.exists(file.path(data_folder, "analytics"))){
  dir.create(file.path(data_folder, "analytics"))
}


#### Data Prep/Loading #### 

# Load eBird data 
ebird_data <- read_csv(file.path(data_folder, "input data", "ebd_processed_output.csv")) %>% 
  # Make the 'protocol type' column categorical with 2 options 
  mutate(protocol_type = factor(protocol_type, 
                                levels = c("Stationary" , "Traveling"))) %>%
  # remove observations with no count
  filter(!is.na(observation_count))

# Load covariate data
covariate_data <- read_csv(file.path(data_folder, "input data", "landcover_elevation_and_pressure_checklists.csv")) %>% 
  mutate(year = as.integer(year))

# Combine eBird and covariate data 
checklist_data <- inner_join(ebird_data, covariate_data, by = c("locality_id", "year"))

# Load prediction surface covariate data
pred_covariate_data <- read_csv(file.path(data_folder, "input data", "landcover_elevation_and_pressure_prediction.csv"))

# Load prediction area raster
pred_raster <- raster(file.path(data_folder, "input data", "prediction-surface.tif"))

# Load GIS data for making maps (range, land and country boundaries)
range <- read_sf(file.path(data_folder, "input data", "gis-data.gpkg"), "range") %>% 
  st_geometry()
ne_land <- read_sf(file.path(data_folder, "input data", "gis-data.gpkg"), "ne_land") %>% 
  st_geometry()
ne_country_lines <- read_sf(file.path(data_folder, "input data", "gis-data.gpkg"), "ne_country_lines") %>% 
  st_geometry()

#### end ####


#### Spatiotemporal Subsampling ####

# This step helps reduce temporal and spatial biases in eBird checklist data
# and reduces class imbalance since relatively more non-sightings are removed

# Generate hexagonal mosaic with ~ 5 km between cells 
dggs <- dgconstruct(spacing = 5)

# Assign hexagonal cell id and week number for each checklist in dataset
checklist_cell <- checklist_data %>% 
  mutate(cell = dgGEO_to_SEQNUM(dggs, longitude, latitude)$seqnum,
         week = week(observation_date))

# Add one checklist sample per grid cell, per week to new dataset
checklist_data_ss <- checklist_cell %>% 
  group_by(species_observed, year, week, cell) %>% 
  sample_n(size = 1) %>% 
  ungroup() %>% 
  select(-cell, -week)

#### end ####


#### Baseline Model Fit ####

# Keep only the data to be used in the model
baseline_model_data <- checklist_data_ss %>% 
  select(observation_count,
         # Effort covariates
         day_of_year, time_observations_started, duration_minutes,
         effort_distance_km, number_observers, protocol_type,
         # Habitat covariates
         habitat_covariates)

# Split the model data into a testing and training set
baseline_model_data <- baseline_model_data %>% 
  split(if_else(runif(nrow(.)) <= train_prop, "train", "test"))

# Print number of datapoints in each set
map_int(baseline_model_data, nrow)

# Get names of continuous covariates 
continuous_covs <- baseline_model_data$train %>% 
  select(-observation_count, -protocol_type, -time_observations_started) %>% 
  names()

# Create GAM formula for count predictions
gam_formula_rhs <- str_glue("s({var}, k = {k})", 
                            var = continuous_covs, k = k) %>% 
  str_flatten(collapse = " + ") %>% 
  str_glue(" ~ ", .,
           " + protocol_type + ",
           "s(time_observations_started, bs = \"cc\", k = {k})", 
           k = k_time) %>% 
  as.formula()

# Add response variable (count) to GAM formula 
gam_formula <- update.formula(observation_count ~ ., gam_formula_rhs)

# Required to ensure hours 0 and 24 are mapped to each other 
time_knots <- list(time_observations_started = seq(0, 24, length.out = k_time))

# Fit the baseline GAM model using a 'negative binomial' distribution 
NB_baseline_model <- gam(gam_formula, data = baseline_model_data$train, family = "nb", knots = time_knots)

#### end ####


#### Advanced Model Fit ####

# Keep only the data to be used in the model
adv_model_data <- checklist_data_ss %>% 
  select(observation_count,
         # Effort covariates
         day_of_year, time_observations_started, duration_minutes,
         effort_distance_km, number_observers, protocol_type,
         # Habitat covariates
         habitat_covariates, 
         # Non-habitat covariates 
         non_habitat_covariates) 

# Split the model data into a testing and training set
adv_model_data <- adv_model_data %>% 
  split(if_else(runif(nrow(.)) <= train_prop, "train", "test"))

# Print number of datapoints in each set
map_int(adv_model_data, nrow)

# Get names of continuous covariates 
continuous_covs <- adv_model_data$train %>% 
  select(-observation_count, -protocol_type, -time_observations_started) %>% 
  names()

# Create GAM formula for count predictions
gam_formula_rhs <- str_glue("s({var}, k = {k})", 
                            var = continuous_covs, k = k) %>% 
  str_flatten(collapse = " + ") %>% 
  str_glue(" ~ ", .,
           " + protocol_type + ",
           "s(time_observations_started, bs = \"cc\", k = {k})", 
           k = k_time) %>% 
  as.formula()

# Add response variable (count) to GAM formula 
gam_formula <- update.formula(observation_count ~ ., gam_formula_rhs)

# Required to ensure hours 0 and 24 are mapped to each other 
time_knots <- list(time_observations_started = seq(0, 24, length.out = k_time))

# Fit the advanced GAM model using a 'negative binomial' distribution 
NB_adv_model <- gam(gam_formula, data = adv_model_data$train, family = "nb", knots = time_knots)

#### end ####


#### Model Comparison ####

# Get the count data from the test sets
baseline_count_test_data <- select(baseline_model_data$test, obs = observation_count)
adv_count_test_data <- select(adv_model_data$test, obs = observation_count)

# Make predictions for test set using both GAM models
baseline_model_predictions <- predict(NB_baseline_model, baseline_model_data$test, type = "response") %>% 
  tibble(model_type = "Baseline", pred = .) %>% 
  bind_cols(baseline_count_test_data)

adv_model_predictions <- predict(NB_adv_model, adv_model_data$test, type = "response") %>% 
  tibble(model_type = "Advanced", pred = .) %>% 
  bind_cols(adv_count_test_data)

# combine predictions from each of the models
test_predictions <- bind_rows(baseline_model_predictions, adv_model_predictions) %>% 
  mutate(model_type = as_factor(model_type))

# Compare models using Spearman’s rank correlation (higher = better)
spearman_ranks <- test_predictions %>% 
                  group_by(model_type) %>% 
                  summarise(rank_cor = cor.test(obs, pred, 
                                                method = "spearman", 
                                                exact = FALSE)$estimate) %>% 
                  ungroup()

# Compare models using mean absolute deviation (smaller = better)
mad_ranks <- test_predictions %>% 
             group_by(model_type) %>% 
             summarise(mad = mean(abs(obs - pred), na.rm = TRUE)) %>% 
             ungroup()

# Save the comparison metrics
write.csv(spearman_ranks, file.path(data_folder, "analytics", "spearman_ranks.csv"), row.names = FALSE)
write.csv(mad_ranks, file.path(data_folder, "analytics", "mad_ranks.csv"), row.names = FALSE)

#### end #### 


#### Baseline Model Analytics ####

## Count distributions (applies to both models)
# Setup plot(s) of distribution of count numbers 
pdf(file = file.path(data_folder, "analytics", "count_distributions.pdf"))
p <- par(mfrow = c(1, 2))

# Plot distribution with zero counts included
hist(checklist_data_ss$observation_count, main = "Histogram of counts", 
     xlab = "Observed count")

# Plot distribution without zero counts included
pos_counts <- keep(checklist_data_ss$observation_count, ~ . > 0)
hist(pos_counts, main = "Histogram of counts > 0", 
     xlab = "Observed non-zero count")
dev.off()


## Optimal time of day for observations  
# Get latest year of landcover data available 
max_lc_year <- pred_covariate_data$year[1]

# Find optimal time of day for observing the species 
seq_tod <- seq(0, 24, length.out = 300)
tod_df <- baseline_model_data$train %>% 
  # Find average of continuous covariates 
  select(habitat_covariates) %>% 
  summarize_all(mean, na.rm = TRUE) %>% 
  ungroup() %>% 
  # Use a 'standard' checklist
  mutate(day_of_year = yday(ymd(str_glue("{max_lc_year}-06-15"))),
         duration_minutes = 60,
         effort_distance_km = 1,
         number_observers = 1,
         protocol_type = "Traveling") %>% 
  cbind(time_observations_started = seq_tod)

# Predict counts at start times throughout the day
pred_tod <- predict(NB_baseline_model, newdata = tod_df, 
                    type = "link", 
                    se.fit = TRUE) %>% 
  as_tibble() %>% 
  # Calculate backtransformed confidence limits
  transmute(time_observations_started = seq_tod,
            pred = NB_baseline_model$family$linkinv(fit),
            pred_lcl = NB_baseline_model$family$linkinv(fit - 1.96 * se.fit),
            pred_ucl = NB_baseline_model$family$linkinv(fit + 1.96 * se.fit))

# Optimal time of day is when lower conf interval is maximised
t_peak <- pred_tod$time_observations_started[which.max(pred_tod$pred_lcl)]

# Plot the findings to illustrate the result 
pdf(file = file.path(data_folder, "analytics", "baseline_optimal_time_of_day.pdf"))
ggplot(pred_tod) +
  aes(x = time_observations_started, y = pred,
      ymin = pred_lcl, ymax = pred_ucl) +
  geom_ribbon(fill = "grey80", alpha = 0.5) +
  geom_line() +
  geom_vline(xintercept = t_peak, color = "blue", linetype = "dashed") +
  labs(x = "Hours since midnight",
       y = "Predicted relative abundance",
       title = paste("Effect of observation start time on", species_name, "reporting"),
       subtitle = "Peak detectability shown as dashed blue line")
dev.off()
  
  
## Abundance maps 
# Add effort covariates to those already have
pred_covariate_effort_data <- pred_covariate_data %>% 
  mutate(day_of_year = yday(ymd(str_glue("{max_lc_year}-", time_of_year))),
         time_observations_started = t_peak,
         duration_minutes = 60,
         effort_distance_km = 1,
         number_observers = 1,
         protocol_type = "Traveling")
  
  
# Now make predictions over the complete surface 
pred <- predict(NB_baseline_model, newdata = pred_covariate_effort_data, 
                type = "link", 
                se.fit = TRUE) %>% 
  as_tibble() %>% 
  # Get abundance and standard error in abundance 
  transmute(abd = NB_baseline_model$family$linkinv(fit),
            abd_se = NB_baseline_model$family$linkinv(se.fit)) %>%
  # Add to lat/lon from prediction surface
  bind_cols(pred_covariate_effort_data, .) %>% 
  select(latitude, longitude, abd, abd_se)
  
# Convert predictions into a raster file
pred_map <- pred %>% 
  # convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  select(abd, abd_se) %>% 
  st_transform(crs = projection(pred_raster)) %>% 
  rasterize(pred_raster)
pred_map <- pred_map[[-1]]
  
# Save prediction maps to output folder
writeRaster(pred_map[["abd"]], 
            filename = file.path(data_folder, "output data", "baseline_abundance_map.tif"),
            overwrite = TRUE)
writeRaster(pred_map[["abd_se"]], 
            filename = file.path(data_folder, "output data", "baseline_abundance_map_uncertainty.tif"), 
            overwrite = TRUE)
  
  
# Change the CRS of the raster to match that of the GIS data for plotting
pred_map_proj <- projectRaster(pred_map, crs = 4326, method = "ngb")

# Plot the data for each layer (abd and abd_se)
pdf(file = file.path(data_folder, "analytics", "baseline_abundance_maps.pdf"))
par(mfrow=c(1,2))
for (nm in names(pred_map)) {
  r_plot <- pred_map_proj[[nm]]
    
  # Set the plot area and add GIS components
  par(mar = c(4.5, 0.25, 0.25, 0.25))
  plot(ne_land, xlim = st_bbox(range)[c(1,3)], 
       ylim = st_bbox(range)[c(2,4)])
  
  # Create a colour palette
  plasma_rev <- rev(plasma(25, end = 0.9))
  gray_int <- colorRampPalette(c("#dddddd", plasma_rev[1]))
  pal <- c(gray_int(4)[2], plasma_rev)
  
  # Setup plot characetristics depending on whether abd or abd_se
  if (nm == "abd") {
    title <- paste(species_name, "Relative Abundance")
    # Set very low values to zero
    r_plot[r_plot <= zero_threshold] <- NA
    # Log transform
    r_plot <- log10(r_plot)
    # Breaks and legend
    mx <- ceiling(100 * cellStats(r_plot, max)) / 100
    mn <- floor(100 * cellStats(r_plot, min)) / 100
    brks <- seq(mn, mx, length.out = length(pal) + 1)
    lbl_brks <- sort(c(-2:2, mn, mx))
    lbls <- round(10^lbl_brks, 2)
  } else {
    title <- paste(species_name, "Abundance Uncertainty (SE)")
    # Hard code a few uncertainty values are massive 
    mx <- 4 
    mn <- 0
    brks <- seq(mn, mx, length.out = length(pal) + 1)
    lbl_brks <- seq(mn, mx, length.out = 5)
    lbls <- round(lbl_brks, 2)
  }
  
  # Plot the data
  plot(r_plot, col = pal, breaks = brks, maxpixels = ncell(r_plot),
          legend = FALSE, add = TRUE, border='black')
  
  # Add country lines for GIS
  plot(ne_country_lines, add=TRUE)
  box()
  
  # Legend
  par(new = TRUE, mar = c(0, 0, 0, 0))
  image.plot(zlim = range(brks), legend.only = TRUE, col = pal,
             smallplot = c(0.25, 0.75, 0.06, 0.09),
             horizontal = TRUE,
             axis.args = list(at = lbl_brks, 
                              labels = lbls,
                              fg = "black", col.axis = "black",
                              cex.axis = 0.75, lwd.ticks = 0.5,
                              padj = -1.5),
             legend.args = list(text = title,
                                side = 3, col = "black",
                                cex = 0.8, line = 0))
}
dev.off()

#### end ####


#### Advanced Model Analytics ####

## Optimal time of day for observations  
# Get latest year of landcover data available 
max_lc_year <- pred_covariate_data$year[1]

# Find optimal time of day for observing the species 
seq_tod <- seq(0, 24, length.out = 300)
tod_df <- adv_model_data$train %>% 
  # Find average of continuous covariates 
  select(habitat_covariates, non_habitat_covariates) %>% 
  summarize_all(mean, na.rm = TRUE) %>% 
  ungroup() %>% 
  # Use a 'standard' checklist
  mutate(day_of_year = yday(ymd(str_glue("{max_lc_year}-06-15"))),
         duration_minutes = 60,
         effort_distance_km = 1,
         number_observers = 1,
         protocol_type = "Traveling") %>% 
  cbind(time_observations_started = seq_tod)

# Predict counts at start times throughout the day
pred_tod <- predict(NB_adv_model, newdata = tod_df, 
                    type = "link", 
                    se.fit = TRUE) %>% 
  as_tibble() %>% 
  # Calculate backtransformed confidence limits
  transmute(time_observations_started = seq_tod,
            pred = NB_adv_model$family$linkinv(fit),
            pred_lcl = NB_adv_model$family$linkinv(fit - 1.96 * se.fit),
            pred_ucl = NB_adv_model$family$linkinv(fit + 1.96 * se.fit))

# Optimal time of day is when lower conf interval is maximised
t_peak <- pred_tod$time_observations_started[which.max(pred_tod$pred_lcl)]

# Plot the findings to illustrate the result 
pdf(file = file.path(data_folder, "analytics", "advanced_optimal_time_of_day.pdf"))
ggplot(pred_tod) +
  aes(x = time_observations_started, y = pred,
      ymin = pred_lcl, ymax = pred_ucl) +
  geom_ribbon(fill = "grey80", alpha = 0.5) +
  geom_line() +
  geom_vline(xintercept = t_peak, color = "blue", linetype = "dashed") +
  labs(x = "Hours since midnight",
       y = "Predicted relative abundance",
       title = paste("Effect of observation start time on", species_name, "reporting"),
       subtitle = "Peak detectability shown as dashed blue line")
dev.off()


## Abundance maps 
# Add effort covariates to those already have
pred_covariate_effort_data <- pred_covariate_data %>% 
  mutate(day_of_year = yday(ymd(str_glue("{max_lc_year}-", time_of_year))),
         time_observations_started = t_peak,
         duration_minutes = 60,
         effort_distance_km = 1,
         number_observers = 1,
         protocol_type = "Traveling")


# Now make predictions over the complete surface 
pred <- predict(NB_adv_model, newdata = pred_covariate_effort_data, 
                type = "link", 
                se.fit = TRUE) %>% 
  as_tibble() %>% 
  # Get abundance and standard error in abundance 
  transmute(abd = NB_adv_model$family$linkinv(fit),
            abd_se = NB_adv_model$family$linkinv(se.fit)) %>%
  # Add to lat/lon from prediction surface
  bind_cols(pred_covariate_effort_data, .) %>% 
  select(latitude, longitude, abd, abd_se)

# Convert predictions into a raster file
pred_map <- pred %>% 
  # convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  select(abd, abd_se) %>% 
  st_transform(crs = projection(pred_raster)) %>% 
  rasterize(pred_raster)
pred_map <- pred_map[[-1]]

# Save prediction maps to output folder
writeRaster(pred_map[["abd"]], 
            filename = file.path(data_folder, "output data", "advanced_abundance_map.tif"),
            overwrite = TRUE)
writeRaster(pred_map[["abd_se"]], 
            filename = file.path(data_folder, "output data", "advanced_abundance_map_uncertainty.tif"), 
            overwrite = TRUE)


# Change the CRS of the raster to match that of the GIS data for plotting
pred_map_proj <- projectRaster(pred_map, crs = 4326, method = "ngb")

# Plot the data for each layer (abd and abd_se)
pdf(file = file.path(data_folder, "analytics", "advanced_abundance_maps.pdf"))
par(mfrow=c(1,2))
for (nm in names(pred_map)) {
  r_plot <- pred_map_proj[[nm]]
  
  # Set the plot area and add GIS components
  par(mar = c(4.5, 0.25, 0.25, 0.25))
  plot(ne_land, xlim = st_bbox(range)[c(1,3)], 
       ylim = st_bbox(range)[c(2,4)])
  
  # Create a colour palette
  plasma_rev <- rev(plasma(25, end = 0.9))
  gray_int <- colorRampPalette(c("#dddddd", plasma_rev[1]))
  pal <- c(gray_int(4)[2], plasma_rev)
  
  # Setup plot characetristics depending on whether abd or abd_se
  if (nm == "abd") {
    title <- paste(species_name, "Relative Abundance")
    # Set very low values to zero
    r_plot[r_plot <= zero_threshold] <- NA
    # Log transform
    r_plot <- log10(r_plot)
    # Breaks and legend
    mx <- ceiling(100 * cellStats(r_plot, max)) / 100
    mn <- floor(100 * cellStats(r_plot, min)) / 100
    brks <- seq(mn, mx, length.out = length(pal) + 1)
    lbl_brks <- sort(c(-2:2, mn, mx))
    lbls <- round(10^lbl_brks, 2)
  } else {
    title <- paste(species_name, "Abundance Uncertainty (SE)")
    mx <- 4
    mn <- 0
    brks <- seq(mn, mx, length.out = length(pal) + 1)
    lbl_brks <- seq(mn, mx, length.out = 5)
    lbls <- round(lbl_brks, 2)
  }
  
  # Plot the data
  plot(r_plot, col = pal, breaks = brks, maxpixels = ncell(r_plot),
       legend = FALSE, add = TRUE, border='black')
  
  # Add country lines for GIS
  plot(ne_country_lines, add=TRUE)
  box()
  
  # Legend
  par(new = TRUE, mar = c(0, 0, 0, 0))
  image.plot(zlim = range(brks), legend.only = TRUE, col = pal,
             smallplot = c(0.25, 0.75, 0.06, 0.09),
             horizontal = TRUE,
             axis.args = list(at = lbl_brks, 
                              labels = lbls,
                              fg = "black", col.axis = "black",
                              cex.axis = 0.75, lwd.ticks = 0.5,
                              padj = -1.5),
             legend.args = list(text = title,
                                side = 3, col = "black",
                                cex = 0.8, line = 0))
}
dev.off()

#### end ####

