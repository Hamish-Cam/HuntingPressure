# Hamish Campbell

#> Creates abundance maps such as those given by the eBird guide at:
#> https://cornelllabofornithology.github.io/ebird-best-practices/abundance.html

# Core Processes:
#> 1. Data Prep/Loading
#> 2. Spatiotemporal Subsampling
#> 3. Exploratory Data Analysis
#> 4. Abundance Model (GAM)
#> 5. Model Critique
#> 6. Abundance Prediction Mapping

# Need to ensure that we have the correct packages installed 
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

# resolve namespace conflicts (multiple packages use same names for functions)
select <- dplyr::select
map <- purrr::map
projection <- raster::projection

# set random number seed to insure fully repeatable results
set.seed(1)

# Load the variable values from the config file
source("config_abun.R")

# setup output directory for saved results (abundance tifs)
if (!dir.exists(file.path(data_folder, "output data"))){
  dir.create(file.path(data_folder, "output data"))
}

if (!dir.exists(file.path(data_folder, "analytics"))){
  dir.create(file.path(data_folder, "analytics"))
}


#### Data Prep/Loading #### 

# ebird data
ebird <- read_csv(file.path(data_folder, "input data", "ebd_processed_output.csv")) %>% 
  # Make the 'protocol type' column categorical with 2 options 
  mutate(protocol_type = factor(protocol_type, 
                                levels = c("Stationary" , "Traveling"))) %>%
  # remove observations with no count
  filter(!is.na(observation_count))

# modis habitat covariates (landcover proportion and elevation information 
# summarized by neighbourhoods around checklist locations)
habitat <- read_csv(file.path(data_folder, "input data", "landcover_and_elevation_checklists.csv")) %>% 
  mutate(year = as.integer(year))

# combine ebird and habitat data (add env data for each observation)
ebird_habitat <- inner_join(ebird, habitat, by = c("locality_id", "year"))

# prediction surface - this is a regular grid of landcover type/elevation (our
# habitat covariates) that we can use to make predictions for the whole prediction region
# NOTE: this is NOT the same thing as the habitat covariate data seen
pred_surface <- read_csv(file.path(data_folder, "input data", "landcover_and_elevation_prediction.csv"))

# latest year of landcover data available 
max_lc_year <- pred_surface$year[1]

# This is a raster file defining the prediction area (uses the same crs/grid as 
# the rest of the prediction surface data)
r <- raster(file.path(data_folder, "input data", "prediction-surface.tif"))

# load gis data for making maps (land boundary and possibly country lines)
ne_land <- read_sf(file.path(data_folder, "input data", "gis-data.gpkg"), "ne_land") %>% 
  st_geometry()
#ne_country_lines <- read_sf(file.path(data_folder, "gis-data.gpkg"), "ne_country_lines") %>% 
  #st_geometry()

#### end ####


#### Spatiotemporal Subsampling ####
# This step is necessary to help remove the temporal and spatial biases in the 
# eBird data (more sightings during nesting, and western society e.g.). 
# Also helps remove class imbalance (many more non-sightings reported)

# generate hexagonal grid with ~ 5 km between cells (mosaic)
dggs <- dgconstruct(spacing = 5)

# get hexagonal cell id (in mosaic) and week number for each checklist
# Note: uses longitude and latitude columns from ebird_habitat and produces a 
# new dggs object, then extract the seqnum column since we want this. 'week' is 
# a function that converts dates to week number in year
checklist_cell <- ebird_habitat %>% 
  mutate(cell = dgGEO_to_SEQNUM(dggs, longitude, latitude)$seqnum,
         week = week(observation_date))

# subsample one checklist per grid cell, per week
# sample detection/non-detection independently (don't include in grouping)
# Note: groups table entries using these columns, then takes one sample from each group,
# ungroups again and then removes the cell and week columns again
ebird_ss <- checklist_cell %>% 
  group_by(species_observed, year, week, cell) %>% 
  sample_n(size = 1) %>% 
  ungroup() %>% 
  select(-cell, -week)

# Now we have a subsampled dataset that accounts for spatial, temporal biases.
# My belief is that since we didn't sepearte the detection/non-detection that the 
# class imbalance is still present since we have chosen samples randomly

#### end ####


#### Exploratory Data Analysis ####
# Plotting distributions of data indicates which distributions may be appropriate 
# for modeling the counts of this species. We will plot histograms of numbers of 
# counts registered (one including 0's and one not)

# Plot hists side by side 
png(file = file.path(data_folder, "analytics", "count_distributions.png"))
p <- par(mfrow = c(1, 2))

# counts with zeros
hist(ebird_ss$observation_count, main = "Histogram of counts", 
     xlab = "Observed count")

# counts without zeros
pos_counts <- keep(ebird_ss$observation_count, ~ . > 0)
hist(pos_counts, main = "Histogram of counts > 0", 
     xlab = "Observed non-zero count")
dev.off()

#### end ####


#### Abundance Model (GAM) ####

# Note: This code is built such that only the selected covariates from the previous
# section are used and it is not hard coded so changes to these for another species
# will be implemented automatically. 
#
# GAMs are a method of approximating a non-linear relationship between covariates (x)
# and a response variables (y). It does so by using each covariate to make a set of
# basis functions (various forms possible e.g. polynomial: x, x^2, ...). Then coefficients
# are fit to these basis functions to get a best fit. The degrees of freedom (k) is 
# the number of basis functions used for each covariate (higher k = more wiggly possible
# function).
#
# Note: we treat the 'start time' covariate differently to the others since we want
# it to be cyclic (0 hours should be mapped to midnight = 24). 

# Keep only columns used for analysis 
ebird_split <- ebird_ss %>% 
  select(observation_count,
         # Keep columns used to estimate 'effort' of detection
         day_of_year, time_observations_started, duration_minutes,
         effort_distance_km, number_observers, protocol_type,
         # Keep relevant habitat covariates (chosen in config)
         habitat_covariates) 

# select continuous predictors (remove response and time start)
continuous_covs <- ebird_split %>% 
  select(-observation_count, -protocol_type, -time_observations_started) %>% 
  names()

# create model formula for predictors - just a set of string manipulations to get
# all covariates and corresponding degrees of freedom into correct form for GAM 
# formula. 's' means a smooth term (sum of k basis functions)
gam_formula_rhs <- str_glue("s({var}, k = {k})", 
                            var = continuous_covs, k = k) %>% 
  str_flatten(collapse = " + ") %>% 
  # Add ~ to front and covariates we had to treat seperately to end. For 'time started'
  # use a cubic cyclic spline, for protocol type treat as categorical
  str_glue(" ~ ", .,
           " + protocol_type + ",
           "s(time_observations_started, bs = \"cc\", k = {k})", 
           k = k_time) %>% 
  as.formula()

# Add response (observation count) to GAM formula
gam_formula <- update.formula(observation_count ~ ., gam_formula_rhs)

# Example Result: 
#> observation_count ~ s(day_of_year, k = 5) + s(duration_minutes, 
#>     k = 5) + s(effort_distance_km, k = 5) + s(number_observers, 
#>     k = 5) + s(pland_04_deciduous_broadleaf, k = 5) + s(pland_05_mixed_forest, 
#>     k = 5) + s(pland_12_cropland, k = 5) + s(pland_13_urban, 
#>     k = 5) + protocol_type + s(time_observations_started, bs = "cc", 
#>     k = 7)


# explicitly specify where the knots should occur for time_observations_started
# this ensures that the cyclic spline joins the variable at midnight
# this won't happen by default if there are no data near midnight
time_knots <- list(time_observations_started = seq(0, 24, length.out = k_time))

# Specify the error distribution/link function to be used in our GAM model, then fit
# the model using 'gam'. This dist should be similar to the distribution of our 
# response variable (without mean) i.e. in a linear regression estimate we assume 
# a normal distribution of points above/below our line

# From the histogram plots obtained previously we identify:
# negative binomial - useful for data exhibiting 'over-dispersion' that is very 
# common in ecological count data
m_nb <- gam(gam_formula,
            data = ebird_split, 
            family = "nb",
            knots = time_knots)

# NOTE: In general it seems that the Negative Binomial gives the best result. 
# However, should we end up studying only a single species, we should revisit all
# 3 distributions suggested to see if it is the best in the specific case

#### end ####


#### Model Critique ####

# It is recommended for us to check how the observed count predictions vary with
# each covariate - if the relationship is 'too wiggly' to be biologically realistic
# then it is suggested that the k value is lowered.

# Write a ggplot function
plot_gam <- function(m, title = NULL) {
  
  # capture plot
  png(file = file.path(data_folder, "analytics", "model_covariates_dependency.png"))
  p <- plot(m, pages = 1)
  
  # extract data
  p_df <- map_df(p, ~ tibble(cov = rep(.$xlab, length(.$x)),
                             x = .$x, fit = .$fit, se = .$se))
  
  # plot
  g <- ggplot(p_df) +
    aes(x = x, y = fit,
        ymin = fit - se, ymax = fit + se) +
    geom_ribbon(fill = "grey80") +
    geom_line(col = "blue") +
    facet_wrap(~ cov, scales = "free_x") +
    labs(x = NULL,
         y = "Smooth function",
         title = title)
  print(g)
  invisible(p_df)
  dev.off()
}

# Call the plot function for the NB distribution
plot_gam(m_nb, title = "Negative Binomial GAM")

# Finally, given the plots for the NB function select the
# Negative Binomial model as the one to use 
pred_model <- m_nb

#### end ####


#### Abundance Prediction Mapping ####

# So far we have used eBird data to create an abundance model that can predict the 
# count number, given all the effort and habitat covariates we have used. Given that
# we want to predict over a whole region: we already made a prediction surface for 
# the area of interest that contains landcover types (using MODIS). However, the effort
# covariates do not exist in this way - they are dependent on an individual outing. 
# Instead, to make predictions over the surface, we will assume that for any point
# the effort covariates take the values given by a 'standard checklist': 1km travelling,
# 1 hour outing, 1 person. The final component is to find the optimal time of day to start 
# a checklist (results in most accurate abundance predictions).

# Create a df of covariates that take their mean values, or are 'standard' checklist
# values but have a range of start times throughout the day
seq_tod <- seq(0, 24, length.out = 300)
tod_df <- ebird_split %>% 
  select(starts_with("pland")) %>% 
  summarize_all(mean, na.rm = TRUE) %>% 
  ungroup() %>% 
  # Can change date depending on your time of interest
  mutate(day_of_year = yday(ymd(str_glue("{max_lc_year}-06-15"))),
         duration_minutes = 60,
         effort_distance_km = 1,
         number_observers = 1,
         protocol_type = "Traveling") %>% 
  cbind(time_observations_started = seq_tod)

# Predict at the different start times
pred_tod <- predict(pred_model, newdata = tod_df, 
                    type = "link", 
                    se.fit = TRUE) %>% 
  as_tibble() %>% 
  # calculate backtransformed confidence limits
  transmute(time_observations_started = seq_tod,
            pred = pred_model$family$linkinv(fit),
            pred_lcl = pred_model$family$linkinv(fit - 1.96 * se.fit),
            pred_ucl = pred_model$family$linkinv(fit + 1.96 * se.fit))

# Optimal time of day is that with max lower bound of 95% conf interval (deals with if 
# limited data, we select a time we are confident in)
t_peak <- pred_tod$time_observations_started[which.max(pred_tod$pred_lcl)]

# Plot these predicted abundance values and conf interval to visualise 
png(file = file.path(data_folder, "analytics", "optimal_time_of_day.png"))
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

# Now we have everything we need to add effort covariates to the prediction surface
# giving us a complete prediction surface with all the covariates used in our model
pred_surface_eff <- pred_surface %>% 
  # Can change date depending on your time of interest
  mutate(day_of_year = yday(ymd(str_glue("{max_lc_year}-", time_of_year))),
         time_observations_started = t_peak,
         duration_minutes = 60,
         effort_distance_km = 1,
         number_observers = 1,
         protocol_type = "Traveling")


# Make predictions over the complete surface 
pred <- predict(pred_model, newdata = pred_surface_eff, 
                type = "link", 
                se.fit = TRUE) %>% 
  as_tibble() %>% 
  # calculate confidence limits and back transform (use inverse link function and 
  # apply to the link fit we have just obtained)
  transmute(abd = pred_model$family$linkinv(fit),
            abd_se = pred_model$family$linkinv(se.fit),
            abd_lcl = pred_model$family$linkinv(fit - 1.96 * se.fit),
            abd_ucl = pred_model$family$linkinv(fit + 1.96 * se.fit)) %>%
  # add to lat/lon from prediction surface
  bind_cols(pred_surface_eff, .) %>% 
  select(latitude, longitude, abd, abd_se, abd_lcl, abd_ucl)

# Now we want to convert these predictions into a raster file
r_pred <- pred %>% 
  # convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  select(abd, abd_se) %>% 
  st_transform(crs = projection(r)) %>% 
  rasterize(r)
# Not sure why we remove a layer?
r_pred <- r_pred[[-1]]

# save the rasters
writeRaster(r_pred[["abd"]], 
            filename = file.path(data_folder, "output data", "abundance_map.tif"),
            overwrite = TRUE)
writeRaster(r_pred[["abd_se"]], 
            filename = file.path(data_folder, "output data", "abundance_map_uncertainty.tif"), 
            overwrite = TRUE)

# Now lets plot the results. The values obtained on this map are the expected number 
# of species seen by an average eBirder conducting a 1 hour, 1 km checklist for 
# which counting started at the optimal time. As detectability is not perfect, we expect
# true species abundance to be higher than these values, but without estimating 
# the detection rate directly it’s difficult to say how much higher. However, the 
# relative sizes should be okay (hence the name relative abundance)

# Change the CRS of the raster to match that of the GIS data for plotting
r_pred_proj <- projectRaster(r_pred, crs = 4326, method = "ngb")

# Plot the data for each layer (abd and abd_se)
png(file = file.path(data_folder, "analytics", "abundance_maps.png"))
par(mfrow=c(1,2))
for (nm in names(r_pred)) {
  r_plot <- r_pred_proj[[nm]]
  
  # Set the plot area and add GIS components
  par(mar = c(4.5, 0.25, 0.25, 0.25))
  plot(ne_land, col = "#dddddd", border = "#888888", lwd = 0.5)
  
  # Modified plasma (colour) palette
  plasma_rev <- rev(plasma(25, end = 0.9))
  gray_int <- colorRampPalette(c("#dddddd", plasma_rev[1]))
  pal <- c(gray_int(4)[2], plasma_rev)
  
  # Setup plot characetristics depending on whether abd or abd_se
  if (nm == "abd") {
    title <- paste(species_name, "Relative Abundance")
    # set very low values to zero
    r_plot[r_plot <= zero_threshold] <- NA
    # log transform
    r_plot <- log10(r_plot)
    # breaks and legend
    mx <- ceiling(100 * cellStats(r_plot, max)) / 100
    mn <- floor(100 * cellStats(r_plot, min)) / 100
    brks <- seq(mn, mx, length.out = length(pal) + 1)
    lbl_brks <- sort(c(-2:2, mn, mx))
    lbls <- round(10^lbl_brks, 2)
  } else {
    title <- paste(species_name, "Abundance Uncertainty (SE)")
    # breaks and legend
    mx <- ceiling(1000 * cellStats(r_plot, max)) / 1000
    mn <- floor(1000 * cellStats(r_plot, min)) / 1000
    brks <- seq(mn, mx, length.out = length(pal) + 1)
    lbl_brks <- seq(mn, mx, length.out = 5)
    lbls <- round(lbl_brks, 2)
  }
  
  # Plot the data
  plot(r_plot, 
       col = pal, breaks = brks, 
       maxpixels = ncell(r_plot),
       legend = FALSE, add = TRUE)
  
  # Add country lines for GIS
  #plot(ne_country_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
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
                                cex = 1, line = 0))
}
dev.off()

#### end ####
