# Hamish Campbell
# Load, plot and export threat maps for birds as raster files
# Uses the work by: https://www.nature.com/articles/s41559-021-01542-9

# Load the required packages
library(sf)
library(tidyverse)
library(viridis)

# Run checks to see if the relevant output directories exist
if (!dir.exists("output data")){
  dir.create("output data")
}

if (!dir.exists("analytics")){
  dir.create("analytics")
}

# Load the threat data and keep only non-habitat threats for birds
threats <- st_read("input data/Grid_mammal_amph_threat_predictions.shp")  # type=1
threats_birds <- select(threats, c("B5_1","B8_1","B9","B11","geometry"))

#### Save a plot of each threat map ####

# Hunting
setEPS()
postscript("analytics/hunting.eps")
plot(threats_birds["B5_1"], border = 'transparent', 
        pal=plasma, nbreaks=8, key.pos=4, main=NULL, key.length=lcm(8))
dev.off()

# Invasives
setEPS()
postscript("analytics/invasives.eps")
plot(threats_birds["B8_1"], border = 'transparent', 
     pal=plasma, nbreaks=8, key.pos=4, main=NULL, key.length=lcm(8))
dev.off()

# Pollution
setEPS()
postscript("analytics/pollution.eps")
plot(threats_birds["B9"], border = 'transparent', 
     pal=plasma, nbreaks=8, key.pos=4, main=NULL, key.length=lcm(8))
dev.off()

# Climate Change
setEPS()
postscript("analytics/climate_change.eps")
plot(threats_birds["B11"], border = 'transparent', 
     pal=plasma, nbreaks=8, key.pos=4, main=NULL, key.length=lcm(8))
dev.off()

#### end ####

# Convert CRS 
threats_birds <- st_transform(threats_birds, crs="EPSG:4326")

# Rasterise and save as raster (stack if possible)





