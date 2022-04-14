library(ebirdst)
library(raster)
library(sf)
library(fields)
library(rnaturalearth)

# Store the Access Key obtained from the eBird website
# THIS REQUIRES FILLING IN BY USER
set_ebirdst_access_key("XXXXX")

# download example data, yellow-bellied sapsucker in michigan
dl_path <- ebirdst_download(species = "example_data")
# load relative abundance raster stack with 52 layers, one for each week
abd <- load_raster("abundance", path = dl_path)

# load species specific mapping parameters
pars <- load_fac_map_parameters(dl_path)
# custom coordinate reference system
crs <- pars$custom_projection
# legend breaks
breaks <- pars$abundance_bins
# legend labels for top, middle, and bottom
labels <- attr(breaks, "labels")

# get a date vector specifying which week each raster layer corresponds to
weeks <- parse_raster_dates(abd)
print(weeks)
#>  [1] "2019-01-04" "2019-01-11" "2019-01-18" "2019-01-25" "2019-02-01" "2019-02-08"
#>  [7] "2019-02-15" "2019-02-22" "2019-03-01" "2019-03-08" "2019-03-15" "2019-03-22"
#> [13] "2019-03-29" "2019-04-05" "2019-04-12" "2019-04-19" "2019-04-26" "2019-05-03"
#> [19] "2019-05-10" "2019-05-17" "2019-05-24" "2019-05-31" "2019-06-07" "2019-06-14"
#> [25] "2019-06-21" "2019-06-28" "2019-07-06" "2019-07-13" "2019-07-20" "2019-07-27"
#> [31] "2019-08-03" "2019-08-10" "2019-08-17" "2019-08-24" "2019-08-31" "2019-09-07"
#> [37] "2019-09-14" "2019-09-21" "2019-09-28" "2019-10-05" "2019-10-12" "2019-10-19"
#> [43] "2019-10-26" "2019-11-02" "2019-11-09" "2019-11-16" "2019-11-23" "2019-11-30"
#> [49] "2019-12-07" "2019-12-14" "2019-12-21" "2019-12-28"

# select a week in the middle of the year
abd <- abd[[26]]

# project to species specific coordinates
# the nearest neighbor method preserves cell values across projections
abd_prj <- projectRaster(abd, crs = crs, method = "ngb")

# get reference data from the rnaturalearth package
# the example data currently shows only the US state of Michigan
wh_states <- ne_states(country = c("United States of America", "Canada"),
                       returnclass = "sf") %>% 
  st_transform(crs = crs) %>% 
  st_geometry()

# start plotting
par(mfrow = c(1, 1), mar = c(0, 0, 0, 0))

# use raster bounding box to set the spatial extent for the plot
bb <- st_as_sfc(st_bbox(trim(abd_prj)))
plot(bb, col = "white", border = "white")
# add background reference data
plot(wh_states, col = "#cfcfcf", border = NA, add = TRUE)

# plot zeroes as light gray
plot(abd_prj, col = "#e6e6e6", maxpixels = ncell(abd_prj),
     axes = FALSE, legend = FALSE, add = TRUE)

# define color palette
pal <- abundance_palette(length(breaks) - 1, "weekly")
# plot abundance
plot(abd_prj, col = pal, breaks = breaks, maxpixels = ncell(abd_prj),
     axes = FALSE, legend = FALSE, add = TRUE)

# state boundaries
plot(wh_states, add = TRUE, col = NA, border = "white", lwd = 1.5)

# legend
label_breaks <- seq(0, 1, length.out = length(breaks))
image.plot(zlim = c(0, 1), breaks = label_breaks, col = pal,
           smallplot = c(0.90, 0.93, 0.15, 0.85),
           legend.only = TRUE,
           axis.args = list(at = c(0, 0.5, 1), 
                            labels = round(labels, 2),
                            cex.axis = 0.9, lwd.ticks = 0))


# MY TESTS
print(ebirdst_runs)  # this lists all species available and associated info

