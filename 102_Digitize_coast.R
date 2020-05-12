
#
# Make 'coast' used in script 103
# 'coast' is a segmented line along the Norwegian coast
# See graphs in 103
#
# Needs only to be done once, can just use thesaved result thereafter
#

library(ggplot2)
library(dplyr)

library(sp)
crs_longlat <- "+proj=longlat +ellps=WGS84 +datum=WGS84"
crs_utm <- "+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m"

if (FALSE){
  # Click to mark locations of "coastal current"
  windows()
  maps::map("world", "Norway", exact = TRUE)   # mainland only
  X <- locator()
  # saveRDS(X, "Data/102_coast_coordinates_list.rmd")
  # X <- readRDS("Data/102_coast_coordinates.rmd")
  
  coast <- data.frame(Longitude = X$x, Latitude = X$y)

  # note that you might want to change zone=32 to e.g. 33
  # - from Long,Lat to UTM zone 32
  SP <- SpatialPoints(coast[,c("Longitude", "Latitude")],
                      proj4string=CRS(crs_longlat)
  )
  SP.UTM <- spTransform(SP, CRS(crs_utm))
  # Add transformed coords to data set
  coast$x <- SP.UTM@coords[,1]
  coast$y <- SP.UTM@coords[,2]
  
  # saveRDS(coast, "Data/102_coast_coordinates.rmd")

  }

