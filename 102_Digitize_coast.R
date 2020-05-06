
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

#
# Coordinates along coast
#
coast <- readRDS("Data/102_coast_coordinates.rmd")

#
# Get Norway map data
#
test <- maps::map("world", "Norway", plot = FALSE)   # map data for Norway - this is just to get region names
sel <- grepl("Svalbard", test$names) | test$names == "Norway:Jan Mayen"  # select Svalbard + Jan Mayen
# test$names[!sel]
map <- maps::map("world", test$names[!sel], exact = TRUE, plot = FALSE)  # Norway w/o Svalbard + Jan Mayen
runs <- rle(x)
map_df <- data.frame(Longitude = map$x, Latitude = map$y)


#
# Test plot (long-lat)
#

ggplot(map_df, aes(Longitude, Latitude)) + 
  geom_path(size = 0.4)


#
# Add UTM coordinates (x and y) to map
#
coordinate_exists <- !is.na(map_df$Longitude)   # sp doesn't like NAs
SP <- SpatialPoints(map_df[coordinate_exists, c("Longitude", "Latitude")],
                    proj4string=CRS(crs_longlat)
)
SP.UTM <- spTransform(SP, CRS(crs_utm))
# Add transformed coords to data set
map_df$x[coordinate_exists] <- SP.UTM@coords[,1]
map_df$y[coordinate_exists] <- SP.UTM@coords[,2]

# Plot map 
ggplot(map_df, aes(x, y)) + 
  geom_path() +
  coord_fixed()


#
# Find point on segment i where a perpendicular
#   line from 'point' ends up
#

# From https://stackoverflow.com/a/12499474/1734247
coordinate_on_segment <- function(i, point, df_segments){
  x1 <- df_segments$x[i]
  x2 <- df_segments$x[i+1]
  y1 <- df_segments$y[i]
  y2 <- df_segments$y[i+1]
  px = x2-x1
  py = y2-y1
  dAB = px*px + py*py
  u = ((point$x - x1) * px + (point$y - y1) * py) / dAB
  coord <- data.frame(x = x1 + u * px, 
             y = y1 + u * py)
  on_segment <- coord$x >= min(df_segments$x[i:(i+1)]) &
    coord$x <= max(df_segments$x[i:(i+1)]) &
    coord$y >= min(df_segments$y[i:(i+1)]) &
    coord$y <= max(df_segments$y[i:(i+1)])
  list(coord = coord, 
       on_segment = on_segment)
}

#
# Test
#
i <- 6
p <- list(x = 500000, y = 6550000)
q <- coordinate_on_segment(i, p, coast)

#
# Plot result
#
gg <- ggplot(map_df, aes(x, y)) + 
  geom_path() +
  coord_fixed(xlim = c(400000, 650000), ylim = c(6400000, 6700000))
# gg

gg + 
  geom_path(data = coast, color = "blue") +
  geom_point(data = coast, color = "blue") +
  geom_point(data = coast[i:(i+1),], color = "blue", shape = 1, size = 4) +
  geom_point(data = data.frame(p), color = "red") +
  geom_point(data = q$coord, color = "green3")

q$coord
dx <- (q$coord$x - coast$x[i])/1000
dy <- (q$coord$y - coast$y[i])/1000
sqrt(dx^2 + dy^2)

#
# Distances along coast segments
#
segment_dx <- diff(coast$x)/1000
segment_dy <- diff(coast$y)/1000
segment_distance <- sqrt(segment_dx^2 + segment_dy^2) %>% cumsum()
segment_distance <- c(0, segment_distance)
segment_distance

#
# Find best segment and find distance from start of that segment
#
i <- 5:7
check_segments <- purrr::map(i, coordinate_on_segment, point = p, df_segments = coast)
on_segment <- purrr::transpose(check_segments)$on_segment %>% unlist()
i_sel <- i[on_segment]
i_sel
coord <- check_segments[i_sel] %>% purrr::transpose() %>% .$coord %>% .[[1]]
coord
dx <- (coord$x - coast$x[i_sel])/1000
dy <- (coord$y - coast$y[i_sel])/1000
sqrt(dx^2 + dy^2)

#
# Distances for start points of the coast segments
#
segment_dx <- diff(coast$x)/1000
segment_dy <- diff(coast$y)/1000
coastsegment_distance <- sqrt(segment_dx^2 + segment_dy^2) %>% cumsum()
coastsegment_distance <- c(0, coastsegment_distance)
coastsegment_distance

get_distance_along_coast_all <- function(point, df_segments, df_segments_distances){
  i <- 1:(nrow(df_segments)-1)
  check_segments <- purrr::map(i, coordinate_on_segment, point = point, df_segments = df_segments)
  on_segment <- purrr::transpose(check_segments)$on_segment %>% unlist()
  # if (sum(on_segment, na.rm = TRUE) > 1){
  #   stop("More than one segment found")
  if (sum(on_segment, na.rm = TRUE) == 0){
    stop("No segments found")
  }
  i_sel <- i[on_segment]
  coord <- check_segments[i_sel] %>% purrr::transpose() %>% .$coord %>% purrr::map_df(~as.data.frame(.))
  dx1 <- (coord$x - point$x)/1000
  dy1 <- (coord$y - point$y)/1000
  dx2 <- (coord$x - coast$x[i_sel])/1000
  dy2 <- (coord$y - coast$y[i_sel])/1000
  list(
    coord = coord,
    segment_no = i_sel,
    distance_to_segment = sqrt(dx1^2 + dy1^2),
    distance_from_segment_start = sqrt(dx2^2 + dy2^2),
    distance = sqrt(dx2^2 + dy2^2) + df_segments_distances[i_sel]
  )
}

# debugonce(get_distance_along_coast)
# Test
p <- list(x = 475000, y = 6450000)
p <- list(x = 500000, y = 6550000)
get_distance_along_coast_all(point = p, df_segments = coast, df_segments_distances = coastsegment_distance)

get_distance_along_coast <- function(point, df_segments, df_segments_distances){
  result <- get_distance_along_coast_all(point = point, df_segments = df_segments, df_segments_distances = df_segments_distances)
  i <- which.min(result$distance_to_segment)  # closest segment
  list(
    coord = result$coord[i,],
    segment_no = result$segment_no[i],
    distance_to_segment = result$distance_to_segment[i],
    distance_from_segment_start = result$distance_from_segment_start[i],
    distance = result$distance[i]
  )
  }

p <- list(x = 475000, y = 6450000)
p <- list(x = 500000, y = 6550000)
get_distance_along_coast(point = p, df_segments = coast, df_segments_distances = coastsegment_distance)

check_distance_along_coast <- function(point, df_segments, df_segments_distances, mapdata = map_df){
  result <- get_distance_along_coast(point = point, df_segments = df_segments, df_segments_distances = df_segments_distances)
  i <- result$segment_no
  gg <- ggplot(mapdata, aes(x, y)) + 
    geom_path() +
    coord_fixed(xlim = point$x + c(-100000, 100000), ylim = point$y + c(-100000, 100000)) +
    geom_path(data = df_segments, color = "blue") +
    geom_point(data = df_segments, color = "blue") +
    geom_point(data = df_segments[i:(i+1),], color = "blue", shape = 1, size = 4) +
    geom_point(data = data.frame(point), color = "red") +
    geom_point(data = result$coord, color = "green3")
  print(result)
  print(gg)
  invisible(gg)
}

check_distance_along_coast(point = list(x = 500000, y = 6550000), 
                           df_segments = coast, 
                           df_segments_distances = coastsegment_distance)

check_distance_along_coast(point = list(x = 475000, y = 6450000), 
                           df_segments = coast, 
                           df_segments_distances = coastsegment_distance)

check_distance_along_coast(point = list(x = 500000, y = 6450000), 
                           df_segments = coast, 
                           df_segments_distances = coastsegment_distance)

check_distance_along_coast(point = list(x = 300000, y = 6550000), 
                           df_segments = coast, 
                           df_segments_distances = coastsegment_distance)

check_distance_along_coast(point = list(x = 600000, y = 6550000), 
                           df_segments = coast, 
                           df_segments_distances = coastsegment_distance)

check_distance_along_coast(point = list(x = 600000, y = 6550000), 
                           df_segments = coast, 
                           df_segments_distances = coastsegment_distance)
