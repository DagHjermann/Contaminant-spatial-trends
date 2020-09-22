

crs_longlat <- "+proj=longlat +ellps=WGS84 +datum=WGS84"
crs_utm <- "+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m"


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
if (FALSE){
  i <- 6
  p <- list(x = 500000, y = 6550000)
  q <- coordinate_on_segment(i, p, coast)
}


get_distance_to_point <- function(point1, point2){
  dx <- (point1$x - point2$x)/1000
  dy <- (point1$y - point2$y)/1000
  sqrt(dx^2 + dy^2)
}


get_distance_along_coast_all <- function(point, df_segments, df_segments_distances){
  i <- 1:(nrow(df_segments)-1)
  check_segments <- purrr::map(i, coordinate_on_segment, point = point, df_segments = df_segments)
  on_segment <- purrr::transpose(check_segments)$on_segment %>% unlist()
  # result1 = all segments which have a normal line that hits our point 
  # I.e. closest point along the segment line(s)
  # Note: if 'point' is "outside" an concave edge, the normal line 
  #   from no segments may hits our point)
  # Note 2: may also be more than one point; in get_distance_along_coast() we select the closest one
  if (sum(on_segment, na.rm = TRUE) > 0){
    # Usual case: 'the'point' is along the normal line from at least one segment 
    i_sel <- i[on_segment]
    coord <- check_segments[i_sel] %>% purrr::transpose() %>% .$coord %>% purrr::map_df(~as.data.frame(.))
    dx1 <- (coord$x - point$x)/1000
    dy1 <- (coord$y - point$y)/1000
    dx2 <- (coord$x - coast$x[i_sel])/1000
    dy2 <- (coord$y - coast$y[i_sel])/1000
    result1 <-   list(
      coord = coord,
      segment_no = i_sel,
      distance_to_segment = sqrt(dx1^2 + dy1^2),
      distance_from_segment_start = sqrt(dx2^2 + dy2^2),
      distance = sqrt(dx2^2 + dy2^2) + df_segments_distances[i_sel]
    )
  } 
  # Then also find closest corner of df_segments (result2)
  # I.e. closest corner of the segment line(s)
  dist <- 1:nrow(df_segments) %>% map_dbl(~get_distance_to_point(df_segments[.,], point))
  i_sel <- which.min(dist)
  result2 <-   list(
    coord = df_segments[i_sel, c("x", "y")],
    segment_no = i_sel,
    distance_to_segment = dist[i_sel],
    distance_from_segment_start = 0,
    distance = df_segments_distances[i_sel]
  )
  # Same if as above: if points along segment li
  if (sum(on_segment, na.rm = TRUE) > 0){
    result <- list(
      coord = bind_rows(result1$coord, 
                        result2$coord),
      segment_no = c(result1$segment_no, 
                     result2$segment_no),
      distance_to_segment = c(result1$distance_to_segment, 
                              result2$distance_to_segment),
      distance_from_segment_start = c(result1$distance_from_segment_start, 
                                      result2$distance_from_segment_start),
      distance = c(result1$distance, result2$distance)
    )
  } else {
    result <- result2
  }
  result
}

# debugonce(get_distance_along_coast)

# TEST
if (FALSE){
  p <- list(x = 475000, y = 6450000)  # usual case
  p <- list(x = 500000, y = 6550000)  # usual case
  p <- list(x = 1285975, y = 7880861) # unusual case: no normal lines hits our point
  get_distance_along_coast_all(point = p, df_segments = coast, df_segments_distances = coastsegment_distance)
  
}

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

# TEST
if (FALSE){
  p <- list(x = 475000, y = 6450000)
  p <- list(x = 500000, y = 6550000)
  get_distance_along_coast(point = p, df_segments = coast, df_segments_distances = coastsegment_distance)
}

check_distance_along_coast <- function(point, 
                                       df_segments, 
                                       df_segments_distances, 
                                       df_map = mapdata){
  get_distance_along_coast_s <- safely(get_distance_along_coast)
  result <- get_distance_along_coast_s(point = point, df_segments = df_segments, df_segments_distances = df_segments_distances)
  gg <- ggplot(df_map, aes(x, y)) + 
    geom_path() +
    coord_fixed(xlim = point$x + c(-100000, 100000), ylim = point$y + c(-100000, 100000)) +
    geom_path(data = df_segments, color = "blue") +
    geom_point(data = df_segments, color = "blue") +
    geom_point(data = data.frame(point), color = "red")
  if (is.null(result$error)){
    result <- result$result
    i <- result$segment_no
    gg <- gg + 
      geom_point(data = df_segments[i:(i+1),], color = "blue", shape = 1, size = 4) +
      geom_point(data = result$coord, color = "green3")
    print(result)
  } else {
    
  }
  print(gg)
  invisible(gg)
}


# TEST
if (FALSE){
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
}



#
# Give distance along the coastline, get the coordinates in return  
# Only works for a single distance  
# Note that the '
#
get_point_on_coastline <- function(distance = 1000, 
                                   df_segments = coast,
                                   df_segments_distances = coastsegment_distance){
  # First, which segment are we on;
  segment_no <- which(df_segments_distances < distance) %>% tail(1)
  if (length(segment_no) > 0){
    # Distance from start of segment to the coordinate
    distanceleft <- distance - df_segments_distances[segment_no]
    # SAme, as fraction of the distance length
    distanceleft_fraction <- distanceleft/diff(df_segments_distances[segment_no:(segment_no+1)])
    # Start and end coordinates of the segment
    segment_start <- df_segments[segment_no, c("x", "y")]
    segment_end <- df_segments[segment_no + 1, c("x", "y")]
    # The coordinate we want
    x <- segment_start["x"] + distanceleft_fraction*(segment_end["x"] - segment_start["x"])
    y <- segment_start["y"] + distanceleft_fraction*(segment_end["y"] - segment_start["y"])
    # Include 'distance' in the result (for use in map_df)
    result <- data.frame(distance = distance, x = as.numeric(x), y = as.numeric(y))
  } else {
    # If the input is distance = 0 (start of the 'df_segments' data)
    result <- data.frame(distance = distance, x = df_segments$x[1], y = df_segments$y[1])
  }
  result
}

# Test
if (FALSE){
  # debugonce(get_point_on_coastline)
  get_point_on_coastline(1000)
  get_point_on_coastline(0) %>% str()
}
