---
title: "103_Distance_along_coast"
author: "DHJ"
date: "12 5 2020"
output:
  html_document:
    toc: true    
    toc_float: true
    keep_md: true  
---

## Make station data  
Makes station data with **distance along coast* (Dist_along_coast) as one of the variables.   
* The station data is based on the ICES station dictionary   
* Also includes MSTAT, with codes for station type:   
    + IH: impacted by harmful substances (harbours, industry sites)   
    + RH: representative stations     
    + B: background stations   
    

## 1. Packages + definitions  

```r
library(ggplot2)
library(dplyr)
library(purrr)
library(sp)
library(readxl)
library(safejoin)   # github hrbrmstr

crs_longlat <- "+proj=longlat +ellps=WGS84 +datum=WGS84"
crs_utm <- "+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m"
```

Functions  

```r
source("103_Distance_along_coast_functions.R")
```


## 2. Data

### Station data  
From Kartbase.xlsx   

```r
# For check
df_stations <- read_excel(
  "Input_data/Kartbase.xlsx",
  range = "$AM$1:$AP$72") %>%
  rename()
colnames(df_stations) <- c("STATION_CODE", "Lat", "Lon", "Station_Name")
```

### ICES stations   
For MSTAT  

```r
df_stations_ices <- readRDS("Data/101_Selected_stations.rds") %>%
  group_by(STATION_CODE) %>%
  mutate(StartYear_last = max(StartYear)) %>%
  filter(StartYear_last == StartYear)

# Check that we have all stations
n1 <- unique(readRDS("Data/101_Selected_stations.rds")$STATION_CODE)
n2 <- unique(df_stations_ices$STATION_CODE)
cat("Number of stations lost in filtering:", sum(!n1 %in% n2), "\n")
```

```
## Number of stations lost in filtering: 0
```
### Add MSTAT  

```r
df_stations <- df_stations %>%
  filter(!grepl("G", STATION_CODE)) %>%
  safe_left_join(
    df_stations_ices %>% select(STATION_CODE, MSTAT),
    by = "STATION_CODE",
    na_matches = "never",
    check = "bCV"
  ) %>%
  mutate(MSTAT = case_when(
    STATION_CODE == "19B" ~ "B",
    STATION_CODE == "19N" ~ "B",
    STATION_CODE == "97A3" ~ "IH",
    STATION_CODE == "28A2" ~ "RH",
    TRUE ~ MSTAT)
  )

cat("Stations lacking MSTAT:", sum(is.na(df_stations$MSTAT)), "\n")  
```

```
## Stations lacking MSTAT: 0
```


### Add stations


```r
df_stations_extra <- read.csv(textConnection("
STATION_CODE, Station_Name, Lat, Lon, MSTAT,
24B, 24B Bergen harbour, 60.39664, 5.27069, IH
36A1, 36A1 Færder, 59.07357, 10.42522, RH
28A2, 28A2 Vågsvåg, 61.93622,  5.048783333, RH
97A3, 97A3 Bodø harbour, 67.2963, 14.3956, IH
I714, I714 Sylterøya (Langesundfjord), 59.0514, 9.7038, IH
"),
stringsAsFactors = FALSE
)

# OLD  
# df_stations <- df_stations %>%
#   bind_rows(df_stations_extra)
```

### Coordinates along the coastal current  
Made in script 102  

```r
coast <- readRDS("Data/102_coast_coordinates.rmd")
```

### Norway map data

```r
#
# Get Norway map data
#
test <- maps::map("world", "Norway", plot = FALSE)   # map data for Norway - this is just to get region names
sel <- grepl("Svalbard", test$names) | test$names == "Norway:Jan Mayen"  # select Svalbard + Jan Mayen
# test$names[!sel]
map <- maps::map("world", test$names[!sel], exact = TRUE, plot = FALSE)  # Norway w/o Svalbard + Jan Mayen
mapdata <- data.frame(Longitude = map$x, Latitude = map$y)
```


### Plot coordinates along coast  
Long-lat coordinates  

```r
ggplot(mapdata, aes(Longitude, Latitude)) + 
  geom_path() +
  coord_map("lambert", lat0 = 64, lon0 = 11) +
  geom_path(data = coast, color = "blue") +
  geom_point(data = coast, color = "blue") 
```

![](103_Distance_along_coast_files/figure-html/unnamed-chunk-9-1.png)<!-- -->



## 4. Add UTM coordinates (x and y) to map

```r
#
# Add UTM coordinates (x and y) to map
#
coordinate_exists <- !is.na(mapdata$Longitude)   # sp doesn't like NAs
SP <- SpatialPoints(mapdata[coordinate_exists, c("Longitude", "Latitude")],
                    proj4string=CRS(crs_longlat)
)
SP.UTM <- spTransform(SP, CRS(crs_utm))
# Add transformed coords to data set
mapdata$x[coordinate_exists] <- SP.UTM@coords[,1]
mapdata$y[coordinate_exists] <- SP.UTM@coords[,2]

# Plot map 
ggplot(mapdata, aes(x, y)) + 
  geom_path() +
  coord_fixed()
```

![](103_Distance_along_coast_files/figure-html/unnamed-chunk-10-1.png)<!-- -->


## 5. Distance from point along coast  

### Make 'coastsegment_distance'   
Distances along coast segments

```r
#
# Distances for start points of the coast segments
#
segment_dx <- diff(coast$x)/1000
segment_dy <- diff(coast$y)/1000
coastsegment_distance <- sqrt(segment_dx^2 + segment_dy^2) %>% cumsum()
coastsegment_distance <- c(0, coastsegment_distance)
coastsegment_distance
```

```
##  [1]    0.00000   35.57563   69.07444   77.13420  104.65681  139.71662
##  [7]  201.95271  261.89173  333.57830  386.84844  480.28367  566.00644
## [13]  642.95392  723.92513  794.51827  870.45524  939.27022 1063.94254
## [19] 1170.08556 1267.64526 1499.11265 1639.76846 1736.21744 1818.39978
## [25] 1882.21931 1982.03134 2054.87673 2142.02076 2215.68761 2294.55571
## [31] 2371.42529 2444.72861 2496.01659 2561.46976 2612.91895 2665.37530
## [37] 2685.72980
```

### Show some "fixed km" along the coast on map    
Using function get_point_on_coastline()   

```r
# Get positions for these km's:
points <- 
  c(0, 500, 1000, 1500, 2000, 2500, 2685) %>% map_df(~get_point_on_coastline(.))

# Direction of text labels:
points$Text_direction <- "West"
points$Text_direction[c(1,6,7)] <- "East"

# PLot
ggplot(mapdata, aes(x, y)) +
  geom_path() +
  coord_fixed() +
  geom_path(data = coast, color = "red") +
  geom_point(data = points, color = "blue") +
  geom_text(data = points %>% filter(Text_direction == "West"), 
            aes(x = x - 20000, label = paste(distance, "km")), color = "blue", hjust = 1) +
  geom_text(data = points %>% filter(Text_direction == "East"), 
            aes(x = x + 20000, label = paste(distance, "km")), color = "blue", hjust = 0) +
  expand_limits(x = c(min(mapdata$x, na.rm = TRUE) - 150000,
                      max(mapdata$x, na.rm = TRUE) + 200000)) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )
```

![](103_Distance_along_coast_files/figure-html/unnamed-chunk-12-1.png)<!-- -->

### Plot 

```r
ggplot(mapdata, aes(x, y)) + 
  geom_path() +
  coord_fixed() +
  geom_path(data = coast, color = "red")
```

![](103_Distance_along_coast_files/figure-html/unnamed-chunk-13-1.png)<!-- -->


### Test 'get_distance_along_coast' and 'check_distance_along_coast'  

```r
# Usual case: the point is on the normal to at least one of the edges
check_distance_along_coast(point = list(x = 600000, y = 6550000), 
                           df_segments = coast, 
                           df_segments_distances = coastsegment_distance,
                           df_map = mapdata)
```

```
## $coord
##          x       y
## 1 610086.5 6559883
## 
## $segment_no
## [1] 2
## 
## $distance_to_segment
## [1] 14.12165
## 
## $distance_from_segment_start
## [1] 10.59748
## 
## $distance
## [1] 46.1731
```

![](103_Distance_along_coast_files/figure-html/unnamed-chunk-14-1.png)<!-- -->

```r
# Unusual case: the point is by a "concave" edge on the coast
check_distance_along_coast(point = list(x = 1285975, y = 7880861), 
                           df_segments = coast, 
                           df_segments_distances = coastsegment_distance)
```

```
## $coord
##          x       y
## 36 1341902 7887778
## 
## $segment_no
## [1] 36
## 
## $distance_to_segment
## [1] 56.35326
## 
## $distance_from_segment_start
## [1] 0
## 
## $distance
## [1] 2665.375
```

![](103_Distance_along_coast_files/figure-html/unnamed-chunk-14-2.png)<!-- -->

```r
if (FALSE){
  check_distance_along_coast(point = list(x = 500000, y = 6550000), 
                             df_segments = coast, 
                             df_segments_distances = coastsegment_distance)
  
  check_distance_along_coast(point = list(x = 500000, y = 6450000), 
                             df_segments = coast, 
                             df_segments_distances = coastsegment_distance)
  
}
```


## 6. Add Dist_along_coast to station data  

### Add UTM coordinates (x,y)  

```r
df_stations <- as.data.frame(df_stations)

coordinate_exists <- !is.na(df_stations$Lon)   # sp doesn't like NAs
SP <- SpatialPoints(df_stations[coordinate_exists, c("Lon", "Lat")],
                    proj4string=CRS(crs_longlat)
)
SP.UTM <- spTransform(SP, CRS(crs_utm))
# Add transformed coords to data set
df_stations$x[coordinate_exists] <- SP.UTM@coords[,1]
df_stations$y[coordinate_exists] <- SP.UTM@coords[,2]
```

### Plot  

```r
ggplot(mapdata, aes(x, y)) + 
  geom_path() +
  coord_fixed() +
  geom_path(data = coast, color = "blue") +
  geom_point(data = df_stations, color = "red3") 
```

![](103_Distance_along_coast_files/figure-html/unnamed-chunk-16-1.png)<!-- -->


### Get Dist_along_coast   
For all stations  

```r
get_distance_along_coast_s <- safely(get_distance_along_coast)

rows <- 1:nrow(df_stations)

result_list <- rows %>%
  map(~ get_distance_along_coast_s(
    point = list(x = df_stations$x[.], y = df_stations$y[.]), 
    df_segments = coast, 
    df_segments_distances = coastsegment_distance)) %>%
  purrr::transpose()

# Sum up results
ok <- result_list$error %>% map_lgl(is.null)
dist <- result_list$result[ok] %>% map_dbl(~.$distance)

# Number that didnæt work
cat("Dist_along_coast found for", sum(ok), "stations \n")
cat("Dist_along_coast not found for", sum(!ok), "stations \n")
```

```
## Dist_along_coast found for 62 stations 
## Dist_along_coast not found for 0 stations
```

### Check out the ones with error, if necessary

```r
if (sum(!ok) > 0){

    i <- which(!ok)[1]
  
  # debugonce(get_distance_along_coast_all)
  check_distance_along_coast(
    point = list(x = df_stations$x[i], y = df_stations$y[i]), 
    df_segments = coast, 
    df_segments_distances = coastsegment_distance)
  
}
```

### Add Dist_along_coast  


```r
df_stations$Dist_along_coast <- NA
df_stations$Dist_along_coast[ok] <- dist
```

## 7. Save  

```r
saveRDS(df_stations, "Data/103_Selected_stations.rds")
```

