# setup ----
library(tidyverse)  # process dataframe
library(lubridate)  # process date time 
library(readr)      # read csv
library(sf)         # process GIS files
library(stars)      # process raster files
library(viridis)    # color scale
library(terra)
library(tmap)

# read data ----
dir_oras <- "./Env/ORAS5_0.25deg_1958-2022"
file_oras <- "oras5.tif"

# load file
oras <- read_stars(file.path(dir_oras, file_oras))

# plot
tm_shape(oras[,,,1]) + tm_raster()

# process data - ices 4bc ----
#### summarize by ices area
# load ices gis data
dir_gis <- "./Admin" 
file_gis <- "ices_areas_sub_4bc_7a_8ab_4326.gpkg" 

ices_area <- read_sf(file.path(dir_gis, file_gis))
ices_area <- ices_area %>% select(Area_27) #to have ICES Area only

# define area of interest
bb <- st_bbox(ices_area)

# visualise
tm_shape(oras[,,,1]) +
  tm_raster() +
  tm_shape(ices_area) +
  tm_borders()

#### extract mean temp
oras_mean <- aggregate(oras, ices_area, FUN = mean, na.rm = T)
oras_df <- left_join(as.data.frame(oras_mean), as.data.frame(ices_area), by = join_by(geom))
oras_df <- oras_df %>% 
  select(band, `oras5.tif`, Area_27) %>%
  mutate(IcesArea = Area_27,
         oras_sbt = `oras5.tif`,
         Date = as_date(band),
         Year = year(band)) %>% 
  select(IcesArea:Year) %>% 
  arrange(IcesArea, Year)

#### save file
write_rds(oras_df, "./Env/ORAS5_0.25deg_1958-2022/oras5_ices.rds")

# process data - datras ----
#### summarize by ices area
# load ices gis data
dir_gis <- "./ICES/DATRAS" 
file_gis <- "hl_loc_4abc7a8ab.gpkg" 

ices_area <- read_sf(file.path(dir_gis, file_gis))
ices_area <- ices_area %>% select(Area_27) #to have ICES Area only

# define area of interest
bb <- st_bbox(ices_area)

# visualise
tm_shape(oras[,,,1]) +
  tm_raster() +
  tm_shape(ices_area) +
  tm_borders()

#### extract mean temp
oras_mean <- aggregate(oras, ices_area, FUN = mean, na.rm = T)
oras_df <- left_join(as.data.frame(oras_mean), as.data.frame(ices_area), by = join_by(geom))
oras_df <- oras_df %>% 
  select(band, `oras5.tif`, Area_27) %>%
  mutate(IcesArea = Area_27,
         oras_sbt = `oras5.tif`,
         Date = as_date(band),
         Year = year(band)) %>% 
  select(IcesArea:Year) %>% 
  arrange(IcesArea, Year)

#### save file
write_rds(oras_df, "./Env/ORAS5_0.25deg_1958-2022/oras5_datras.rds")
