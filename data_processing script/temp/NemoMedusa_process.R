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
dir_nm <- "./Env/NEMO MEDUSA"

# load file
nm <- read_stars(file.path(dir_nm, "nemomedusa.tif"))
nm <- nm[,,,2:517]
#note: 1st layer is bathymetry, subsequent layers are monthly bottom temperature

# plot
tm_shape(nm[,,,1]) + tm_raster()

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
tm_shape(nm[,,,1]) +
  tm_raster() +
  tm_shape(ices_area) +
  tm_borders()

#### extract mean temp
nm_mean <- aggregate(nm, ices_area, FUN = mean, na.rm = T)
nm_df <- left_join(as.data.frame(nm_mean), as.data.frame(ices_area), by = join_by(geom))
nm_df <- nm_df %>% 
  select(band, `nemomedusa.tif`, Area_27) %>%
  mutate(IcesArea = Area_27,
         nemomedusa_sbt = `nemomedusa.tif`,
         Date = as.Date(paste0("01", "-", str_sub(band, start = 11, end = 12), "-", str_sub(band, start = 7, end = 10)),format = "%d-%m-%Y"),
         Year = str_sub(band, start = 7, end = 10)) %>% 
  select(IcesArea:Year) %>% 
  arrange(IcesArea, Year)

#### save file
write_rds(nm_df, file.path(dir_nm, "nemomedusa_ices.rds"))

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
tm_shape(nm[,,,1]) +
  tm_raster() +
  tm_shape(ices_area) +
  tm_borders()

#### extract mean temp
nm_mean <- aggregate(nm, ices_area, FUN = mean, na.rm = T)
nm_df <- left_join(as.data.frame(nm_mean), as.data.frame(ices_area), by = join_by(geom))
nm_df <- nm_df %>% 
  select(band, `nemomedusa.tif`, Area_27) %>%
  mutate(IcesArea = Area_27,
         nemomedusa_sbt = `nemomedusa.tif`,
         Date = as.Date(paste0("01", "-", str_sub(band, start = 11, end = 12), "-", str_sub(band, start = 7, end = 10)),format = "%d-%m-%Y"),
         Year = str_sub(band, start = 7, end = 10)) %>% 
  select(IcesArea:Year) %>% 
  arrange(IcesArea, Year)

#### save file
write_rds(nm_df, file.path(dir_nm, "nemomedusa_datras.rds"))



