# ********************************************
# ICES HIGH RESOLUTION CTD DATA EXPLORATION **
# ********************************************

# INFO --------------------------------------------------------------------

# PROJECT: FWO PHD - WP1
# Tuan Anh Bui (20/02/2023)

# Explore and process ICES High resolution CTD data
# Processing procedure:
# 1. Process CTD data: for each unique CTD event, get the maximum depth if there are multiple depth values and 
# get the lowest temperature values (assuming that lowest temperature value is at near-bottom depth as CTD gets deeper during measurement)
# 2. Validate if the CTD measurement is actually at near-bottom,
# by comparing depth from CTD data to depth from Emodnet digital bathymetry (100m resolution)

# 1. SETUP ----------------------------------------------------------------

library(readxl)     # read excel
library(tidyverse)  # process dataframe
library(lubridate)  # process date time 
library(readr)      # read csv
library(writexl)    # write csv
library(sf)         # process GIS files
library(tmap)       # plot map

# 2. READ AND PROCESS DATA ---------------------------------------------------------------

# ices area data
dir_gis <- "./Admin"
file_gis <- "ices_areas_sub_group_4abc_4326_new.gpkg" #"ices_areas_sub_group_4326_new.gpkg"

ices_area <- read_sf(file.path(dir_gis, file_gis))
ices_area <- ices_area %>% select(Area_27) #to have ICES Area only

# ctd data
dir_ctd <- "./Env/Sea Bottom Temperature_ICES CTD_1970-2020"
file_ctd <- "01171629.txt"

ctd <- read_delim(file = file.path(dir_ctd, file_ctd), 
                  delim = "\t", escape_double = FALSE, 
                  trim_ws = TRUE)

ctd <- ctd %>% select(Cruise, Station, Type, Year, Month, Day, Hour, Minute, 
                      `Latitude [degrees_north]`, `Longitude [degrees_east]`,
                      `Bot. Depth [m]`, `Temperature [degC]`)


# Summary 
summary(ctd %>% select(Year, `Latitude [degrees_north]`, `Longitude [degrees_east]`,
                       `Bot. Depth [m]`, `Temperature [degC]`))

## Process ctd ----
ctd <- ctd %>% filter(is.na(`Temperature [degC]`) == F) # Remove NA Temp
ctd <- ctd %>% mutate(sample_date = as_date(paste(Year, Month, Day, sep = "-")),
                      lat_deg = `Latitude [degrees_north]`,
                      lon_deg = `Longitude [degrees_east]`,
                      depth_m =  `Bot. Depth [m]`,
                      temp_degC = `Temperature [degC]`)


ctd_sub <- ctd %>% select(Cruise, Station, Type, Year, Month, Day, 
                          sample_date, lat_deg, lon_deg, depth_m, temp_degC)

# Normally, the CTD temperature decreases when the depth increase (from surface to bottom). 
# However, most or almost all records, there is only 1 record of depth (bottom depth??) 
# so we assume the lowest temperature is the temperature at bottom.


# To get temperature at depth we need to do the following steps:

# 0. Create ID: 1 measurement position has 1 compound of Cruise, Station, sample_date, lat_deg, lon_deg
# 1. See which positions have 1 depth, which have multiple depth
# 2. Summarize with depth, split into those with 1 depth and with multiple depth
# Final: ctd data with depth

# 0. create ID compound: Cruise - Station - sample_date - lat_deg - lon_deg
ctd_sub <- ctd_sub %>% mutate(ID = paste(Cruise, Station, sample_date, lat_deg, lon_deg))

# 1. see which positions have 1 depth, which has multiple depth
ctd_depth <- ctd_sub %>% 
  group_by(ID, depth_m) %>% summarize(n = n()) %>% 
  group_by(ID) %>% summarize(n = n())

# 2. summarize with depth
# 1 depth
ctd_depth1 <- ctd_depth %>% filter(n == 1)

ctd_sub_depth1 <- ctd_sub %>% filter(ID %in% ctd_depth1$ID) 

ctd_sbt_depth1 <- ctd_sub_depth1 %>% 
  group_by(ID, Cruise, Station, Year, Month, Day, 
           sample_date, lat_deg, lon_deg, depth_m) %>% 
  summarize(temp_degC = min(temp_degC))

# multiple depth
ctd_depth_multi <- ctd_depth %>% filter(n != 1)

ctd_sub_depth_multi <- ctd_sub %>% filter(ID %in% ctd_depth_multi$ID) 

# add ID depth to split the maximum depth
ctd_sub_depth_multi <- ctd_sub_depth_multi %>% mutate(ID_depth = paste(ID, depth_m))

# find max depth
ctd_sub_depth_multi_depth <- ctd_sub_depth_multi %>% group_by(ID) %>% summarize(depth_m = max(depth_m))
ctd_sub_depth_multi_depth <- ctd_sub_depth_multi_depth %>% mutate(ID_depth = paste(ID, depth_m))

ctd_sub_depth_multi_sub <- ctd_sub_depth_multi %>% filter(ID_depth %in% ctd_sub_depth_multi_depth$ID_depth)

# summarize
ctd_sbt_depth_multi <- ctd_sub_depth_multi_sub %>% 
  group_by(ID, Cruise, Station, Year, Month, Day, 
           sample_date, lat_deg, lon_deg, depth_m) %>% 
  summarize(temp_degC = min(temp_degC))

# merge depth 1 and depth multi
ctd_sbt <- bind_rows(ctd_sbt_depth1, ctd_sbt_depth_multi)

# transform to sf format
ctd_sbt_sf <- st_as_sf(ctd_sbt, 
                       coords = c("lon_deg","lat_deg"), # c(lon,lat) 
                       crs = 4326) 

# filter only ctd points within ices_area
ctd_sbt_ices <- st_intersection(ctd_sbt_sf, ices_area)

## compare depth from ctd data and bathymetry data ----
# bathymetry data from Emodnet digital bathymetry (DTM 2020)  
dir_bathy <- "./Env/Bathymetry_Emodnet_100m_2020"

emodnet <- read_stars("./Env/Bathymetry_Emodnet_100m_2020/Bathymetry_Emodnet_100m_2020.tif")
tm_shape(emodnet) +
  tm_raster(style = "quantile")  

# extract bathymetry from emodnet at ctd locations
ctd_emodnet_depth <- st_extract(emodnet, ctd_sbt_ices)

ctd_emodnet <- left_join(as.data.frame(ctd_sbt_ices), 
                      as.data.frame(ctd_emodnet_depth),
                      by = join_by(geom)) %>%
  unique() %>%
  st_as_sf(crs = 4326)

# rename and convert bathymetry to positive values
ctd_emodnet <- ctd_emodnet %>% 
  rename(bathy_emodnet = Bathymetry_Emodnet_100m_2020.tif) %>%
  mutate(bathy_emodnet = 0 - bathy_emodnet)

# check negative values (points on land) -> convert to NA
tm_shape(ctd_emodnet %>% filter(bathy_emodnet < 0)) +
  tm_dots() +
  tm_shape(emodnet) +
  tm_raster(style = "quantile")

tm_shape(ctd_emodnet %>% filter(depth_m <= 0)) +
  tm_dots() +
  tm_shape(emodnet) +
  tm_raster(style = "quantile")

# check extreme values (> 5000m) due to NA in extraction -> convert to NA
max(ctd_emodnet$depth_m, na.rm = T) #maximum depth from ctd is > 4000
# no points with depth within 5-10000m -> extreme values are those >= 5000

tm_shape(ctd_emodnet %>% filter(bathy_emodnet >= 5000)) +
  tm_dots() +
  tm_shape(emodnet) +
  tm_raster(style = "quantile")

# convert negative and extreme values (> 5000m) to NA
ctd_emodnet <- ctd_emodnet %>%
  mutate(bathy_emodnet = if_else(bathy_emodnet < 0, NA, 
                                 if_else(bathy_emodnet >= 5000, NA, bathy_emodnet)))

#### save file
st_write(ctd_emodnet, file.path(dir_ctd, "ices_ctd_sea-bottom-temperature.gpkg"), delete_dsn = T)

#### correlation depth ctd and emodnet
ggplot(data = ctd_emodnet, aes(x = bathy_emodnet, y = depth_m)) +
  geom_point() +
  geom_smooth(method = 'lm')

lm <- lm(depth_m ~ bathy_emodnet, data = ctd_emodnet)
nrow(ctd_emodnet %>% filter(is.na(bathy_emodnet) == F))/nrow(ctd_emodnet) #99.5% available data
round(summary(lm)$adj.r.squared, 2) #0.98

