# setup ----
library(tidyverse)
library(lubridate)  # process date time 
library(readr)      # read csv
library(sf)         # process GIS files
library(stars)      # process raster files
library(viridis)    # color scale
library(terra)
library(tmap)

# ices data ----
# load ices gis data
dir_icesdiv <- "D:/OneDrive - UGent/data/Admin" 
file_icesdiv <- "ices_areas_sub_group_4326_new.gpkg" 
path_icesdiv <- file.path(dir_icesdiv, file_icesdiv)

icesdiv <- read_sf(path_icesdiv)
icesdiv_area <- icesdiv$Area_27

ggplot(data = icesdiv) +
  geom_sf()

#### get bounding box 
bbox <- st_bbox(icesdiv)
#xmin      ymin      xmax      ymax 
#-8.746612 43.314507 12.005942 57.500000 

# temp data ----
#### load data
dir_data <- "./NEMO MEDUSA/Sea floor"
nm <- read_rds(file.path(dir_data, "NM.01.1997.rds"))

#### filter by icesdiv bbox
nm_ices <- nm %>% 
  filter(longitude >= bbox[1], longitude <= bbox[3],
         latitude >= bbox[2], latitude <= bbox[4])

ggplot() +
  geom_point(data = nm_ices, aes(x = longitude, y = latitude), 
             color = "grey", alpha = 0.5) +
  geom_sf(data = icesdiv, alpha = 0, color = "black") +
  coord_sf() 

#### convert to sf
# assuming crs WGS84 (EPSG 4326)
projcrs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
nm_ices_sf <- st_as_sf(nm_ices, 
                       coords = c("longitude", "latitude"),
                       crs = projcrs)
ggplot() +
  geom_sf(data = nm_ices_sf, color = "grey", alpha = 0.5) +
  geom_sf(data = icesdiv, alpha = 0, color = "black") +
  coord_sf()

#### intersect nm_ices_sf with icesdiv to keep only points within icesdiv
nm_aoi <- st_intersection(nm_ices_sf, icesdiv)

ggplot() +
  geom_sf(data = nm_aoi, color = "grey", alpha = 0.5) +
  geom_sf(data = icesdiv, alpha = 0, color = "black") +
  coord_sf()

#### turn sf to df and summarize
nm_aoi_df <- as.data.frame(nm_aoi)
nm_aoi_df %>% 
  group_by(Area_27) %>%
  summarize(Temperature = mean(Temperature),
            DIN = mean(DIN),
            Detritus = mean(Detritus),
            Phytoplankton = mean(Phytoplankton),
            Year = unique(Year),
            Month = unique(Month))

# combine temp data -----
dir_data <- "./NEMO MEDUSA/Sea floor"
list_file <- list.files(dir_data)

# create an empty dataframe to save results
nm_sub <- tibble()

# for loop to combine NEMO MEDUSA data
for(i in 1:length(list_file)) { 
  #### note
  print(paste("processing", list_file[i]))
  
  #### load data month by month
  nm <- read_rds(file.path(dir_data, list_file[i]))

  #### filter by icesdiv bbox
  nm_ices <- nm %>% 
    filter(longitude >= bbox[1], longitude <= bbox[3],
           latitude >= bbox[2], latitude <= bbox[4])
  
  #### convert to sf 
  # assuming crs WGS84 (EPSG 4326)
  projcrs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  nm_ices_sf <- st_as_sf(nm_ices, 
                         coords = c("longitude", "latitude"),
                         crs = projcrs)
  
  #### intersect nm_ices_sf with icesdiv to keep only points within icesdiv
  nm_aoi <- st_intersection(nm_ices_sf, icesdiv)
  
  #### turn sf to df and summarize
  nm_aoi_df <- as.data.frame(nm_aoi) %>%
    mutate(Population = Area_27)
  
  nm_aoi_df_sum <- nm_aoi_df %>% 
    group_by(Population) %>%
    summarize(Temperature = mean(Temperature),
              DIN = mean(DIN),
              Detritus = mean(Detritus),
              Phytoplankton = mean(Phytoplankton),
              Year = unique(Year),
              Month = unique(Month))
  
  #### bind all data together
  nm_sub <- bind_rows(nm_sub, nm_aoi_df_sum)
}

#### save file
write_rds(nm_sub, file.path(dir_data, "nemomedusa_ices.rds"))

#### plot
nm_sub_sum <- nm_sub %>% 
  group_by(Population, Year) %>% 
  summarize(Temperature = mean(Temperature)) %>%
  filter(Population %in% c("4bc", "7a", "8ab"),
         Year <= 2021)

ggplot(data = nm_sub_sum, aes(x = Year, y = Temperature, color = Population)) +
  geom_line() 


