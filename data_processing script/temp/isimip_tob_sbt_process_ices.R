# **************************************
# MPI-ESM TEMPERATURE DATA PROCESSING **
# **************************************

# INFO --------------------------------------------------------------------

# PROJECT: FWO PHD - WP1
# Tuan Anh Bui (15/02/2023)

# Process the ISIMIP3b MPI-ESM1-2HR Sea Water Potential Temperature at Sea Floor (tob) data (0.5 deg)
# Input: netcdf files of ISIMIP3b MPI-ESM1-2HR: historical (1850-2014) and prediction (2015-2100) at 3 scenarios ssp126, 370, 585
# Output: 
# 1. dataframe of temperature
# 2. raster stack/netcdf file with ISIMIP3b MPI-ESM1-2HR sea bottom temperature over years

# 1. SETUP ----------------------------------------------------------------

library(sf)         # process GIS files
library(stars)      # process raster
library(tidyverse)  # process dataframe
library(lubridate)  # process date time 
library(readr)      # read/write rds file

# 2. READ DATA ---------------------------------------------------------------

# ICES AREA
dir_gis <- "D:/OneDrive - UGent/data/Admin"
file_gis <- "ices_areas_sub_group_4326_new.gpkg"

ices_area <- read_sf(file.path(dir_gis, file_gis))
ices_area <- ices_area %>% select(Area_27) #to have ICES Area only

# define area of interest
bb <- st_bbox(ices_area)

# ISIMIP3b MPI-ESM1-2HR (tob - Sea Water Potential Temperature at Sea Floor)

dir_tob <- "D:/OneDrive - UGent/data/Env/Sea Bottom Temperature_ISIMIP3b_MPI-ESM1-2HR_0.5deg_1850-2100"

file_tob_hist <- "mpi-esm1-2-hr_r1i1p1f1_historical_tob_30arcmin_global_monthly_1850_2014.nc"
file_tob_ssp126 <- "mpi-esm1-2-hr_r1i1p1f1_ssp126_tob_30arcmin_global_monthly_2015_2100.nc" 
file_tob_ssp370 <- "mpi-esm1-2-hr_r1i1p1f1_ssp370_tob_30arcmin_global_monthly_2015_2100.nc"
file_tob_ssp585 <- "mpi-esm1-2-hr_r1i1p1f1_ssp585_tob_30arcmin_global_monthly_2015_2100.nc"

# 3. PROCESS DATA --------------------------------------------------------------

# 3.1. Historical data ----------------------------------------------------

# 1. read file
nc_hist <- read_ncdf(file.path(dir_tob, file_tob_hist), proxy = TRUE) #proxy = TRUE to avoid crash when file is too large

#NOTE: to plot nc data, we can subset with 4 inputs: var, lon, lat, time
plot(nc_hist[,,,1]) 

# 2. crop by area of interest (identified by comparing cropped results and ices area gis data)
nc_hist = nc_hist[bb]

# 3. extract sbt (tob) by ices areas (aggregate function)
nc_mean_hist <- aggregate(nc_hist, ices_area, FUN = mean, na.rm = T)

# convert to data frame
nc_df_hist <- left_join(as.data.frame(nc_mean_hist), as.data.frame(ices_area), by = c("geometry" = "geom"))
nc_df_hist <- nc_df_hist %>% select(time, tob, Area_27) #remove geometry column
nc_df_hist <- nc_df_hist %>% 
  mutate(IcesArea = Area_27,
         isimip_sbt = tob,
         Date = as_date(time),
         Year = year(time)) %>% 
  select(IcesArea:Year) %>% 
  arrange(IcesArea, Year)

# 3.2. ssp126 data ----------------------------------------------------

# 1. read file
nc_ssp126 <- read_ncdf(file.path(dir_tob, file_tob_ssp126), proxy = TRUE) #proxy = TRUE to avoid crash when file is too large

#NOTE: to plot nc data, we can subset with 4 inputs: var, lon, lat, time
plot(nc_ssp126[,,,1]) 

# 2. crop by area of interest (identified by comparing cropped results and ices area gis data)
nc_ssp126 = nc_ssp126[bb]

# 3. extract sbt (tob) by ices areas (aggregate function)
nc_mean_ssp126 <- aggregate(nc_ssp126, ices_area, FUN = mean, na.rm = T)

# convert to data frame
nc_df_ssp126 <- left_join(as.data.frame(nc_mean_ssp126), as.data.frame(ices_area), by = c("geometry" = "geom"))
nc_df_ssp126 <- nc_df_ssp126 %>% select(time, tob, Area_27) #remove geometry column
nc_df_ssp126 <- nc_df_ssp126 %>% 
  mutate(IcesArea = Area_27,
         isimip_sbt = tob,
         Date = as_date(time),
         Year = year(time)) %>% 
  select(IcesArea:Year) %>% 
  arrange(IcesArea, Year)

# 3.3. ssp370 data ----------------------------------------------------

# 1. read file
nc_ssp370 <- read_ncdf(file.path(dir_tob, file_tob_ssp370), proxy = TRUE) #proxy = TRUE to avoid crash when file is too large

#NOTE: to plot nc data, we can subset with 4 inputs: var, lon, lat, time
plot(nc_ssp370[,,,1]) 

# 2. crop by area of interest (identified by comparing cropped results and ices area gis data)
nc_ssp370 = nc_ssp370[bb]

# 3. extract sbt (tob) by ices areas (aggregate function)
nc_mean_ssp370 <- aggregate(nc_ssp370, ices_area, FUN = mean, na.rm = T)

# convert to data frame
nc_df_ssp370 <- left_join(as.data.frame(nc_mean_ssp370), as.data.frame(ices_area), by = c("geometry" = "geom"))
nc_df_ssp370 <- nc_df_ssp370 %>% select(time, tob, Area_27) #remove geometry column
nc_df_ssp370 <- nc_df_ssp370 %>% 
  mutate(IcesArea = Area_27,
         isimip_sbt = tob,
         Date = as_date(time),
         Year = year(time)) %>% 
  select(IcesArea:Year) %>% 
  arrange(IcesArea, Year)

# 3.2. ssp585 data ----------------------------------------------------

# 1. read file
nc_ssp585 <- read_ncdf(file.path(dir_tob, file_tob_ssp585), proxy = TRUE) #proxy = TRUE to avoid crash when file is too large

#NOTE: to plot nc data, we can subset with 4 inputs: var, lon, lat, time
plot(nc_ssp585[,,,1]) 

# 2. crop by area of interest (identified by comparing cropped results and ices area gis data)
nc_ssp585 = nc_ssp585[bb]

# 3. extract sbt (tob) by ices areas (aggregate function)
nc_mean_ssp585 <- aggregate(nc_ssp585, ices_area, FUN = mean, na.rm = T)

# convert to data frame
nc_df_ssp585 <- left_join(as.data.frame(nc_mean_ssp585), as.data.frame(ices_area), by = c("geometry" = "geom"))
nc_df_ssp585 <- nc_df_ssp585 %>% select(time, tob, Area_27) #remove geometry column
nc_df_ssp585 <- nc_df_ssp585 %>% 
  mutate(IcesArea = Area_27,
         isimip_sbt = tob,
         Date = as_date(time),
         Year = year(time)) %>% 
  select(IcesArea:Year) %>% 
  arrange(IcesArea, Year)

# 4. SAVE DATA ---------------------------------------------------------------
df_hist_ssp126 <- bind_rows(nc_df_hist, nc_df_ssp126)
df_hist_ssp370 <- bind_rows(nc_df_hist, nc_df_ssp370)
df_hist_ssp585 <- bind_rows(nc_df_hist, nc_df_ssp585)

write_rds(df_hist_ssp126, file.path(dir_tob, "isimip_sbt_ices_hist_ssp126.rds"))
write_rds(df_hist_ssp370, file.path(dir_tob, "isimip_sbt_ices_hist_ssp370.rds"))
write_rds(df_hist_ssp585, file.path(dir_tob, "isimip_sbt_ices_hist_ssp585.rds"))

