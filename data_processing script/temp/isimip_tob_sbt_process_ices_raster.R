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
library(terra)      # read/write raster file

# 2. READ DATA ---------------------------------------------------------------

# ICES AREA - 4abc
dir_gis <- "./Admin"
file_gis <- "ices_areas_sub_group_4abc_4326_new.gpkg"

ices_area <- vect(file.path(dir_gis, file_gis))

# define area of interest
bb <- ext(ices_area)

# ISIMIP3b MPI-ESM1-2HR (tob - Sea Water Potential Temperature at Sea Floor)

dir_tob <- "D:/OneDrive - UGent/data/Env/Sea Bottom Temperature_ISIMIP3b_MPI-ESM1-2HR_0.5deg_1850-2100"

file_tob_hist <- "mpi-esm1-2-hr_r1i1p1f1_historical_tob_30arcmin_global_monthly_1850_2014.nc"
file_tob_ssp126 <- "mpi-esm1-2-hr_r1i1p1f1_ssp126_tob_30arcmin_global_monthly_2015_2100.nc" 
file_tob_ssp585 <- "mpi-esm1-2-hr_r1i1p1f1_ssp585_tob_30arcmin_global_monthly_2015_2100.nc" 

# 3. PROCESS DATA --------------------------------------------------------------

# 3.1. Historical data ----------------------------------------------------

# 1. read file
nc_hist <- rast(file.path(dir_tob, file_tob_hist)) 

# 2. crop by area of interest (identified by comparing cropped results and ices area gis data)
nc_hist <- crop(nc_hist, bb)

# 3. change date time file - see file reading as stars
nc_hist_star <- read_ncdf(file.path(dir_tob, file_tob_hist))

year <- seq(1850,2014) 
month <- seq(1,12)  
year_all <- rep(year, each = 12)
month_all <- rep(month, length(year))
date <- as.Date(paste0("01", "-", month_all, "-", year_all),format = "%d-%m-%Y")  

time(nc_hist) <- date
names(nc_hist) <- date

# 3.2. ssp126 data ----------------------------------------------------

# 1. read file
nc_ssp126 <- rast(file.path(dir_tob, file_tob_ssp126)) 

# 2. crop by area of interest (identified by comparing cropped results and ices area gis data)
nc_ssp126 <- crop(nc_ssp126, bb)

# 3. change date time file - see file reading as stars
nc_ssp126_star <- read_ncdf(file.path(dir_tob, file_tob_ssp126))

year <- seq(2015,2100) 
month <- seq(1,12)  
year_all <- rep(year, each = 12)
month_all <- rep(month, length(year))
date <- as.Date(paste0("01", "-", month_all, "-", year_all),format = "%d-%m-%Y")  

time(nc_ssp126) <- date
names(nc_ssp126) <- date

# 3.2. ssp126 data ----------------------------------------------------

# 1. read file
nc_ssp585 <- rast(file.path(dir_tob, file_tob_ssp585)) 

# 2. crop by area of interest (identified by comparing cropped results and ices area gis data)
nc_ssp585 <- crop(nc_ssp585, bb)

# 3. change date time file - see file reading as stars
nc_ssp585_star <- read_ncdf(file.path(dir_tob, file_tob_ssp585))

year <- seq(2015, 2100) 
month <- seq(1, 12)  
year_all <- rep(year, each = 12)
month_all <- rep(month, length(year))
date <- as.Date(paste0("01", "-", month_all, "-", year_all),format = "%d-%m-%Y")  

time(nc_ssp585) <- date
names(nc_ssp585) <- date

# 3.3. Combine and save file ----
# combine file
nc_hist_ssp126 <- c(nc_hist, nc_ssp126)
nc_hist_ssp585 <- c(nc_hist, nc_ssp585)

# save file
writeRaster(nc_hist_ssp126, file.path(dir_tob, "isimip_sbt_hist_ssp126.tif"))
writeRaster(nc_hist_ssp585, file.path(dir_tob, "isimip_sbt_hist_ssp585.tif"))



