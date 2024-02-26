

# INFO --------------------------------------------------------------------

# PROJECT: FWO PHD - WP1
# Tuan Anh Bui (02/12/2022)

# Rename ifremer data to be readable by Smartdot 


# 1. SETUP -------------------------------------------------------------------

library(tidyverse)
library(readr)


# 2. LOAD DATA -------------------------------------------------------------

dir <- "D:/OneDrive - UGent/data/WP1/otolith/@processed"
otl_full <- read_rds(file.path(dir, "sol_select_full.rds"))

# filter ifremer data (no TRI.Code)
otl_ifremer <- otl_full %>% filter(HAU.IcesAreaGroup == "8ab" & TRI.Year >= 2008 & is.na(TRI.Code) == T )


# 3. PROCESS DATA ---------------------------------------------------------

list_year <- sort(unique(otl_ifremer$TRI.Year))

for (i in 1:length(list_year)) {
  
  # select year
  print(paste("processing year", list_year[i]))
  otl_ifremer_year <- otl_ifremer %>% filter(TRI.Year == list_year[i])
  
  # setup path_old and path_new
  dir_old <- file.path("D:/OneDrive - UGent/data/WP1/otolith/ifremer",list_year[i])
  name_old <- paste0(otl_ifremer_year$UniqueID,".tif")
  path_old <- file.path(dir_old, name_old)
  
  dir_new <- file.path("D:/OneDrive - UGent/data/WP1/otolith/ifremer", paste0("new_",list_year[i]))
  name_new <- paste0("ifremer_", list_year[i], "_", otl_ifremer_year$Smartlab.Number,".tif")
  path_new <- file.path(dir_new, name_new)
  
  #copy file
  file.copy(path_old, path_new)
  
}

table(otl_ifremer$TRI.Year)
otl_ifremer %>% group_by(TRI.Year)


