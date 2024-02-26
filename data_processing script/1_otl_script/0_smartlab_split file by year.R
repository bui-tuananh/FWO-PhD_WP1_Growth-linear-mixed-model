# ***************************
# OTOLITH DATA EXPLORATION **
# ***************************

# INFO --------------------------------------------------------------------

# PROJECT: FWO PHD - WP1
# Tuan Anh Bui (19/06/2022)

# Split selected otoliths into files to import into SmartLab (using smartlab template)

# 1. SETUP ----------------------------------------------------------------

library(readxl)     # read excel
library(tidyverse)  # process dataframe
library(lubridate)  # process date time 
library(readr)      # read csv
library(writexl)    # write csv


# 2. READ AND PROCESS DATA ---------------------------------------------------------------

dir_otolith <- "D:/OneDrive - UGent/data/WP1/otolith/@processed"
otolith <- read_rds(file.path(dir_otolith, "sol_select_full.rds"))

# 2.1. 7A FIRST ATTEMPT ------------------------------------------
# First attempt sampling 7A with a subset of otoliths

# sol selected 7a
dir_otolith_old <- "D:/OneDrive - UGent/data/WP1/otolith/ilvo_old"

otolith_old <- read_rds(file.path(dir_otolith_old, "sol_select_new.rds"))
otolith_7a_old <- otolith_old %>% filter(HAU.IcesArea == "7a")

# template smartlab
dir_smartlab <- "D:/OneDrive - UGent/data/WP1/otolith/SmartLab"
file_smartlab <- "Sample Sjabloon.xlsx"
path_smartlab <- file.path(dir_smartlab, file_smartlab)

df_smartlab <- read_excel(path_smartlab)
names(df_smartlab)

# CONVERT TO SMARTLAB TEMPLATE --------------------------------------------

# Conver SPE.Length (before 2000 in cm, after 2000 in mm)
otolith_7a_old$SPE.Length <- if_else(otolith_7a_old$SPE.Length < 100, otolith_7a_old$SPE.Length*10, otolith_7a_old$SPE.Length)

# Convert HAU.IcesArea to format of SmartLab
otolith_7a_old$HAU.IcesArea <- "7A"

# convert otolith fields to be the same as df_smartlab
otolith_7a_old <- otolith_7a_old %>% mutate(Number = Smartlab.Number,
                                    `Number External` = UniqueID,
                                    Species = SAM.SpeciesFaoCode,
                                    `Ices Division` = HAU.IcesArea,
                                    `Statistical Rectangle` = NA,
                                    `Catch Date (dd/mm/jjjj)` = format(as_date(SAM.Date), "%d/%m/%Y"), #dd/mm/yyyy
                                    `Length (in mm)` = SPE.Length,
                                    `Weight (in gram)` = SPE.Weight,
                                    `Maturity` = NA,
                                    `Sex` = "F",
                                    `Position Receipt` = NA,
                                    `Description` = NA,
                                    `Comment` = NA)

otolith_7a_sub <- otolith_7a_old %>% select(TRI.Year, Number:`Weight (in gram)`,Maturity, Sex, `Position Receipt`, `Description`, `Comment`)

# split data by year
list_year <- sort(unique(otolith_7a_sub$TRI.Year))

dir_name <- "D:/OneDrive - UGent/data/WP1/otolith/SmartLab/tbui_7a_split"

for (i in 1:length(list_year)) {
  
  print(paste0("Processing Year ", list_year[i]))
  df <- otolith_7a_sub %>% filter(TRI.Year == list_year[i])
  df <- df %>% select(-TRI.Year)
  
  
  # save file
  file_name <- paste0("tbui_7a_SOL_",list_year[i],".xlsx")
  path_name <- file.path(dir_name, file_name)
  
  write_xlsx(list(`_SampleInvoer` = df), path = path_name) # save with sheet name: _SampleInvoer (SmartLab sample)
}


# 2.2. 7A ADDITION -------------------------------------------------------------------------

# Metadata of otoliths that were not included in the first sampling attempt

# WORKFLOW
# 1. Subset sol_select_full without otoliths of the 1st attemp (sol_select_new)
# 2. Split data by Year


# sol selected 7a

# only select additional 
otolith_add7a <- otolith %>% filter(!(UniqueID %in% otolith_7a_old$UniqueID), HAU.IcesAreaGroup == "7a")


# CONVERT TO SMARTLAB TEMPLATE --------------------------------------------

# template smartlab
dir_smartlab <- "D:/OneDrive - UGent/data/WP1/otolith/SmartLab"
file_smartlab <- "Sample Sjabloon.xlsx"
path_smartlab <- file.path(dir_smartlab, file_smartlab)

df_smartlab <- read_excel(path_smartlab)
names(df_smartlab)

# Convert HAU.IcesArea to format of SmartLab
otolith_add7a$HAU.IcesArea <- "7A"

# convert otolith fields to be the same as df_smartlab
otolith_add7a <- otolith_add7a %>% mutate(Number = Smartlab.Number,
                                            `Number External` = UniqueID,
                                            Species = SAM.SpeciesFaoCode,
                                            `Ices Division` = HAU.IcesAreaGroup,
                                            `Statistical Rectangle` = NA,
                                            `Catch Date (dd/mm/jjjj)` = format(as_date(SAM.Date), "%d/%m/%Y"), #dd/mm/yyyy
                                            `Length (in mm)` = SPE.Length,
                                            `Weight (in gram)` = SPE.Weight,
                                            `Maturity` = NA,
                                            `Sex` = "F",
                                            `Position Receipt` = NA,
                                            `Description` = NA,
                                            `Comment` = NA)


otolith_7a_sub <- otolith_add7a %>% select(TRI.Year, Number:`Weight (in gram)`,Maturity, Sex, `Position Receipt`, `Description`, `Comment`)

# split data by year
list_year <- sort(unique(otolith_7a_sub$TRI.Year))

dir_name <- "D:/OneDrive - UGent/data/WP1/otolith/SmartLab/tbui_7a_split_addition"

for (i in 1:length(list_year)) {
  
  print(paste0("Processing Year ", list_year[i]))
  df <- otolith_7a_sub %>% filter(TRI.Year == list_year[i])
  df <- df %>% select(-TRI.Year)
  
  
  # save file
  file_name <- paste0("tbui_7a_SOL_",list_year[i],".xlsx")
  path_name <- file.path(dir_name, file_name)
  
  write_xlsx(list(`_SampleInvoer` = df), path = path_name) # save with sheet name: _SampleInvoer (SmartLab sample)
}



# 2.3. 4BC ------------------------------------------

# sol selected 4bc
otolith_4bc <- otolith %>% filter(HAU.IcesAreaGroup == "4bc")

# template smartlab
dir_smartlab <- "D:/OneDrive - UGent/data/WP1/otolith/SmartLab"
file_smartlab <- "Sample Sjabloon.xlsx"
path_smartlab <- file.path(dir_smartlab, file_smartlab)

df_smartlab <- read_excel(path_smartlab)
names(df_smartlab)

# CONVERT TO SMARTLAB TEMPLATE --------------------------------------------

# Convert HAU.IcesArea to format of SmartLab
otolith_4bc$HAU.IcesArea <- toupper(otolith_4bc$HAU.IcesArea) 

# convert otolith fields to be the same as df_smartlab
otolith_4bc <- otolith_4bc %>% mutate(Number = Smartlab.Number,
                                            `Number External` = UniqueID,
                                            Species = SAM.SpeciesFaoCode,
                                            `Ices Division` = HAU.IcesArea,
                                            `Statistical Rectangle` = NA,
                                            `Catch Date (dd/mm/jjjj)` = format(as_date(SAM.Date), "%d/%m/%Y"), #dd/mm/yyyy
                                            `Length (in mm)` = SPE.Length,
                                            `Weight (in gram)` = SPE.Weight,
                                            `Maturity` = NA,
                                            `Sex` = "F",
                                            `Position Receipt` = NA,
                                            `Description` = NA,
                                            `Comment` = NA)

otolith_4bc_sub <- otolith_4bc %>% select(TRI.Year, Number:`Weight (in gram)`,Maturity, Sex, `Position Receipt`, `Description`, `Comment`)

# split data by year
list_year <- sort(unique(otolith_4bc$TRI.Year))

dir_name <- "D:/OneDrive - UGent/data/WP1/otolith/SmartLab/tbui_4bc_split"

for (i in 1:length(list_year)) {
  
  print(paste0("Processing Year ", list_year[i]))
  df <- otolith_4bc_sub %>% filter(TRI.Year == list_year[i])
  df <- df %>% select(-TRI.Year)
  
  
  # save file
  file_name <- paste0("tbui_4bc_SOL_",list_year[i],".xlsx")
  path_name <- file.path(dir_name, file_name)
  
  write_xlsx(list(`_SampleInvoer` = df), path = path_name) # save with sheet name: _SampleInvoer (SmartLab sample)
}

# 2.4. 8AB ------------------------------------------

# sol selected 8ab
otolith_8ab <- otolith %>% filter(HAU.IcesAreaGroup == "8ab")

# template smartlab
dir_smartlab <- "D:/OneDrive - UGent/data/WP1/otolith/SmartLab"
file_smartlab <- "Sample Sjabloon.xlsx"
path_smartlab <- file.path(dir_smartlab, file_smartlab)

df_smartlab <- read_excel(path_smartlab)
names(df_smartlab)

# CONVERT TO SMARTLAB TEMPLATE --------------------------------------------

# Convert HAU.IcesArea to format of SmartLab
otolith_8ab$HAU.IcesArea <- "8A" # the Smartlab template allow either 8A/8B so select 8A for the sake of sampling 

# convert otolith fields to be the same as df_smartlab
otolith_8ab <- otolith_8ab %>% mutate(Number = Smartlab.Number,
                                      `Number External` = UniqueID,
                                      Species = SAM.SpeciesFaoCode,
                                      `Ices Division` = HAU.IcesArea,
                                      `Statistical Rectangle` = NA,
                                      `Catch Date (dd/mm/jjjj)` = format(as_date(SAM.Date), "%d/%m/%Y"), #dd/mm/yyyy
                                      `Length (in mm)` = SPE.Length,
                                      `Weight (in gram)` = SPE.Weight,
                                      `Maturity` = NA,
                                      `Sex` = "F",
                                      `Position Receipt` = NA,
                                      `Description` = NA,
                                      `Comment` = NA)

otolith_8ab_sub <- otolith_8ab %>% select(TRI.Year, Number:`Weight (in gram)`,Maturity, Sex, `Position Receipt`, `Description`, `Comment`)

# split data by year
list_year <- sort(unique(otolith_8ab$TRI.Year))

dir_name <- "D:/OneDrive - UGent/data/WP1/otolith/SmartLab/tbui_8ab_split"

for (i in 1:length(list_year)) {
  
  print(paste0("Processing Year ", list_year[i]))
  df <- otolith_8ab_sub %>% filter(TRI.Year == list_year[i])
  df <- df %>% select(-TRI.Year)
  
  
  # save file
  file_name <- paste0("tbui_8ab_SOL_",list_year[i],".xlsx")
  path_name <- file.path(dir_name, file_name)
  
  write_xlsx(list(`_SampleInvoer` = df), path = path_name) # save with sheet name: _SampleInvoer (SmartLab sample)
}


