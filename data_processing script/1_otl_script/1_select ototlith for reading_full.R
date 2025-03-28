# ***************************
# OTOLITH DATA EXPLORATION **
# ***************************

# INFO --------------------------------------------------------------------

# PROJECT: FWO PHD - WP1
# Tuan Anh Bui (10/10/2022)

# Prepare otolith metadata to picture and read in SmartDots

# INPUT: 2 raw files
# 1. metadata from 1969 to 2004 - test1 - extracted from excel files in ILVO server (further cleaning needed) 
# 2. metadata from 2004 to 2021 - sol_2004_2021.RDS - extracted from SmartFish


# 1969 to 2004 - otoliths can be found using the UniqueID - example: SOL_48_LPB2_Z.19_17-11-2003_2112

# SOL 		    |	species code
# 48   		    |	otolith sequence in the year collection in the archive
# LPB2		    |	fishing location (i.e. Liverpool Bay - Irish Sea)
# Z.19: 		  |	fishing vessel
# 17-11-2003	|	sampling  (fishing/otolith processing) date  
# 2112		    |	NA


# 2004 to 2021 - the following infomartion is needed

# SAM.Date				          |	sample date 
# TRI.Vessel                | vessel 
# SAM.FateCategory	      	|	sample category (discards or landings)
# SPE.Sequence			        |	otolith sequence
# Tri.Code		             	|	trip code - sometimes helpful


# a first attemp to select otolith was made in the otolith_selection for reading.R (3 otolith/age/year class) - file: sol_select_new.csv
# but in this a full set of otolith will be selected so more otolith can be read


# OUTPUT: list of all otolith with metadata to be found in the physical archive

# WORKFLOW: 
# 1. Process and merge pre- and post-2004 full data
# 2. Add processed data from sol_select_new.csv 

# 1. SETUP ----------------------------------------------------------------

library(readxl)     # read excel
library(tidyverse)  # process dataframe
library(lubridate)  # process date time 
library(readr)      # read csv
library(writexl)    # write csv

# 2. READ DATA ---------------------------------------------------------------

## sol pre 2004 ----
# Historical otolith data (1969-2004)
dir_otolith_raw <- "./data/otolith/@raw/ilvo"
ilvo <- read_xlsx(path = file.path(dir_otolith_raw, "test1.xlsx"))

# Filter only sol
sol <- ilvo %>% filter(FAA_species_code == "SOL") 

# Data description

# lengte          - length (cm)
# gewicht         - weight (g)
# jaarklas        - Year class (sometimes noted as age, sometimes noted as Year class)

# sex             -  "/"  "0"  "1"  "2"  "22" "3"  "32" (assume 1 and 2 are female?)
# gewicht_gonaden - gonad weight (g) - also interested, is the gonad weight reliable? (ranging from -2.2 to > 100000). 

# Confirmed
# Zone            - 100-200 - 4c 
#                 - 300-400 - 4b 
#                 - 500-700 - 4a
#                 - 800-1000 - unknown # what are these?
#                 - 1600    - 8ab
#                 - 2000    - 7a

# jaarklas W      - whole otolith method
# jaarklas S      - sectioned otolith method - should use this one


## sol post 2004 ----
sol_post04_raw <- read_rds(file = file.path(dir_otolith_raw, "sol_2004_2021.RDS"))

# 3. PROCESS DATA  ---------------------------------------------------------

# data selection criteria
# age >= 3
# yearclass >= 1957
# IcesArea 7a, 4a, 4b, 4c, 8ab
# sex female

# 3.1. SOL PRE 2004 -----------------------------------------------------------

# 3.1.1. clean data ---------------------------------------------------

# the data uses different codes compared to the current data so some cleaning needed

# Change scientific name 
sol$scientific_name <- "Solea solea"

# Convert zone to ices area (still need to ask zone 800-1000)
sort(unique(sol$Zone))

# Convert gonad weight (divided by 100) # Note 
#need double check because some recent years (2001 does not need division by 100)
sol$gewicht_gonaden <- sol$gewicht_gonaden/100

#problem
#2112-1044 (multiple regions in a trip? -> should not choose)

#View(sol %>% filter(Zone == "2112-1044")) #need to check if this zone can be 7a
sol$zone_certainty <- if_else(sol$Zone == "2112-1044", "uncertain", "certain")
sort(unique(sol$Zone))

# remove - and / in the zone 
sol$Zone_sub <- sub("\\-.*", "", sol$Zone)
sol$Zone_sub <- sub("\\/.*", "", sol$Zone_sub)
sol$Zone_sub <- as.numeric(sol$Zone_sub)

sort(unique(sol$Zone_sub))

sol$ices_area <- if_else(sol$Zone_sub <= 300, "4c", 
                         if_else(sol$Zone_sub > 300 & sol$Zone_sub <= 500, "4b",
                                 if_else(sol$Zone_sub > 500 & sol$Zone_sub <= 800, "4a",
                                         if_else(sol$Zone_sub > 2000, "7a",
                                                 if_else(sol$Zone_sub > 1600 & sol$Zone_sub <= 1700, "8ab", 
                                                         if_else(sol$Zone_sub > 800 & sol$Zone_sub <= 900, "7de",
                                                                 if_else(sol$Zone_sub > 900 & sol$Zone_sub <= 1000, "7fg",NA)))))))

#sol <- sol %>% filter(is.na(ices_area) == F)
#sol <- sol %>% filter(zone_uncertainty == "certain") #for now only consider the certain zones


# convert and add year_sample
sol$date <- dmy(sol$date)
sol$year_catch <- year(sol$date)


# select only data before 2004
sol <- sol %>% filter(year_catch < 2004)

# process year class
sort(unique(sol$jaarklas))

sol %>% filter(is.na(`jaarklas W`) == F | is.na(`jaarklas S`) == F) # 105 records
# for the mommment we will exclude these records

sol$jaarklas_new <- if_else(str_length(sol$jaarklas) != 3, #remove 3 digit numbers
                            sol$jaarklas, 
                            NA)

sol$jaarklas_new <- as.numeric(sol$jaarklas_new) #transform to numeric, strange format will be NA

# There is uncertainty in the jaarklas so here I am trying to get the most certain ones
# Examples of uncertainty:

# Assumption about the Year class
# 0-40 are the age
# 50-99 are the end of year -> 1950, 1999
sort(unique(sol$jaarklas_new))

sol$yearclass <- if_else(sol$jaarklas_new >= 41 & sol$jaarklas_new <= 100, sol$jaarklas_new + 1900, 
                         if_else(sol$jaarklas_new >= 1976 & sol$jaarklas_new <= 2000, sol$jaarklas_new,
                                 if_else(sol$jaarklas == 19996, 1996, 
                                         if_else(sol$jaarklas_new < 40, sol$year_catch - sol$jaarklas_new, NA))))

# set certainty for year class (select 28 to maybe include that fish)
sol$yearclass_certainty <- if_else(sol$jaarklas_new >= 1976 | sol$jaarklas_new <= 28, "certain", "uncertain")

# calculate fish age
sol$age <- sol$year_catch - sol$yearclass
sol$age_certainty <- sol$yearclass_certainty

write_csv(sol, file.path(dir_otolith_raw, "sol_database_pre2004_raw.csv"))

# 3.1.1. process data ---------------------------------------------------

## Process data to the format of ILVO SmartFish database

sol_pre04 <- sol %>% filter(ices_area %in% c("7a","4a","4b","4c","8ab","7de","7fg"), 
                            sol$zone_certainty == "certain",
                            yearclass >= 1957,
                            age >= 3,
                            lengte <= 100,           # remove extreme values
                            gewicht <= 2000,         # remove extreme values
                            sex == 2)                #female


sol_pre04 <- sol_pre04 %>% mutate(TRI.Year = year_catch, 
                                  TRI.Vessel = ship,
                                  
                                  HAU.IcesArea = ices_area,
                                  
                                  SAM.SpeciesFaoCode = FAA_species_code,
                                  SAM.SpeciesNameEnglish = "Common sole",
                                  
                                  SPE.Length = lengte*10, #before 2004 length is in cm
                                  SPE.Weight = gewicht,
                                  SPE.WeightGonads = gewicht_gonaden,
                                  SPE.Sex = "F",
                                  SPA.Number = nummer,
                                  SPA.Age = age,
                                  
                                  SAM.Date = date,
                                  SPA.Yearclass = yearclass,
                                  Yearclass.Certainty = yearclass_certainty,
                                  Ices.Zone = Zone_sub,
                                  UniqueID = unique_id
                                  )

sol_pre04 <- sol_pre04 %>% select(TRI.Year, TRI.Vessel, 
                                  HAU.IcesArea,
                                  SAM.SpeciesFaoCode, SAM.SpeciesNameEnglish,
                                  SPE.Length, SPE.Weight,
                                  SPE.Sex,
                                  SPA.Number,
                                  SPA.Age,
                                  SAM.Date,
                                  SPA.Yearclass,
                                  Yearclass.Certainty,
                                  Ices.Zone, 
                                  File, 
                                  UniqueID
                                  )

# There are duplicate UniquID (unknown reason by the time the data is made)
# so better to not included these duplicates

duplicate_sol_pre04 <- sol_pre04[duplicated(sol_pre04$UniqueID),]
list_duplicate_sol_pre04 <- unique(duplicate_sol_pre04$UniqueID)

sol_pre04 <- sol_pre04 %>% filter(!(UniqueID %in% list_duplicate_sol_pre04))

# 3.2. SOL POST 2004  ------------------------------------------------------

# Post 2004 data
sol_post04 <- sol_post04_raw %>% select(TripID, TRI.Year, TRI.Vessel, TRI.Code,
                                        HaulID, HAU.TimeZone, HAU.ShootTime, HAU.HaulTime, HAU.IcesArea,
                                        HAU.HaulLongitude, HAU.HaulLatitude, HAU.ShootLongitude, HAU.ShootLatitude, 
                                        SegmentID, 
                                        SampleID, SAM.SpeciesFaoCode, SAM.SpeciesNameEnglish, 
                                        SpecimenID, SPE.Length, SPE.Weight, SPE.WeightGonads, SPE.Sex, 
                                        SPA.Number, SPA.Age,
                                        TRI.DepartureDate, TRI.ReturnDate, 
                                        SAM.Observer, SAM.FateCategoryDescription, SPE.Sequence)



# Add/change certain fields
sol_post04 <- sol_post04 %>% mutate(SAM.Date = as_date(HAU.ShootTime),
                                    TRI.Date = paste(as_date(TRI.DepartureDate), as_date(TRI.ReturnDate), sep = "|"),
                                    SPA.Yearclass = TRI.Year - SPA.Age,
                                    Yearclass.Certainty = "certain",
                                    HAU.IcesArea = if_else(HAU.IcesArea %in% c("8a","8b"), "8ab", HAU.IcesArea))


sol_post04 <- sol_post04 %>% filter(SPE.Sex == "F",
                                    SPA.Age >= 3,
                                    HAU.IcesArea %in% c("4b","4c","7a","8ab"))
unique(sol_post04$HAU.IcesArea)

# 3.3. MERGE RAW DATA  -----------------------------------------------------

sol_all <- bind_rows(sol_pre04, sol_post04)

# Add SpecimenID to UniqueID
# Merge 4b and 4c
# Convert SPE.Length (before 2003 was cm)
sol_all <- sol_all %>% mutate(UniqueID = if_else(is.na(UniqueID) == T, SpecimenID, UniqueID),
                              HAU.IcesAreaGroup = if_else(HAU.IcesArea %in% c("4b","4c"), "4bc", HAU.IcesArea))

# 3.4. MISSING DATA ----
# Some issue occuring leading to the lost track of the script and mismatch between sol_all and the final file 13.10.22
# sol_all here has 62886 observations
# final file: sol_select_full (13/10/2022) 66060
# the missing observations include: 
## ifremer observations (will be added below)
## 1406 observation with no ices_area info in original file, but somehow have HAU.IcesArea in the final file

# missing observations
sol_final <- read_rds("./data/otolith/sol_select_full_22.10.13.rds")

## missing in raw data: obs in the final file but not in sol_all
sol_missing_raw <- anti_join(sol_final, sol_all, by = "UniqueID") %>% filter(TRI.Year < 2004)

# add missing observation to sol_all
sol_all <- bind_rows(sol_all, sol_missing_raw)

# save data 
write_csv(sol_all, file.path(dir_otolith_raw, "sol_database_1957-2020_female_age3+.csv"))
#sol_all <- read_csv(file.path("D:/OneDrive - UGent/data/WP1/otolith/","sol_database_1957-2020_female_age3+.csv"))

# select only certain fields for reading in SmartLab
sol_select <- sol_all %>% select(SAM.SpeciesFaoCode,
                                 HAU.IcesAreaGroup, HAU.IcesArea, TRI.Year, SPA.Yearclass, Yearclass.Certainty, SPA.Age, 
                                 TRI.Vessel, TRI.Code,
                                 TRI.Date, SAM.Date, UniqueID, 
                                 SAM.Observer, SAM.FateCategoryDescription, SPE.Sequence,
                                 SPA.Number, SPE.Length, SPE.Weight, SPE.WeightGonads, SPE.Sex)

# 3.6. CREATE Smartlab.Number --------

# previous attemps created 2 system of Smartlab.Number
# 1. values from 1 to N (but lost track of the script -> use final file 13.10.22 for this)
# 2. using SAM.Date/TRI.Date, SPA.Age, SPE.Weight/SPE.Length, SPE.Sequence/SPA.Number

## system 1 ----
sol_final_1 <- sol_final %>% filter(as.numeric(Smartlab.Number) < 2000) %>%
  select(UniqueID, Smartlab.Number)
sol_select_1 <- sol_select %>% 
  filter(UniqueID %in% sol_final_1$UniqueID) %>%
  left_join(sol_final_1, by = "UniqueID")

## system 2 ----
sol_select_2 <- sol_select %>% filter(!(UniqueID %in% sol_final_1$UniqueID))

# formula: SAM.Date_SPA.Number 
sol_select_2 <- sol_select_2 %>% 
  mutate(Smartlab.Number = if_else(TRI.Year >= 2004, 
                                   as.character(SPA.Number), 
                                   as.character(paste0(gsub("[[:punct:]]", "", as.character(SAM.Date)), SPA.Age, SPE.Weight, SPA.Number))))

# remove duplicate values (TRI.Year 1977)
duplicate <- sol_select_2[duplicated(sol_select_2$Smartlab.Number),] %>% filter(is.na(Smartlab.Number) == F)
sol_select_2 <- sol_select_2 %>% filter(!(Smartlab.Number %in% duplicate$Smartlab.Number))

# for Belgica A.962 - there is no SAM.Date due to HAU.Haul/Shoottime so we use TRI.Date instead
sol_select_2 <- sol_select_2 %>% 
  mutate(Smartlab.Number = if_else(is.na(SPA.Number) == T, 
                                   as.character(paste0(gsub("[[:punct:]]", "", substr((TRI.Date),0,10)),  SPA.Age, SPE.Length, SPE.Sequence)),
                                   Smartlab.Number))


# remove duplicate values (A.962 2005 2006)
duplicate <- sol_select_2[duplicated(sol_select_2$Smartlab.Number),]
sol_select_2 <- sol_select_2 %>% filter(!(Smartlab.Number %in% duplicate$Smartlab.Number))

## merge 2 system
sol_select <- bind_rows(sol_select_1, sol_select_2)

# 3.7. IFREMER DATA ------------------------------------------------------------
# # 1. read sol_select_full (13.10.22)
# # select UniqueID and from Status.Picture to the end - including Status.Read, Note, Consistency.test
# 
# sol_select_full_221013 <- read_csv(file.path("D:/OneDrive - UGent/data/WP1/otolith/ilvo_22.10.13","sol_select_full_221013.csv"))
# 
# # check number of otl pictured and read
# sum(sol_select_full_221013$Status.Picture, na.rm = T) # 1745
# sum(sol_select_full_221013$Status.Read, na.rm = T) # 703
# 
# # change SPE.Sex to character (wrong data format (logical) due to xlsx file)
# sol_select_full_221013 <- sol_select_full_221013 %>% mutate(SPE.Sex = as.character("F"),
#                                                             Smartlab.Number = as.character(Smartlab.Number))

# 2. Load and process ifremer data
sol_ifremer <- read_rds(file.path("./data/otolith", "sol_ifremer.rds"))

# filter sol_ifremer for 8AB, Age >= 3, and Sexe F
sol_ifremer <- sol_ifremer %>% filter(CIEM %in% c("8AB", "8A", "8B"), #There are a few (132) images in other regions 
                                      as.numeric(Age) >= 3, 
                                      Sexe == "F"
) 

# change field names sol_ifremer to be consistent with sol_select
names(sol_ifremer)
names(sol_select)

sol_ifremer <- sol_ifremer %>% mutate(SAM.SpeciesFaoCode  = Code_FAO,
                                      HAU.IcesAreaGroup   = "8ab",
                                      HAU.IcesArea        = if_else(CIEM == "8AB", "8ab",
                                                                    if_else(CIEM == "8A", "8a", "8b")),
                                      TRI.Year            = year,
                                      Yearclass           = TRI.Year - as.numeric(Age),
                                      Yearclass.Certainty = "certain",
                                      SPA.Age             = as.numeric(Age),
                                      TRI.Vessel          = Navire,
                                      TRI.Code            = NA,
                                      TRI.Date            = NA,
                                      SAM.Date            = as_date(Date),
                                      UniqueID            = Reference_PC,
                                      SAM.Observer        = NA,
                                      SAM.FateCategoryDescription = NA,
                                      SPE.Sequence        = NA,
                                      SPA.Number          = NA,
                                      SPE.Length          = as.numeric(Taille),
                                      SPE.Weight          = as.numeric(Poids),
                                      SPE.WeightGonads    = NA,
                                      SPE.Sex             = Sexe,
                                      Smartlab.Number = as.character(paste0(gsub("[[:punct:]]", "", as.character(SAM.Date)), 
                                                                            #SPA.Age, 
                                                                            #SPE.Length, 
                                                                            #SPE.Weight, 
                                                                            sapply(regmatches(Reference_PC, gregexpr("[[:digit:]]+", Reference_PC)), paste, collapse="") #extract number from Reference_PC
                                      )
                                      )
                                      )
# retain only new field names
sol_ifremer <- sol_ifremer %>% select(SAM.SpeciesFaoCode:Smartlab.Number)

# check duplicate
duplicate <- sol_ifremer[duplicated(sol_ifremer$Smartlab.Number),]

# 3.8. MERGE ILVO AND IFREMER DATA ----

sol_select_full <- bind_rows(sol_select, sol_ifremer)

# check duplicate
duplicate <- sol_select_full[duplicated(sol_select_full$Smartlab.Number),]

# save file
dir_otolith <- "./data/otolith"
write_csv(sol_select_full, file.path(dir_otolith,"sol_select_full.csv"))
write_xlsx(sol_select_full, file.path(dir_otolith,"sol_select_full.xlsx"))
write_rds(sol_select_full, file.path(dir_otolith,"sol_select_full.rds"))


