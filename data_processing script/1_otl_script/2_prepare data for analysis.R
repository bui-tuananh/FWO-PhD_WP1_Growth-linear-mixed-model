# ****************************************
# OTOLITH DATA PREPARATION FOR ANALYSIS **
# ****************************************

# INFO --------------------------------------------------------------------

# PROJECT: FWO PHD - WP1

# Prepare otolith data for analysis

# Input: 
## otolith data from ILVO (including IFREMER pictures) and WUR
### ILVO: Smartlab
### WUR: dutch_sole_otolith_increments_20090326.xls | in D:/OneDrive - UGent/data/WP1/otolith/@raw

## Temperature and Fishing data
### Temperature: isimip_sbt_ices.rds | in D:/OneDrive - UGent/data/Env/Sea Temperature_ISIMIP3_0.5deg_1850-2014
### Temperature data is processed from script: isimip_thetao_process.R

### Fishing data: ICES_StockAssessment_2021.csv | in D:/OneDrive - UGent/data/ICES
### Fishing data is collated from ICES stock assessment reports (North Sea, Irish Sea, Bay of Biscay)

# Output: otolith data with information of otolith growth, temperature, and fishing data
## the data has potential outliers but are remained to be removed in data exploration step

# 1. SETUP ----------------------------------------------------------------

# load packages
library(tidyverse) # process data frame
library(readr)     # process csv, rds
library(readxl)    # read xlxs
library(lubridate) # process time

# dir otl
dir_otl <- "D:/OneDrive - UGent/data/WP1/otolith"

# load function
dir_function <- "D:/OneDrive - UGent/data/WP1/otolith/function"

## function to get data from SmartLab
source(file.path(dir_function,"getdata_SmartLab.R"))

## function to transform SmartLab data into AnnulusDiameter and AnnulusDiameterIncrement
source(file.path(dir_function,"transformdata_SmartLab.R"))

# 2. LOAD AND PROCESS DATA ------------------------------------------------

## 2.1. ILVO DATA  ----------------------------------------------------------------

### 2.1.1. Extract data from Smartlab ---------------------------------------

# Set variables to get SmartLab data
lab <- 'OTL'  
project <- 'SEAWISE'
year <- 2022

# metadata info
DsOtl <- loadDataOtl(lab, project, year)

# reading info
DsAgeReadingDetail <- loadAgeReadingDetail(lab, project, year)

# link DsOtl and DsAgeReadingDetail via Outcome_ID -> otl_ilvo
otl_smartlab <- left_join(DsOtl, DsAgeReadingDetail, by = c("Outcome_ID" = "OutcomeID") ) %>% select(-LineIndex, -DotShape)
otl_smartlab <- otl_smartlab %>% filter(is.na(Sample_Number) == F, is.na(Outcome_ID) == F, Outcome_Result != 0) # remove NA values

# save otl_smartlab
write_rds(otl_smartlab, file.path(dir_otl, "otl_smartlab.rds"))

### 2.1.2. Process main data (otl_ilvo) -----------------------------------------------------

#### 2.1.2.1. Transform Smartlab data into actual measurements ------------------------

# 1. keep only data read by Tuan Anh (4bc, 7a) and Kelly (8ab)
## list experiment and 8ab sampletset 
list_exp <- grep("Rosa", unique(otl_smartlab$SampleSet_Code), value = TRUE) #list of experiment SampleSet_Code
list_8ab <- grep("8AB", unique(otl_smartlab$SampleSet_Code), value = TRUE)

otl_tbui <- otl_smartlab %>% filter(Outcome_LabTechnician == "CLO\\tbui", 
                                    !SampleSet_Code %in% list_8ab,
                                    !SampleSet_Code %in% list_exp) 

otl_kdiaz <- otl_smartlab %>% filter(Outcome_LabTechnician == "CLO\\kdiaz", 
                                     SampleSet_Code %in% list_8ab) 

otl_ilvo <- bind_rows(otl_tbui, otl_kdiaz)

# 2. transform SmartLab data into AnnulusDiameter and AnnulusDiameterIncrement 
otl_ilvo <- SmartDots_to_OtolithIncrement(otl_ilvo)

#### 2.1.2.2. Add metadata info and change field names ------------------------

# 1. add information from metadata
# note Sample_IcesDivison was set to 8A for all 8AB, for correct area need to link to SmartFish

# Add SAM.Date from metadata because 
# 1) Sample_CatchDate is in unfriendly date format dd/mm/yyyy (that needs transformation)
# 2) Sample_CatchDate has NA values due to NA HAU.HaulTime so use TRI.DepartureDate as sampling date SAM.Date
# 3) SAM.Date is used at Year and Month level so uncertainty in day is acceptable

metadata_otl_ilvo <- read_rds(file.path("D:/OneDrive - UGent/data/WP1/otolith/@processed","sol_select_full.rds"))
metadata_otl_ilvo <- metadata_otl_ilvo %>% filter(UniqueID %in% otl_ilvo$Sample_NumberExternal) 

# change date - NA SAM.Date use TRI.Departure date
metadata_otl_ilvo <- metadata_otl_ilvo %>% mutate(SAM.Date = if_else(is.na(SAM.Date) == T,
                                                                     as_date(sub("\\|.*", "", TRI.Date)),
                                                                     SAM.Date)) 

metadata_otl_ilvo <- metadata_otl_ilvo %>% select(SAM.SpeciesFaoCode:SPA.Age,
                                                  TRI.Date:UniqueID, SPE.Length:SPE.Sex)

# join info from metadata into otl_ilvo
otl_ilvo <- left_join(otl_ilvo, metadata_otl_ilvo, by = c("Sample_NumberExternal" = "UniqueID"))

# 2. arrange field names
## calculate reading age at capture
AgeAtCapture <- otl_ilvo %>% group_by(Sample_NumberExternal) %>% summarize(AgeAtCapture = max(DotIndex))
otl_ilvo <- left_join(otl_ilvo, AgeAtCapture)

## arrange field name
otl_ilvo <- otl_ilvo %>% mutate(FishID = Sample_NumberExternal,
                                SmartlabNumber = Sample_Number,
                                SpeciesFaoCode = SAM.SpeciesFaoCode,
                                OtolithProcessingMethod = Analysis_ProcessingMethod,
                                AQCode = Outcome_Quality,
                                Scale.pixelpermm = File_Scale,
                                IcesArea = HAU.IcesArea,
                                IcesAreaGroup = HAU.IcesAreaGroup,
                                TripDate = TRI.Date,
                                SamplingDate = SAM.Date,
                                SamplingYear = TRI.Year,
                                Cohort = SamplingYear - AgeAtCapture,
                                AgeAtCapture = AgeAtCapture,
                                Length.mm = SPE.Length,
                                Weight.g = SPE.Weight,
                                Sex = SPE.Sex,
                                Age = DotIndex,
                                AgeAtCapture.Database = SPA.Age,
                                GrowingYear = Cohort + Age - 1, #Growing year (Age 1 - Cohort 1970 -> growing period in year 1970)
                                AnnulusDiameter.um = round(AnnulusDiameter.um, 2),
                                AnnulusDiameterIncrement.um = round(AnnulusDiameterIncrement.um, 2),
                                OtolithWidth.um = round(OtolithWidth.um, 2),
                                Reader = if_else(Outcome_LabTechnician == "CLO\\tbui", "tbui", "kdiaz"),
                                DataSource = "ILVO")

otl_ilvo <- otl_ilvo %>% select(FishID:GrowingYear, AgeAtCapture, AnnulusDiameter.um:OtolithWidth.um, Reader, DataSource)

## change OtolithProcessingMethod 
### 1. more explicit: B&B -> broken/burned; S&S -> sectioned/stained
### 2. change ifremer otolith to: sectioned
list_ifremer <- unique(filter(otl_ilvo, str_detect(FishID, "RE|AL|CO") == TRUE)$FishID)

otl_ilvo <- otl_ilvo %>% mutate(OtolithProcessingMethod = if_else(OtolithProcessingMethod == "B&B", "broken/burned", "sectioned/stained"))
otl_ilvo <- otl_ilvo %>% mutate(OtolithProcessingMethod = if_else(FishID %in% list_ifremer, "sectioned", OtolithProcessingMethod))

#### 2.1.2.3. Pre-explore data ------------------------

# 3.0 check negative increment
neg_increment <- otl_ilvo %>% filter(AnnulusDiameterIncrement.um <= 0) %>% 
  select(IcesAreaGroup, SamplingDate, Scale.pixelpermm, SmartlabNumber, Age, AnnulusDiameter.um:OtolithWidth.um)

# 3.1. check for NA increment (due to error in aging annotation, or comment (same LONG/SHORT))
na_increment <- otl_ilvo %>% filter(is.na(AnnulusDiameter.um) == T) %>% 
  select(IcesAreaGroup, SamplingYear, SmartlabNumber)

# 3.2. check for NA File_Scale 
na_scale <- otl_ilvo %>% filter(is.na(Scale.pixelpermm) == T) %>% 
  select(IcesAreaGroup, SamplingYear, SmartlabNumber)

# 3.3. check for unusual scale
## list ifremer data
list_ifremer <- unique(filter(otl_ilvo, str_detect(FishID, "RE|AL|CO") == TRUE)$FishID)

# usual scale zeiss: 134, 172, 217, 219, 269, 270, 272, 345
# usual scale leica: 429, 437, 545, 593
sort(unique(filter(otl_ilvo, !FishID %in% list_ifremer)$Scale.pixelpermm)) 

# usual scale ifremer: 182, 212, 223, 224, 220, 226, 228, 228.5, 232, 260, 267, 279 
sort(unique(filter(otl_ilvo, FishID %in% list_ifremer)$Scale.pixelpermm)) 

# 3.4. check for NA Outcome_Quality
na_quality <- otl_ilvo %>% filter(is.na(AQCode) == T)

# 3.5. Check for outliner (due to wrong assignment of scale/reading error)
ggplot() +
  geom_point(data = otl_ilvo, aes(x = Length.mm, y = AnnulusDiameter.um), alpha = 0.5) + 
  facet_grid(~ Age)

# 1. few fish with very big growth
#age5 > 5000
otl_diameter_age5_gte5000 <- otl_ilvo %>% filter(Age == 5, AnnulusDiameter.um >= 5000)
otl_diameter_age5_gte5000 <- otl_ilvo %>% filter(FishID %in% otl_diameter_age5_gte5000$FishID)

#age3 > 4700
otl_diameter_age3_gte4700 <- otl_ilvo %>% filter(Age == 3, AnnulusDiameter.um >= 4700)
otl_diameter_age3_gte4700 <- otl_ilvo %>% filter(FishID %in% otl_diameter_age3_gte4700$FishID)

# 2. few fish with very small growth
#age4 < 2100
otl_diameter_age4_lt2100 <- otl_ilvo %>% filter(Age == 4, AnnulusDiameter.um < 2100)
otl_diameter_age4_lt2100 <- otl_ilvo %>% filter(FishID %in% otl_diameter_age4_lt2100$FishID)

#age9 < 3300
otl_diameter_age9_lt3300 <- otl_ilvo %>% filter(Age == 9, AnnulusDiameter.um < 3300)
otl_diameter_age9_lt3300 <- otl_ilvo %>% filter(FishID %in% otl_diameter_age9_lt3300$FishID)

# list of otl_ilvo to be checked
otl_ilvo_check <- bind_rows(otl_diameter_age5_gte5000, 
                            otl_diameter_age3_gte4700, 
                            otl_diameter_age4_lt2100, 
                            otl_diameter_age9_lt3300) 

ggplot() +
  geom_point(data = otl_ilvo, aes(x = Length.mm, y = AnnulusDiameter.um), alpha = 0.5) + 
  geom_point(data = otl_ilvo_check %>% filter(IcesAreaGroup == "8ab"), aes(x = Length.mm, y = AnnulusDiameter.um), color = "red", alpha = 0.5) + 
  facet_grid(~ Age)

## list 
otl_ilvo_check_list <- otl_ilvo_check %>%
  select(IcesAreaGroup, SamplingYear, SmartlabNumber, Scale.pixelpermm, OtolithWidth.um) %>%
  unique() %>% 
  arrange(IcesAreaGroup, SamplingYear, SmartlabNumber, OtolithWidth.um)

# note: otoliths to be checked for scale/measurement error
### 4bc
#### 1977 - 197710176660684 | correct
#### 1978 - 197802276535135 | correct
#### 1979 - 197903265553197 | correct
#### 1979 - 197903269696210 | correct
#### 1986 - 19860113651867 | correct
#### 1991 - 199105139465354 |  correct
#### 2005 - 175451 | correct

### 7a 
#### 1977 - 197703146555122 | remove (the edge is not flat -> difficult to see the rings)
#### 1989 - 198904036447288 | correct
#### 2004 - 168875 | wrong scale
#### 2004 - 169032 | wrong scale
#### 2004 - 169162 | wrong scale
#### 2004 - 169231 | wrong scale
#### 2004 - 169212 | wrong scale
#### 2004 - 169214 | remove (cannot see 1st, 2nd rings)
#### 2004 - 169215 | wrong scale

### 8ab
#### 1989 - 19890203325973 | remove (not flat surface)
#### 2013 - 32399 | correct
#### 2013 - 20130214134420005 | seems ok
#### 2013 - 20130215134330006 | seems ok
#### 2013 - 201302151344100003 | seems ok

#### 2.1.2.4. Save file ---------------------------------------------------------------
write_rds(otl_ilvo, file.path(dir_otl, "otl_ilvo.rds"))
# remain all suspected outliters, consider removing at exploration stage later (after check for observation error)

# re-read otl_ilvo for further steps
otl_ilvo <- read_rds(file.path(dir_otl, "otl_ilvo.rds"))

### 2.1.3. Process 8ab consistency test data (otl_ilvo_consistency) --------------------------------
# Get data
## tbui data 
otl_ilvo_consistency_tbui <- otl_smartlab %>% filter(Outcome_LabTechnician == "CLO\\tbui", 
                                                     SampleSet_Code %in% list_8ab)
otl_ilvo_consistency_tbui <- SmartDots_to_OtolithIncrement(otl_ilvo_consistency_tbui) #calculate measurement
### calculate reading age at capture
AgeAtCapture_tbui <- otl_ilvo_consistency_tbui %>% group_by(Sample_NumberExternal) %>% summarize(AgeAtCapture = max(DotIndex))
otl_ilvo_consistency_tbui <- left_join(otl_ilvo_consistency_tbui, AgeAtCapture_tbui)

## kdiaz data
otl_ilvo_consistency_kdiaz <- otl_smartlab %>% filter(Outcome_LabTechnician == "CLO\\kdiaz",
                                                      Sample_NumberExternal %in% otl_ilvo_consistency_tbui$Sample_NumberExternal)
otl_ilvo_consistency_kdiaz <- SmartDots_to_OtolithIncrement(otl_ilvo_consistency_kdiaz) #calculate measurement
### calculate reading age at capture
AgeAtCapture_kdiaz <- otl_ilvo_consistency_kdiaz %>% group_by(Sample_NumberExternal) %>% summarize(AgeAtCapture = max(DotIndex))
otl_ilvo_consistency_kdiaz <- left_join(otl_ilvo_consistency_kdiaz, AgeAtCapture_kdiaz)

## merge data
otl_ilvo_consistency <- bind_rows(otl_ilvo_consistency_tbui, otl_ilvo_consistency_kdiaz)

# Process data
## join info from metadata 
otl_ilvo_consistency <- left_join(otl_ilvo_consistency, metadata_otl_ilvo, by = c("Sample_NumberExternal" = "UniqueID"))

# 6. arrange field names
## arrange field name
otl_ilvo_consistency <- otl_ilvo_consistency %>% 
  mutate(FishID = Sample_NumberExternal,
         SmartlabNumber = Sample_Number,
         SpeciesFaoCode = SAM.SpeciesFaoCode,
         OtolithProcessingMethod = Analysis_ProcessingMethod,
         AQCode = Outcome_Quality,
         Scale.pixelpermm = File_Scale,
         IcesArea = HAU.IcesArea,
         IcesAreaGroup = HAU.IcesAreaGroup,
         TripDate = TRI.Date,
         SamplingDate = SAM.Date,
         SamplingYear = TRI.Year,
         Cohort = SamplingYear - AgeAtCapture,
         AgeAtCapture = AgeAtCapture,
         Length.mm = SPE.Length,
         Weight.g = SPE.Weight,
         Sex = SPE.Sex,
         Age = DotIndex,
         AgeAtCapture.Database = SPA.Age,
         GrowingYear = Cohort + Age - 1, #Growing year (Age 1 - Yearclass 1970 -> growing period in year 1970)
         AnnulusDiameter.um = round(AnnulusDiameter.um, 2),
         AnnulusDiameterIncrement.um = round(AnnulusDiameterIncrement.um, 2),
         OtolithWidth.um = round(OtolithWidth.um, 2),
         Reader = Outcome_LabTechnician,
         DataSource = "ILVO")

otl_ilvo_consistency <- otl_ilvo_consistency %>% select(FishID:GrowingYear, AgeAtCapture, AnnulusDiameter.um:OtolithWidth.um, Reader, DataSource)

## change OtolithProcessingMethod 
### 1. more explicit: B&B -> broken/burned; S&S -> sectioned/stained
### 2. change ifremer otolith to: sectioned
otl_ilvo_consistency <- otl_ilvo_consistency %>% mutate(OtolithProcessingMethod = if_else(OtolithProcessingMethod == "B&B", "broken/burned", "sectioned/stained"))
otl_ilvo_consistency <- otl_ilvo_consistency %>% mutate(OtolithProcessingMethod = if_else(FishID %in% list_ifremer, "sectioned", OtolithProcessingMethod))

## change Reader
otl_ilvo_consistency <- otl_ilvo_consistency %>% mutate(Reader = if_else(Reader == "CLO\\tbui", "tbui", "kdiaz"))

# Save file
write_rds(otl_ilvo_consistency, file.path(dir_otl, "otl_ilvo_consistency.rds"))

### 2.1.4. Process re-aging data (otl_ilvo_reage) --------------------------------------------

# load data
otl_ilvo_reage <- otl_smartlab %>% filter(Outcome_LabTechnician == "CLO\\imaertens")

# process data
## join info from metadata 
otl_ilvo_reage <- left_join(otl_ilvo_reage, metadata_otl_ilvo, by = c("Sample_NumberExternal" = "UniqueID"))

# 6. arrange field names
## calculate reading age at capture
AgeAtCapture <- otl_ilvo_reage %>% group_by(Sample_NumberExternal) %>% summarize(AgeAtCapture = max(DotIndex))
otl_ilvo_reage <- left_join(otl_ilvo_reage, AgeAtCapture)

## arrange field name
otl_ilvo_reage <- otl_ilvo_reage %>% 
  mutate(FishID = Sample_NumberExternal,
         SmartlabNumber = Sample_Number,
         SpeciesFaoCode = SAM.SpeciesFaoCode,
         OtolithProcessingMethod = Analysis_ProcessingMethod,
         AQCode = Outcome_Quality,
         Scale.pixelpermm = File_Scale,
         IcesArea = HAU.IcesArea,
         IcesAreaGroup = HAU.IcesAreaGroup,
         TripDate = TRI.Date,
         SamplingDate = SAM.Date,
         SamplingYear = TRI.Year,
         Cohort = SamplingYear - AgeAtCapture,
         AgeAtCapture = AgeAtCapture,
         Length.mm = SPE.Length,
         Weight.g = SPE.Weight,
         Sex = SPE.Sex,
         Age = DotIndex,
         AgeAtCapture.Database = SPA.Age,
         GrowingYear = Cohort + Age - 1, #Growing year (Age 1 - Yearclass 1970 -> growing period in year 1970)
         Reader = Outcome_LabTechnician,
         DataSource = "ILVO")

otl_ilvo_reage <- otl_ilvo_reage %>% select(FishID:GrowingYear, AgeAtCapture, Reader, DataSource)

## change OtolithProcessingMethod 
### 1. more explicit: B&B -> broken/burned; S&S -> sectioned/stained
### 2. change ifremer otolith to: sectioned
otl_ilvo_reage <- otl_ilvo_reage %>% mutate(OtolithProcessingMethod = if_else(OtolithProcessingMethod == "B&B", "broken/burned", "sectioned/stained"))
otl_ilvo_reage <- otl_ilvo_reage %>% mutate(OtolithProcessingMethod = if_else(FishID %in% list_ifremer, "sectioned", OtolithProcessingMethod))

## merge data
otl_ilvo_reage$Scale.pixelpermm <- as.numeric(otl_ilvo_reage$Scale.pixelpermm) #somehow the scale was character
otl_ilvo_reage <- bind_rows(otl_ilvo, otl_ilvo_reage)

## change Reader
otl_ilvo_reage <- otl_ilvo_reage %>% mutate(Reader = if_else(Reader == "CLO\\imaertens", "imaertens", Reader))

# save file
write_rds(otl_ilvo_reage, file.path(dir_otl, "otl_ilvo_reage.rds"))


## 2.2. WUR DATA ---------------------------------------------------------------------

### 2.2.1. Load and Process data ---------------------------------------------------

## LOAD DATA
# otolith data
dir_otl_wur <- "D:/OneDrive - UGent/data/WP1/otolith/@raw"
otl_wur_raw <- read_excel(path = file.path(dir_otl_wur, "dutch_sole_otolith_increments_20090326.xls"), sheet = "sol")

## PROCESS DATA 

#Data note
# edge value: ror - ring op rand → the last annulus diameter = otolith width (the same as ILVO 1st quarter) and the value of last annulus diameter is not noted. So backcal_age = 6 has 5 annulus diameter and the last diameter is otolith width
# edge value: blank -  all diameters of readable annulus noted, no change needed. 
# age > number of annulus diameters: the diameters of the latest years of life were difficult to measure, so only the diameters of the last readable annulus and the final otolith width. - follow the approach mentioned in Millner and Whitting 1996 "In old specimens it was not possible to separate the outer annual growth rings and measurements were made to the last reliable ring and to the outside edge"

#logic to process data 
#edge = ror 
#   - age = increment + 1 → add last diameter as otolith width
#   - age > increment + 1 (large age) → no change
#other: → nochange

# 1. up pivot to have each row is 1 increment record
otl_wur <- pivot_longer(otl_wur_raw, cols = starts_with("_"), names_to = "age_increment", values_to = "diameter")

otl_wur <- otl_wur %>% filter(is.na(diameter) == F) # remove NA diameter values 
otl_wur <- otl_wur %>% mutate(age_increment = as.numeric(str_remove(age_increment,"_"))) #remove _ in age_increment

# 3. Extract increment from diameter values 
# late months (sep - dec): need to reduce 1 age
# early months (jan - april): 
# - age = age_increment: keep age -> need to add increment
# - age > age_increment + 1: no change

otl_wur_list <- unique(otl_wur$image_ID) # create list of image_ID
otl_wur_increment = otl_wur[FALSE,]        # create an empty df

for (i in 1:length(otl_wur_list)) {
  
  print(paste("processing", otl_wur_list[i], sep = " "))
  
  otl <- otl_wur %>% filter(image_ID == otl_wur_list[i])
  
  if(is.na(unique(otl$edge)) == F) {
    print("edge = ror - ring op rand")
    
    if(unique(otl$backcal_age) == max(otl$age_increment) + 1) {
      print("age = age_increment + 1 -> add last diameter width = otolith width")
      a <- nrow(otl)
      otl[a+1,] <- otl[a,]
      otl[a+1,]$age_increment <- unique(otl$backcal_age)
      otl[a+1,]$diameter <- unique(otl$width)
    }
    
    else {
      print("age > age_increment + 1 -> no change")
      otl <- otl
    }
  } 
  
  else {
    print("edge = blank -> no change")
    }
  
  otl$increment <- c(otl$diameter[1], diff(otl$diameter))
  otl_wur_increment <- rbind(otl_wur_increment, otl)
  
   
}

# 4. arrange field names
otl_wur <-  otl_wur_increment %>%  mutate(FishID = image_ID,
                                          SmartlabNumber = NA,
                                          SpeciesFaoCode = "SOL",
                                          OtolithProcessingMethod = "sectioned/stained",
                                          AQCode = NA,
                                          Scale.pixelpermm = NA,
                                          IcesArea = if_else(lat_deg_min <= 53.5, "4c", "4b"), 
                                          IcesAreaGroup = "4bc",
                                          TripDate = NA,
                                          SamplingDate = as_date(paste(Year, month, day, sep = "-")),
                                          SamplingYear = Year,
                                          Cohort = yc,
                                          AgeAtCapture = backcal_age,
                                          Length.mm = length*10,
                                          Weight.g = weight, #weight otl_wur has problem (160 fish with weight < 250, even very old/long fish)
                                          Sex = sex,
                                          Age = age_increment,
                                          GrowingYear = Cohort + Age - 1, #Growing year (Age 1 - Yearclass 1970 -> growing period in year 1970)
                                          AnnulusDiameter.um = round(diameter, 2),
                                          AnnulusDiameterIncrement.um = round(increment, 2),
                                          OtolithWidth.um = round(width, 2),
                                          DataSource = "WUR")
                                
otl_wur <- otl_wur %>% 
  select(FishID:DataSource) %>% 
  filter(Sex == "F", Cohort >= 1957) #only female + yearclass from 1957 for this study

### 2.2.2. Pre-explore data ------------------------------------------------------------
#### 2.2.2.1. OtolithWidth.um -------------------------------------------------------------

# OtolithWidth.um vs AgeAtCapture

ggplot(data = otl_wur, aes(x = AgeAtCapture, y = OtolithWidth.um)) + 
  geom_point() 

# 2 otoliths have very small width (< 1500) 
otl_width_lt1500 <-  otl_wur %>% filter(OtolithWidth.um < 1500)

ggplot() +
  geom_point(data = otl_wur, aes(x = Length.mm, y = AnnulusDiameter.um), alpha = 0.5) + 
  geom_point(data = otl_width_lt1500, aes(x = Length.mm, y = AnnulusDiameter.um), alpha = 0.5, color = "red") +
  facet_grid(~ Age)

# width < diamater -> maybe typo 
# the growth curves looks fine so does not need to be removed

#### 2.2.2.2. AnnulusDiameter.um  -------------------------------------------------------------

# AnnulusDiameter.um vs Age and Length.mm (both ilvo and wur data)
ggplot() +
  geom_point(data = otl_ilvo, aes(x = Length.mm, y = AnnulusDiameter.um), alpha = 0.5, color = "grey") +
  geom_point(data = otl_wur, aes(x = Length.mm, y = AnnulusDiameter.um), alpha = 0.5, color = "black") + 
  facet_grid(~ Age)

# 1. 1 otolith with very small dimater at age 1
otl_diameter_age1_lt200 <- otl_wur %>% filter(AnnulusDiameter.um < 200)
otl_diameter_age1_lt200 <- otl_wur %>% filter(FishID %in% otl_diameter_age1_lt200$FishID)

# 2. 1 otolith with AnnulusDiameter.um < 1100 at age 2
otl_diameter_age2_lt1100 <- otl_wur %>% filter(Age == 2, AnnulusDiameter.um < 1100)
otl_diameter_age2_lt1100 <- otl_wur %>% filter(FishID %in% otl_diameter_age2_lt1100$FishID)

# 3. 2 otoliths with AnnulusDiameter.um < 600 at age 1  
otl_diameter_age1_lt600 <- otl_wur %>% filter(Age == 1, AnnulusDiameter.um < 600)
otl_diameter_age1_lt600 <- otl_wur %>% filter(FishID %in% otl_diameter_age1_lt600$FishID)

# 4. 1 otolith with AnnulusDiameter.um > 5200 at age 7
otl_diameter_age7_gte5200 <- otl_wur %>% filter(Age == 7, AnnulusDiameter.um >= 5200)
otl_diameter_age7_gte5200 <- otl_wur %>% filter(FishID %in% otl_diameter_age7_gte5200$FishID)

# list of otl_wur to be checked
otl_wur_check <- bind_rows(otl_diameter_age1_lt200, 
                           otl_diameter_age2_lt1100, 
                           otl_diameter_age1_lt600, 
                           otl_diameter_age7_gte5200) 

ggplot() +
  geom_point(data = otl_ilvo, aes(x = Length.mm, y = AnnulusDiameter.um), alpha = 0.5, color = "grey") +
  geom_point(data = otl_wur, aes(x = Length.mm, y = AnnulusDiameter.um), alpha = 0.5, color = "black") + 
  geom_point(data = otl_wur_check, aes(x = Length.mm, y = AnnulusDiameter.um), color = "red", alpha = 0.5) + 
  facet_grid(~ Age)

## list 
otl_wur_check_list <- otl_wur_check %>%
  select(FishID) %>%
  unique() %>% 
  arrange(FishID)

#### 2.2.2.3. AnnulusDiameterIncrement.um  -------------------------------------------------------------

# increment distribution -> normal - reduce eventually by age
ggplot() + geom_histogram(data = otl_wur, aes(x = AnnulusDiameterIncrement.um), binwidth = 30) +
  facet_wrap(~ Age) # by age

# AnnulusDiameterIncrement.um vs Age
ggplot() + 
  geom_point(data = otl_ilvo, aes(x = Age, y = AnnulusDiameterIncrement.um), alpha = 0.5, color = "grey") +
  geom_point(data = otl_wur, aes(x = Age, y = AnnulusDiameterIncrement.um), alpha = 0.5, color = "black") +
  facet_grid(DataSource ~ IcesAreaGroup) 

# similar in 4bc in ILVO and WUR

#### 2.2.2.4. Length.g -----------------------------------------------------------

#Length.mm vs OtolithWidth.um
ggplot() + 
  geom_point(data = otl_ilvo, aes(x = OtolithWidth.um, y = Length.mm, color = DataSource), alpha = 0.5, color = "grey") + 
  geom_point(data = otl_wur, aes(x = OtolithWidth.um, y = Length.mm, color = DataSource), alpha = 0.5, color = "black") + 
  facet_grid(~ IcesAreaGroup) 

# The relationship between length and otolith width is similar in 4bc in ILVO and WUR data (except for the suspected outliers)

#### 2.2.2.5. Weight.g ------------------------------------------------------------

# Weight.g vs Length.mm

ggplot() + 
  geom_point(data = otl_wur, aes(x = Length.mm, y = Weight.g), alpha = 0.5) 

# WUR: 4bc - small weight (< 250) at length > 300 mm
otl_weight_lt250 <- otl_wur %>% filter(Weight.g < 250, Length.mm > 300)
otl_weight_lt250$Weight.g_new <- otl_weight_lt250$Weight.g*10

ggplot() + 
  geom_point(data = otl_wur %>% filter(!(FishID %in% otl_weight_lt250$FishID)), aes(x = Length.mm, y = Weight.g), color = "grey") +
  geom_point(data = otl_weight_lt250, aes(x = Length.mm, y = Weight.g), color = "black") +
  geom_point(data = otl_weight_lt250, aes(x = Length.mm, y = Weight.g_new), color = "red") 

#Very likely that there was a typo in weight and the weight should be *10 for those with very small weight
otl_wur <- otl_wur %>% mutate(Weight.g = if_else(FishID %in% otl_weight_lt250$FishID, Weight.g*10, Weight.g))

### 2.2.3. Save file --------------------------------------------------------

write_rds(otl_wur, file.path(dir_otl, "otl_wur.rds"))
# remain all suspected outliters, consider removing at exploration stage later (after check for observation error)

## 2.3. MERGE ILVO AND WUR DATA ----------------------------------------------------
# re read the files
otl_ilvo <- read_rds(file.path(dir_otl, "otl_ilvo.rds"))
otl_wur <- read_rds(file.path(dir_otl, "otl_wur.rds"))

otl_full <- bind_rows(otl_ilvo, otl_wur)

## 2.4. REMOVE INCOMPLETE GROWTH INCREMENT ----------------------------------------------------

# incomplete growth increments are increment that are the increment in a year which is 
# measured but not completed as the winter ring of that growth year

# incomplete growth increments are the last increment in fish sampled in the first quarter (Jan - March)
# in April, it is very likely that the increment is incomplete so will not be included in analysis
# in May and June, it is ambiguous so the otoliths will be checked carefully at image reading stage
# it is also assumed that growth is completed in May and June

# workflow to remove incomplete growth increment
# create increment id -> extract incomplete growth -> remove incomplete growth

# 1. create increment id
otl_full <- otl_full %>% rowid_to_column("IncrementID")

# 2. extract incomplete growth increments
otl_incomplete <- otl_full %>% filter(month(SamplingDate) <= 4, AgeAtCapture == Age) 

# 3. remove incomplete growth
otl_full <- otl_full %>% filter(!(IncrementID %in% otl_incomplete$IncrementID)) 
otl_full <- otl_full %>% select(-IncrementID)

## 2.5. ADD TEMPERATURE AND FISHING DATA ----------------------------------------------------

### 2.5.1. TEMPERATURE (SEA BOTTOM TEMPERATURE) -----------------------------
# ISIMIP Sea Bottom Temperature - processed in isimip_thetao_process.R

dir_sbt <- "D:/OneDrive - UGent/data/Env/Sea Temperature_ISIMIP3_0.5deg_1850-2014"
sbt_ices <- read_rds(file.path(dir_sbt, "isimip_sbt_ices.rds"))

# pre-process sbt_ices to merge 
sbt_ices <- sbt_ices %>% group_by(IcesArea, Year) %>% summarize(SeaBottomTemperature.degC = mean(isimip_sbt, na.rm = T))

# add SeaBottomTemperature to otl data
otl_full <- left_join(otl_full, sbt_ices, by = c("IcesAreaGroup" = "IcesArea", "GrowingYear" = "Year"))

### 2.5.2. FISHING MORTALITY, STOCK BIOMASS -----------------------------

# load ices data
dir_ices <- "D:/OneDrive - UGent/data/ICES"

ices_4 <- read_excel(file.path(dir_ices, "WGNSSK 2021_4_summary.xlsx")) %>% mutate(IcesAreaGroup = "4bc")
ices_7a <- readxl::read_excel(file.path(dir_ices, "WGCSE2021_7a_summary.xlsx")) %>% mutate(IcesAreaGroup = "7a")
ices_8ab <- readxl::read_excel(file.path(dir_ices, "WGBIE 2021_8ab_summary.xlsx")) %>% mutate(IcesAreaGroup = "8ab")

ices <- bind_rows(ices_4, ices_7a, ices_8ab)

# rename filds
ices <- ices %>% 
  mutate(SpawningStockBiomass.1000t = SSB_t/1000,
         FishingMortality = as.numeric(Fbar),
         GrowingYear = Year) %>%
  select(IcesAreaGroup, GrowingYear, SpawningStockBiomass.1000t, FishingMortality)

# Add FishingMortality and SpawningStockBiomass to otl_full 
otl_full <- left_join(otl_full, ices, by = c("IcesAreaGroup", "GrowingYear"))

# 3. SAVE DATA ---------------------------------------------------------------
write_rds(otl_full, file.path(dir_otl, "otl_full.rds"))

