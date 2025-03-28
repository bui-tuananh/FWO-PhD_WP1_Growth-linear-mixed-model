# Process IFREMER metadata

# Workflow
# 1. Create a English description field 
# 2. Check duplicate in UniqueID - Reference_PC
# 3. Create list of pictures by years - change Reference_PC data by year (certain years has differnece: 17 vs 2017)
# 4. Filter metadata by list of pictures


# 1. SETUP ----------------------------------------------------------------

library(tidyverse)  # process data frame
library(readr)      # process rds, csv 
library(lubridate)  # process time
library(readxl)     # read xlsx

# 2. LOAD DATA ------------------------------------------------------------

# IFREMER metadata
dir <- "D:/OneDrive - UGent/data/WP1/otolith/ifremer"
ifremer_raw <- readxl::read_xlsx(file.path(dir, "french data sole otolith.xlsx"))

str(ifremer_raw)
names(ifremer_raw)

# reference dataset
otl <- read_rds(file.path("D:/OneDrive - UGent/data/WP1/otolith/@processed/sol_select_full.rds"))
names(otl)

# certain fields that need to be asked
names(ifremer_raw)


# 3. PROCESS DATA ---------------------------------------------------------


# 3.1. Create a English description field --------------------------------

# Navire            - Vessel
# Engin             - Gear
# Code_Espece       - SpeciesCode
# Code_Fao          - SpeciesFaoCode
# Nom_Scientifique  - SpeciesName
# Type_Longueur     - LengthType
# Increment 
# Unite_Taille      - LengthUnit
# Unit_Poids        - WeightUnit
# Presentation      - Presentation (Entier, Vide - Whole, ?)
# Maturite_Echelle  - MaturityScale
# Type_PC           
# Preparation_PC    - OtolithProcessingMethod (sectioned, sectioned/stained, burnt whole)
# Poids_PC          
# year              - SamplingYear
# Date              - SamplingDate
# Trimestre         - SamplingQuarter
# Numero_Trait
# Lieu 
# Zone              
# CIEM              - IcesArea
# Taille            - Length.mm
# Poids             - Weight.g
# Sexe              - Sex
# Maturite          - Maturity
# Age               - AgeAtCapture
# Reference_Prelevement
# Reference_PC      - UniqueID
# Observations


# 3.2. Check duplicate in UniqueID - Reference_PC -------------------------------------------------------------------

UniqueID <- unique(ifremer_raw$Reference_PC) # no duplicate in UniqueID - Reference_PC

# 3.3. Create list of pictures by year -----------------------------------
# change Reference_PC data by year (certain years has differnece: 17 vs 2017)

# list year
list_year <- as.numeric(list.files(dir))
list_year <- list_year[is.na(list_year) == F]

# list of pictures by year
ifremer_pic <- tibble(NULL)

for (i in 1:length(list_year)) {
  
  print(paste0("Processing year ", list_year[i]))
  
  df_temp <- tibble(SamplingYear = list_year[i],
                    UniqueID = list.files(file.path(dir, list_year[i])))
  
  ifremer_pic <- bind_rows(ifremer_pic, df_temp)
  
  rm(df_temp)
}

ifremer_pic <- ifremer_pic %>% mutate(UniqueID = str_remove(UniqueID, ".tif"))

#remove a names of folders containing picture of unclear metadata (e.g. pictures in folder 2009 but has metadata of 2008) 
ifremer_pic <- ifremer_pic %>% filter(!UniqueID %in% c("maybe 2008", "maybe 2011")) 

# 3.4. Filter metadata by list of pictures --------------------------------

ifremer_filter <- ifremer_raw %>% filter(Reference_PC %in% ifremer_pic$UniqueID)

metadata <- ifremer_filter %>% group_by(year) %>% summarize(no_metadata = n())
picture <- ifremer_pic %>% group_by(SamplingYear) %>% summarize(no_picture = n())

check <- left_join(metadata, picture, by = c("year" = "SamplingYear"))

# there is mismatch between metadata (Reference_PC) and picture name (UniqueID) in 2011-2017 (except for 2016)

# 2011 --------------------------------------------------------------------

data11 <- ifremer_raw %>% filter(year == 2011)

# 1. check all names start with AL
data11_AL <- data11 %>% filter(str_detect(Reference_PC, regex("AL", ignore_case = TRUE)) == T )
# no finding of AL_solper + images seem to include only small ages (<= 3) -> ignore AL_solper 

# difference: metadata "_", picture "-" -> replace "_" by "-" 
data11_AL <- data11_AL %>% mutate(Reference_PC = str_replace_all(Reference_PC, "_", "-"))

data11_AL_filter <- data11_AL %>% filter(Reference_PC %in% ifremer_pic$UniqueID)
# no match 

# 2. check all names start with RE
data11_RE <- data11 %>% filter(str_detect(Reference_PC, regex("RE", ignore_case = TRUE)) == T )

# pictures: RE_11_b79, RE_12_b80
data11_RE %>% filter(str_detect(Reference_PC, regex("RE_11_b79", ignore_case = TRUE)) == T )
data11_RE %>% filter(str_detect(Reference_PC, regex("RE_12_b80", ignore_case = TRUE)) == T )
# no match 

# unmatch pictures
pic11_unmatch <- ifremer_pic %>% filter(SamplingYear == 2011)
pic_unmatch <- pic11_unmatch

# 2012 --------------------------------------------------------------------

data12 <- ifremer_raw %>% filter(year == 2012)

# 1. check all names start with AL
data12_AL <- data12 %>% filter(str_detect(Reference_PC, regex("AL", ignore_case = TRUE)) == T )

# difference: 
# metadata "_", picture "-" -> replace "_" by "-" 
# metadata ends with 3 digits (eg 005), picture ends with 4 digits (eg 0005) -> replace O-0 with O-00

data12_AL <- data12_AL %>% mutate(Reference_PC = str_replace_all(Reference_PC, "_", "-"))
data12_AL <- data12_AL %>% mutate(Reference_PC = str_replace_all(Reference_PC, "O-0", "O-00"))

data12_AL_filter <- data12_AL %>% filter(Reference_PC %in% ifremer_pic$UniqueID) #337 obs
ifremer_pic %>% filter(SamplingYear == 2012, str_detect(UniqueID, regex("AL", ignore_case = TRUE)) == T ) #337 obs
# all names with AL matches

# add to ifremer_filter
ifremer_filter <- ifremer_filter %>% bind_rows(data12_AL_filter)

# 2. check all names start with RE
data12_RE <- data12 %>% filter(str_detect(Reference_PC, regex("RE", ignore_case = TRUE)) == T )

# only check for names that do not match 
pic12_RE_match <- ifremer_filter %>% filter(year == 2012, str_detect(Reference_PC, regex("RE", ignore_case = TRUE)) == T) #matched names
data12_RE <- data12_RE %>% filter(!Reference_PC %in% pic12_RE_match$Reference_PC) #unmatched names

#unmatched names
pic12_RE_unmatch <- ifremer_pic %>% filter(SamplingYear == 2012, 
                                        !UniqueID %in% pic12_RE_match$Reference_PC, 
                                        str_detect(UniqueID, regex("RE", ignore_case = TRUE)) == T)

sort(unique(data12_RE$Reference_PC))
sort(unique(pic12_RE_unmatch$UniqueID)) #299 obs

# no match 

# unmatch pictures
pic_unmatch <- pic_unmatch %>% bind_rows(pic12_RE_unmatch)

# 2013 --------------------------------------------------------------------

data13 <- ifremer_raw %>% filter(year == 2013)
data13 <- data13 %>% filter(!Reference_PC %in% ifremer_filter$Reference_PC) #unmatched names

# 1. check all names start with RE
data13_RE <- data13 %>% filter(str_detect(Reference_PC, regex("RE", ignore_case = TRUE)) == T )

#unmatched names
pic13_RE_match <- ifremer_filter %>% filter(year == 2013, str_detect(Reference_PC, regex("RE", ignore_case = TRUE)) == T) #matched names
pic13_RE_unmatch <- ifremer_pic %>% filter(SamplingYear == 2013, 
                                           !UniqueID %in% pic13_RE_match$Reference_PC, 
                                           str_detect(UniqueID, regex("RE", ignore_case = TRUE)) == T)

sort(unique(data13_RE$Reference_PC))
sort(unique(pic13_RE_unmatch$UniqueID)) #194 obs
# no match

# unmatch pictures
pic_unmatch <- pic_unmatch %>% bind_rows(pic13_RE_unmatch)

# 2014 --------------------------------------------------------------------

data14 <- ifremer_raw %>% filter(year == 2014)
data14 <- data14 %>% filter(!Reference_PC %in% ifremer_filter$Reference_PC) #unmatched names

# 1. check all names start with RE
data14_RE <- data14 %>% filter(str_detect(Reference_PC, regex("RE", ignore_case = TRUE)) == T )

#unmatched names
pic14_RE_match <- ifremer_filter %>% filter(year == 2014, str_detect(Reference_PC, regex("RE", ignore_case = TRUE)) == T) #matched names
pic14_RE_unmatch <- ifremer_pic %>% filter(SamplingYear == 2014, 
                                           !UniqueID %in% pic14_RE_match$Reference_PC, 
                                           str_detect(UniqueID, regex("RE", ignore_case = TRUE)) == T)

sort(unique(data14_RE$Reference_PC))
sort(unique(pic14_RE_unmatch$UniqueID)) #266 obs
# no match

# unmatch pictures
pic_unmatch <- pic_unmatch %>% bind_rows(pic14_RE_unmatch)

# 2015 --------------------------------------------------------------------

data15 <- ifremer_raw %>% filter(year == 2015)
data15 <- data15 %>% filter(!Reference_PC %in% ifremer_filter$Reference_PC) #unmatched names

# 1. check all names start with CO
data15_CO <- data15 %>% filter(str_detect(Reference_PC, regex("CO", ignore_case = TRUE)) == T )

#unmatched names
pic15_CO_match <- ifremer_filter %>% filter(year == 2015, str_detect(Reference_PC, regex("CO", ignore_case = TRUE)) == T) #matched names
pic15_CO_unmatch <- ifremer_pic %>% filter(SamplingYear == 2015, 
                                           !UniqueID %in% pic15_CO_match$Reference_PC, 
                                           str_detect(UniqueID, regex("CO", ignore_case = TRUE)) == T)

sort(unique(data15_CO$Reference_PC))
sort(unique(pic15_CO_unmatch$UniqueID)) #12 obs
# no match

# unmatch pictures
pic_unmatch <- pic_unmatch %>% bind_rows(pic15_CO_unmatch)


# 2017 --------------------------------------------------------------------

data17 <- ifremer_raw %>% filter(year == 2017)
data17 <- data17 %>% filter(!Reference_PC %in% ifremer_filter$Reference_PC) #unmatched names

# 1. check all names start with RE
data17_RE <- data17 %>% filter(str_detect(Reference_PC, regex("RE", ignore_case = TRUE)) == T )

#unmatched names
pic17_RE_match <- ifremer_filter %>% filter(year == 2017, str_detect(Reference_PC, regex("RE", ignore_case = TRUE)) == T) #matched names
pic17_RE_unmatch <- ifremer_pic %>% filter(SamplingYear == 2017, 
                                           !UniqueID %in% pic17_CO_match$Reference_PC, 
                                           str_detect(UniqueID, regex("RE", ignore_case = TRUE)) == T)

sort(unique(data17_RE$Reference_PC))
sort(unique(pic17_RE_unmatch$UniqueID)) #648 obs

# Mismatch:
# 1. metadata: 2017_b5 (51-55) - picture: 17_b5
# 2. metadata: b7 (71-77) - picture: B7 (B77)

data17_b5 <- data17_RE %>% filter(str_detect(Reference_PC, regex("_b5", ignore_case = TRUE)) == T)
data17_b5 <- data17_b5 %>% mutate(Reference_PC = str_replace(Reference_PC, "2017_b5", "17_b5"))
data17_b5_filter <- data17_b5 %>% filter(Reference_PC %in% pic17_RE_unmatch$UniqueID)

data17_b7 <- data17_RE %>% filter(str_detect(Reference_PC, regex("_b7", ignore_case = FALSE)) == T)
data17_b7 <- data17_b7 %>% mutate(Reference_PC = str_replace(Reference_PC, "_b7", "_B7"))
data17_b7_filter <- data17_b7 %>% filter(Reference_PC %in% pic17_RE_unmatch$UniqueID)

# all 2017 pictures are matched - 649 
nrow(pic17_RE_match) + nrow(data17_b5_filter) + nrow(data17_b7_filter) 

# add to ifremer_filter
ifremer_filter <- ifremer_filter %>% bind_rows(data17_b5_filter, data17_b7_filter)


# Summary -------------------------------------------

# Number of matched and mismatched pictures
metadata <- ifremer_filter %>% group_by(year) %>% summarize(no_metadata = n())
picture <- ifremer_pic %>% group_by(SamplingYear) %>% summarize(no_picture = n())

check <- left_join(metadata, picture, by = c("year" = "SamplingYear")) 
check <- check %>% mutate(no_mismatch = no_picture - no_metadata)

# 2011 is not in the list as there is no matched picture (721 pics)

# Save file
#write_rds(ifremer_filter, file.path(dir, "sol_ifremer.rds"))
#write_csv(ifremer_filter, file.path(dir, "sol_ifremer.csv"))

# Questions ---------------------------------------------------------------

# list of unmatch pictures
pic_unmatch
#writexl::write_xlsx(pic_unmatch, file.path(dir, "sol_ifremer_unmatch pictures.xlsx"))

# Is there a lack of metadata for those pictures or the pictures are named wrongly?

# 2011 - 721 pictures
# 2012 - 299
# 2013 - 194 (should be 190 but do not know why there is a difference here - maybe not important)
# 2014 - 266
# 2015 - 12

