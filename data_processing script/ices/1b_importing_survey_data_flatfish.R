# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++                                                          +++
# +++       IMPORTING DATRAS SURVEY DATA FOR BEAM TRAWLS       +++
# +++       FOR THE SPECIFIED PERIOD, AREAS AND SPECIES        +++
# +++                                                          +++
# +++         Jochen Depestele, Charlotte Van Moorleghem       +++
# +++                       October 2022                       +++
# +++                                                          +++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# DATRAS data were retrieved using the ICES DATRAS package, found here: https://github.com/ices-tools-prod/icesDatras 
# https://cran.r-project.org/web/packages/icesDatras/icesDatras.pdf
rm(list=ls())

# install libraries ----
library(data.table)
library(icesDatras)
library(dplyr)
library(ggplot2)
library(sf)
library(mapplots)
library(worms)
library(svMisc)
library(RColorBrewer)
library(plyr)

# ???
tosave = F # tosave=F when we're testing the script, i.e. without saving intermediate results
testrun = T # RUN THE SCRIPT ONLY FOR 2015 or another year 1993 (testing functionalities)
if(testrun)tosave=F
sp <- readRDS("C:/Users/jdepestele/OneDrive - ILVO/gitr/SEAwise/seawise_t5-4_juvsol7afg/input/bts/Worms_SpeciesList_BTS_in1993till2022_in7a7f7g.RDS")
setDT(sp)
# study_species <- sp[scientificname %in% c("Solea solea","Pleuronectes platessa","Merlangius merlangus","Gadus morhua")]$Valid_Aphia
study_species <- c(127160) # vector with valid Aphia IDs
study_area <- c("3.a.20", "4.a", "4.b", "4.c", "7.a", "7.d", "7.e", "7.f", "7.g", "7.h", "7.j.2", "8.a", "8.b") # c("7.a","7.f","7.g")
study_period <- c(1985L:(as.integer(substr(Sys.Date(),1,4))-1)) # DATRAS BTS data started in 1993 for 7afg, North Sea in 1985, Bay of Biscay in 2011
                                                                # get data until last year, as this year may not be uploaded yet
if(testrun)study_period <- 2011L # c(1985L, 2015L, 2021L) # study_period[1]
# ???

# SETUP ####
# set LOCAL folder directory ----
locdir <- sapply(strsplit(rstudioapi::getSourceEditorContext()$path,"/R"),"[[",1) # "C:/Users/jdepestele/OneDrive - ILVO/gitr/SEAwise/seawise_t5-4_juvsol7afg"

inputdir <- paste0(locdir,"/input")
# inputdirbts <- paste0(locdir,"/input/bts")
inputdirflats <- paste0(locdir,"/input/survey_flats")
outputdir <- paste0(locdir,"/output")

# Utilities ----
# Loading DATA from local drivers ====
# Setting formats
theme_set(theme_bw())
# Loading FUNCTIONS from local drivers ====
source("C:/Users/jdepestele/OneDrive - ILVO/gitr/follow_others/ICES/DATRAS_BTS/utilities/getDATRAS_mymodification.R")

# Get ICES divisions 7afg as an sf object
dir <- paste("//clo.be/Home/Home_d1/jdepestele/Documents/000_data/Areas/ICES/ICES_areas_shp")
ICESareas <- sf::st_read(dsn = dir, layer = "ICES_Areas_20160601_cut_dense_3857")
sort(unique(ICESareas$Area_27))
st_crs(ICESareas, parameters = TRUE)$units_gdal # unit= metre
crs_orig <- st_crs(ICESareas)
ICESareas <- ICESareas %>% dplyr::select(Area_27)
ICESareas <- st_transform(ICESareas,crs=4326)
ICESareas_study <- ICESareas[ICESareas$Area_27 %in% study_area,]
# ggplot(ICESareas_study)+ geom_sf()


# IMPORTING DATRAS DATA ####
# # Explore ICES surveys ----
# getSurveyList() # BTS = Beam Trawl Survey
#                 # "BTS-VIII" Beam Trawl Survey in the Bay of Biscay
#                 # "DYFS" Inshore Beam Trawl Survey
#                 # "SNS" Sole Net Survey => EU study in support of the Common Fishery Policy 95/025
# yrs <- getSurveyYearList("BTS") # get all Beam trawl surveys
# sapply(yrs, function(x)getSurveyYearQuarterList("BTS",x)) # Get all surveyed quarters by BTS survey years
# # getDatrasDataOverview("BTS") # start year 1985, Q1+Q3+Q4
# # getDatrasDataOverview("SNS") # start year = 1985, Q3+Q4
# # getDatrasDataOverview("BTS-VIII") # start year = 2011, Q4
# # getDatrasDataOverview("DYFS") # start year = 1985, Q3+Q4

# IMPORTING DATRAS haul (HH) DATA ####
# HH: Haul meta data
HH_BTS = getDATRAS(record = "HH",
                   survey = "BTS",
                   years = min(study_period):max(study_period), quarters = c(1L:4L))
setDT(HH_BTS)
HH_DYFS = getDATRAS(record = "HH",
                    survey = c("DYFS"),
                    years = min(study_period):max(study_period), quarters = c(1L:4L))
setDT(HH_DYFS)
HH_BTS8 = getDATRAS(record = "HH",
                    survey = c("BTS-VIII"),
                    years = min(study_period):max(study_period), quarters = c(1L:4L))
setDT(HH_BTS8)
HH_SNS = getDATRAS(record = "HH",
                   survey = c("SNS"),
                   years = min(study_period):max(study_period), quarters = c(1L:4L))
setDT(HH_SNS)

HHflats <- rbind(HH_BTS, HH_BTS8, HH_DYFS, HH_SNS)
if(tosave)rm(HH_BTS, HH_BTS8, HH_DYFS, HH_SNS)

HHflats[Year %in% c(1985:1989),c("ShootLat","HaulLat","ShootLong","HaulLong")]
# replace -9 with NAs => HaulLat and HaulLong are not available for 1985-1989, only the ShootLat and ShootLong
for (jj in 1:ncol(HHflats))set(HHflats, i = which(HHflats[[jj]]==(-9)), j = jj, v = NA)


#### CVM: the following code sets the geometry inconsistently to either the mean position of the haul or the shoot position.         ####
####      The geometry in dataset sf_HHflats is adjusted in 7_append_zero_estimates_to_HL.R and sf_HHflats file is overwritten there ####
#
# # Add lat and lon ====
# HHflats[, `:=`(lat = rowMeans(as.matrix(.SD), na.rm=T)), # lat is set to the means of the start and end latitude of the haul
#        .SDcols=c("ShootLat","HaulLat")]
# HHflats[, `:=`(lon = rowMeans(as.matrix(.SD), na.rm=T)), # lon is set to the means of the start and end longitude of the haul
#        .SDcols=c("ShootLong","HaulLong")]
# # HHflats[is.na(lat),Year<1990] # There are several hauls with missing data in older years
# # HHflats[is.na(lon),Year<1990] # There are several hauls with missing data in older years
# # table(is.na(HHflats$lat))
# # table(is.na(HHflats$lon))
# # table(HHflats[is.na(HHflats$lon)]$Year)
# HHflats <- HHflats[!is.na(HHflats$lat) & !is.na(HHflats$lon)]

## Add StatRec from latitude and longitude if NA ====
HHflats[is.na(StatRec),StatRec:=ices.rect2(lon, lat)] # ices.rect2 from mapplots

## Assign haulID and clean HH data ====
## Create a unique Haul-ID: Survey * Year * Quarter * Country * Ship * Gear * StNo * HaulNo
HHflats[,haulID:=paste(Survey, Year, Quarter, Country, Ship, Gear, StNo, HaulNo, sep = "%%")]
HHflats <- unique(HHflats, by="haulID") # Remove hauls with duplicated ID's in HH (e.g. differences in WindSpeed)

## Check gear widths ====
unique(HHflats$Gear)
# "BT8"   "BT4A"  "BT7"   "BT4AI" "BT4P"  "BT4S" # these are all beam trawls for the BTS
BTSsurveybeamtrawls <- c("BT8","BT4A","BT7","BT4AI","BT4P","BT4S")
# BT4A = Beam trawl 4m, aft # so, one beam trawl sorted
# BT4AI = Beam trawl 4m aft - in Irish Sea q3 BTS
# BT4P = Beam trawl 4m, port
# BT4S = Beam trawl 4m, Beam trawl 4m, starboard
unique(HHflats$Gear)[!unique(HHflats$Gear) %in% BTSsurveybeamtrawls]
surveybeamtrawls <- c(BTSsurveybeamtrawls,"BT6","BT3")
table(HHflats$Gear, HHflats$GearEx)
table(HHflats$Gear, HHflats$Survey, HHflats$GearEx) #SB=SingleBeam, DB=DoubleBeam
table(HHflats[Country=="GB"]$Month,HHflats[Country=="GB"]$Gear,HHflats[Country=="GB"]$Ship)
table(HHflats$Gear, HHflats$GearEx)
#         DB   SB
# BT3      1 6929
# BT4A     0  125
# BT4AI    0    0
# BT4P     0    0
# BT4S     0    0
# BT6      0 3578
# BT7      0    0
# BT8   1299 3613
table(HHflats[Gear=="BT8"]$GearEx, HHflats[Gear=="BT8"]$Year)
# 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 2004 2005 2006 2007
# DB    0   59   65    0   82   94   99   99  101   91   88  133  126   53   63   69   74    0    0    0    0    2    0
# SB   62    0    0   93    0    1    0    0    0    0    0    0    1   70   94   91   70  162  163  159  163  151  151
# 2008 2009 2010 2011 2012 2013 2014 2015 2016 2017 2018 2019 2020 2021 2022
# DB    0    0    0    0    0    0    1    0    0    0    0    0    0    0    0
# SB  139  151  111  133  165  149  132  156  161  141  157  146  148  144  149
table(HHflats$Year, HHflats$Area_27, HHflats$Quarter)
table(HHflats$Area_27, HHflats$Gear, HHflats$Quarter)

HHflats[,singlegearwidth_m := numeric(length=.N)]
HHflats[Gear=="BT8",singlegearwidth_m :=8]
HHflats[Gear %in% c("BT4A","BT4AI","BT4S"),singlegearwidth_m :=4]
HHflats[Gear=="BT7",singlegearwidth_m :=7]
HHflats[Gear=="BT6",singlegearwidth_m :=6]
HHflats[Gear=="BT3",singlegearwidth_m :=3]
HHflats[singlegearwidth_m==0,singlegearwidth_m:=NA]

table(HHflats$Survey, HHflats$Ship, HHflats$GearEx, useNA="always")

HHflats[,gearwidth_m :=singlegearwidth_m]
HHflats[GearEx=="DB",gearwidth_m :=2*singlegearwidth_m] # Gear width is always single gear width unless otherwise specified

# A single beam has been used unless otherwise specified
# BT4a, BT4AI, BT4P and BT4S are by definition one beam trawl sorted
# BT7 and BT8 are suggested to have one beam trawl sorted (ICES BTS Manual: ICES. 2019. Manual for the Offshore Beam Trawl Surveys, Version 3.4, April 2019, Working Group on Beam Trawl Surveys. 54pp. http://doi.org/10.17895/ices.pub.5353)
# However, BT8 => GearEx=="DB" for some hauls, which means "Double beam, catch from two beam trawls put together" => these hauls use 2*singlegearwidth
# DYFS (BE, NL and GER) and SNS are single beam trawls (Heleen Raat, pers comm)

table(HHflats$GearEx, useNA="always")
table(HHflats$gearwidth_m)

# Add swept areas by haul ====
#	Multiplier = HaulDur/60 if DataType ??? ("C")
table(HHflats$Survey, HHflats$GroundSpeed) # DYFS has a lower fishing speed?
table(HHflats$GroundSpeed, useNA="always")
# BTS: If the groundspeed while fishing was not recorded, it is manually set to 4 knots # WGBEAM ICES Manual 2009 => only France uses trawling speeds = 5 knots, BE+UK+DEU+NL = 4 knots
HHflats[Survey=="BTS" & is.na(GroundSpeed) & Country!= "FR"]$GroundSpeed <- 4
HHflats[Survey=="BTS-VIII" & is.na(GroundSpeed)]$GroundSpeed <- 5
table(HHflats$Survey, HHflats$Country)

# Annex 9 in ICES. 2021. Working Group on Beam Trawl Surveys (WGBEAM). ICES Scientific Reports. 3:46. 89pp. https://doi.org/10.17895/ices.pub.8114
HHflats[Survey=="SNS" & is.na(GroundSpeed) & Country =="NL"]$GroundSpeed <- 3.75 # SNS: If the groundspeed while fishing was not recorded, it is manually set to 3.75 knots for SNS # Annex 9 in ICES. 2021. Working Group on Beam Trawl Surveys (WGBEAM). ICES Scientific Reports. 3:46. 89pp. https://doi.org/10.17895/ices.pub.8114
HHflats[Survey=="DYFS" & is.na(GroundSpeed) & Country =="NL"]$GroundSpeed <- 3 # Annex 9 in ICES. 2021. Working Group on Beam Trawl Surveys (WGBEAM). ICES Scientific Reports. 3:46. 89pp. https://doi.org/10.17895/ices.pub.8114
HHflats[Survey=="DYFS" & is.na(GroundSpeed) & Country =="BE"]$GroundSpeed <- 3.25 
                                                                            # 3.5 knots in Annex 9 in ICES. 2021. Working Group on Beam Trawl Surveys (WGBEAM). ICES Scientific Reports. 3:46. 89pp. https://doi.org/10.17895/ices.pub.8114
                                                                            # 3 knots in Annex 5.2 in ICES. 2021. Working Group on Beam Trawl Surveys (WGBEAM). ICES Scientific Reports. 3:46. 89pp. https://doi.org/10.17895/ices.pub.8114
HHflats[Survey=="DYFS" & is.na(GroundSpeed) & Country =="DE"]$GroundSpeed <- 3 # 3 knots in Annex 5.2 in ICES. 2011. Report of the Working Group on Beam Trawl Surveys (WGBEAM), 7-10 June 2011, Hamburg, Germany. ICES CM 2011/SSGESST:14. 225 pp.
table(HHflats$GroundSpeed, useNA="always")

HHflats$GroundSpeed <- ifelse(is.na(HHflats$GroundSpeed), 4, HHflats$GroundSpeed)

table(HHflats$GroundSpeed, useNA="always")

# The surface area that was swept by each haul is calculated
HHflats <- HHflats %>% mutate(sweptarea_sqkm = (gearwidth_m *((GroundSpeed/(HaulDur/60))*1852))/1000000)
                                                  # Convert HaulDur from  minutes to hours (HaulDur/60)
                                                  # Convert GroundSpeed from knot to m/h (*1852)
                                                  # convert m? to km? by /1000000


# convert to sf back to degrees as coordinate system
sf_HHflats = st_as_sf(HHflats, coords = c("lon","lat"),crs = st_crs(4326))
st_crs(sf_HHflats, parameters = TRUE)$units_gdal # degrees

# Convert crs to planar systems using the British national grid
ICESareas <- st_transform(ICESareas, crs = 27700) # we need planar coordinates for st_join
st_crs(ICESareas, parameters = TRUE)$units_gdal # metres
sf_HHflats <- st_transform(sf_HHflats, crs = 27700) 
st_crs(sf_HHflats, parameters = TRUE)$units_gdal # metres

# Explore hauls across the years and ICES divisions
sf::sf_use_s2(FALSE)
sf_HHflats <- st_join(sf_HHflats, ICESareas)
table(sf_HHflats$Year, sf_HHflats$Area_27)
sf_HHflats <- sf_HHflats %>%
  filter(!is.na(Area_27)) %>% # all points should be within the ICES divisions
  filter(Area_27 %in% study_area) # all points should be within the STUDY AREA
table(sf_HHflats$Year, sf_HHflats$Area_27)

# convert to sf back to degrees as coordinate system
sf_HHflats <- st_transform(sf_HHflats, crs = 4326)


if(tosave)saveRDS(sf_HHflats,file=(paste0(inputdir,"/sf_HHflats_",
                                               paste0(range(study_period),collapse = "till"),"_in",
                                               gsub("\\.","",paste0(study_area,collapse = "")),
                                               ".RDS")))

# if(tosave){ # slow!!!
#   p <- ggplot()+
#     geom_sf(data=ICESareas_study)+
#     geom_sf(data=sf_HHflats, aes(col=Year),alpha=0.15, size=0.75)+
#     facet_grid(Quarter~Survey)+
#     scale_colour_gradientn(colours = brewer.pal(n = 7, name = "Spectral"))
#   print(p)
#   ggsave(paste0(outputdir,"/FlatfishSurveys_HaulsbySurveyQuarter.png"),width = 15,height = 10,dpi = 300)
# }
# if(tosave){ # slow!!!
#   p <- ggplot()+
#     geom_sf(data=ICESareas_study)+
#     geom_sf(data=sf_HHflats, aes(col=as.factor(Quarter)),alpha=0.15, size=0.75)+
#     facet_wrap(Survey~Year,nrow = 6)+
#     scale_color_manual(values = brewer.pal(n = 4, name = "Spectral"))
#   print(p)
#   ggsave(paste0(outputdir,"/FlatfishSurveys_HaulsbyYearQuarter.png"),width = 20,height = 10,dpi = 300)
# }
# rm(p, BTSsurveybeamtrawls, sf_HHflats, yrs, ICESareas_study)

# IMPORTING DATRAS Length (HL) and Age (CA) DATA ####
# Extracting Species Length (HL) and Age (CA) data from the DATRAS BTS datasets (two loops)
# HL: Species length-based information
# CA: Species age-based information
# ICES FAQ: It is possible to relate records by the first 11 fields. 
# So, HL- and CA-records would belong to HHrecord with the same quarter, country, ship, gear, station number, haul number, and year.
#
# the loop will take time depending on the study_area, the study_period and the number of species for which aphiaIDs are requested online (Worms)
#
sf_HHflats <- readRDS(file=(paste0(inputdirflats,"/sf_HHflats_",paste0(c(1985,2022),collapse = "till"),"_in",
                                   gsub("\\.","",paste0(study_area,collapse = "")),
                                   ".RDS")))
HHflats <- cbind(as.data.table(sf_HHflats),st_coordinates(sf_HHflats)) # To get back from sf to dt
HHflats[,geometry:=NULL]
setnames(HHflats,old=c("X","Y"),new=c("lon","lat"))
HHflats <- HHflats[Year %in% study_period]
sort(unique(HHflats$Year))

for(recordtypeLoop in c("HL","CA")){
  if(testrun){
    recordtypeLoop <- "HL"
    study_period <- c(1985L, 2006L, 2015L)
    study_species <- 127160 # vector with valid Aphia IDs
  }
  progress(which(c("HL","CA")==recordtypeLoop), progress.bar=T)
  Sys.sleep(0.01)
  if(recordtypeLoop=="CA")cat("Last run, almost done?!\n")
  ## Import HL or CA data through the R icesDATRAS pacakge, and merge with HH data ====
  rec_BTS = getDATRAS_mymodification(record = recordtypeLoop, # the getDATRAS function is used, but modified to track progress
                     survey = "BTS",
                     years = as.integer(study_period), quarters = c(1L:4L))
  # readDatras gave errors for some URLs, which is why Globa options in RStudio were manually modified following https://stackoverflow.com/questions/22721819/download-file-fails-in-rstudio
  # Go to Tools > Global Options > Packages, and unselect "Use Internet Explorer library/proxy for HTTP"
  setDT(rec_BTS)
  print(recordtypeLoop)
  if(tosave)saveRDS(rec_BTS,paste0(inputdirflats,"/BTS_",recordtypeLoop,paste0(range(study_period),collapse = "till"),"_in",gsub("\\.","",paste0(study_area,collapse = "")),".RDS"))
  rec_BTS8 = getDATRAS_mymodification(record = recordtypeLoop,
                      survey = "BTS-VIII",
                      years = as.integer(study_period), quarters = c(1L:4L))
  setDT(rec_BTS8)
  print(recordtypeLoop)
  if(tosave)saveRDS(rec_BTS8,paste0(inputdirflats,"/BTS8_",recordtypeLoop,paste0(range(study_period),collapse = "till"),"_in",gsub("\\.","",paste0(study_area,collapse = "")),".RDS"))
  rec_SNS = getDATRAS_mymodification(record = recordtypeLoop,
                      survey = "SNS",
                      years = as.integer(study_period), quarters = c(1L:4L))
  setDT(rec_SNS)
  print(recordtypeLoop)
  if(tosave)saveRDS(rec_SNS,paste0(inputdirflats,"/SNS_",recordtypeLoop,paste0(range(study_period),collapse = "till"),"_in",gsub("\\.","",paste0(study_area,collapse = "")),".RDS"))
  rec_DYFS = getDATRAS_mymodification(record = recordtypeLoop,
                      survey = "DYFS",
                      years = as.integer(study_period), quarters = c(1L:4L))
  setDT(rec_DYFS)
  print(recordtypeLoop)
  if(tosave)saveRDS(rec_DYFS,paste0(inputdirflats,"/DYFS_",recordtypeLoop,paste0(range(study_period),collapse = "till"),"_in",gsub("\\.","",paste0(study_area,collapse = "")),".RDS"))
  
  rec_flats <- rbind(rec_BTS, rec_BTS8, rec_DYFS, rec_SNS)
  nrow(rec_flats) # 119916 rows for HL 2015 (all species)
  rec_flats <- rec_flats[Valid_Aphia %in% study_species] # Limit to the study species !!
  nrow(rec_flats) # 4768 rows for HL 2015 (study_species=="Solea solea")
  
  # if(tosave)rm(rec_BTS, rec_BTS8, rec_DYFS, rec_SNS)
  
  ## Assign haulID and clean HL or CA data ====
  rec_flats[,haulID:=paste(Survey, Year, Quarter, Country, Ship, Gear, StNo, HaulNo, sep = "%%")]
  HHflats <- HHflats[haulID %in% rec_flats$haulID] # Remove hauls if haul-ID is not in HL or CA (no length or number or age data available for the haul)
  
  ## Merge HH and HL or CA data ====
  rmcols <- c("RecordType","SweepLngt","DoorType", "DateofCalculation") # exclude SweepLngt and DoorType, because they are not relevant for Beam trawls
  HHflats[,(rmcols):=NULL]
  rec_flats[,(rmcols):=NULL]

  HHHLCAkeys <- c("haulID","Survey","Quarter","Country","Ship","Year","Gear","GearEx","StNo","HaulNo")
  
  setkeyv(HHflats, HHHLCAkeys)
  setkeyv(rec_flats, HHHLCAkeys)
  HHflats[,..HHHLCAkeys] # all BTS, i.e. not only in the study_area  # 29169 records in 1993
  rec_flats[,..HHHLCAkeys]
  
  dt_HHrec_flats_study <- rec_flats[HHflats] # inner join of all hauls
  dt_HHrec_flats_study[,..HHHLCAkeys] # BTS only in the study_area # e.g. 12080 records in 1993 in 7afg
  
  
  ## Only keep hauls that fulfill the criteria: ====
  #	HaulVal == ("V", "A")
  table(dt_HHrec_flats_study$HaulVal)
  dt_HHrec_flats_study <- dt_HHrec_flats_study[HaulVal %in% c("V")] # => # V to only include Valid hauls
  #	DayNight == ????
  table(dt_HHrec_flats_study$DayNight) # 3541 Night hauls which are excluded
  dt_HHrec_flats_study <- dt_HHrec_flats_study[DayNight %in% c("D")] # => All HL hauls are Daylight hauls dt_HHrec_flats_study in 7afg in 1993
  #	StdSpecRecCode == 1
  table(dt_HHrec_flats_study$StdSpecRecCode)
  dt_HHrec_flats_study <- dt_HHrec_flats_study[StdSpecRecCode==1]
  #	SpecVal ??? (1,4,7,10)
  if(recordtypeLoop=="HL"){
    table(dt_HHrec_flats_study$SpecVal)
    dt_HHrec_flats_study <- dt_HHrec_flats_study[SpecVal %in% c(1,4,7,10)]
  }
  #	lat != NA & lon != NA # We need haul locations
  dt_HHrec_flats_study <- dt_HHrec_flats_study[!is.na(lat)|!is.na(lon)]
  
  ##  Replace all -9 values with NAs ====
  # dt_HHrec_flats_study <- dt_HHrec_flats_study %>% dplyr::na_if(-9) # na_if should be applied to vectors but not entire tables/dataframes
  for (jj in 1:ncol(dt_HHrec_flats_study))set(dt_HHrec_flats_study, i = which(dt_HHrec_flats_study[[jj]]==(-9)), j = jj, v = NA)

  ## Assign species to AphiaID ====
  # The species in the hauls are coded using an aphia code; a species list is created from the aphia codes and the species list is tidied up
  species_list = worms::wormsbyid(unique(dt_HHrec_flats_study$Valid_Aphia))
  species_list = species_list %>% 
    dplyr::select(Valid_Aphia = AphiaID, scientificname, phylum, class, order, family, genus)
  # The species name from "worms" are added to the dataset
  setDT(species_list)
  species_list <- unique(species_list[!is.na(Valid_Aphia)])
  # if(tosave)saveRDS(species_list,paste0(inputdirflats,"/Worms_SpeciesList_surveyflats_in",
  #                                       paste0(range(study_period),collapse = "till"),"_in",
  #                                       gsub("\\.","",paste0(study_area,collapse = "")),
  #                                       ".RDS"))
  setkey(species_list, Valid_Aphia)
  setkey(dt_HHrec_flats_study, Valid_Aphia)
  dt_HHrec_flats_study <- species_list[dt_HHrec_flats_study]
  any(is.na(dt_HHrec_flats_study$Valid_Aphia)) # any NAs?
  any(dt_HHrec_flats_study$SpecVal==0) # any invalid info?
  dt_HHrec_flats_study <- dt_HHrec_flats_study %>% 
    filter(!is.na(Valid_Aphia)) %>% # eliminate NAs
    filter(Valid_Aphia!=0) # eliminate zeros

  BTSsurveybeamtrawls <- c("BT8","BT4A","BT7","BT4AI","BT4P","BT4S")
  unique(dt_HHrec_flats_study$Gear)[!unique(dt_HHrec_flats_study$Gear) %in% BTSsurveybeamtrawls]
  surveybeamtrawls <- c(BTSsurveybeamtrawls,"BT6","BT3")
  
  dt_HHrec_flats_study$LngtClass_mm <- NA
  if("LngtClass" %in% names(dt_HHrec_flats_study) & any(dt_HHrec_flats_study$Gear %in% surveybeamtrawls)){
    dt_HHrec_flats_study[,LngtClass_mm := list(ifelse(LngtCode==".",LngtClass, # 1 mm length class, reporting units: mm
                                                    ifelse(LngtCode=="0",LngtClass, #  0.5 cm length class, reporting units: mm
                                                           ifelse(LngtCode=="1",LngtClass*10, # 1 cm length class, reporting units: cm
                                                                  NA))))]
    }
  
  
  if(recordtypeLoop=="HL"){
    dt_HHrec_flats_study <- dt_HHrec_flats_study %>% 
      filter(SpecVal!=0) # remove invalid information
    ## Estimate numbers (Count) ====
    #	If HLNoAgeLngt == NA, use TotalNo and set SubFactor = 1
    dt_HHrec_flats_study$HLNoAtLngt_old = dt_HHrec_flats_study$HLNoAtLngt
    dt_HHrec_flats_study$SubFactor_old = dt_HHrec_flats_study$SubFactor
    changecols <- c("HLNoAtLngt","SubFactor")
    dt_HHrec_flats_study[is.na(HLNoAtLngt),SubFactor:=1]
    dt_HHrec_flats_study[is.na(HLNoAtLngt),HLNoAtLngt:=TotalNo]
    
    selcols <- c("HLNoAtLngt_old","HLNoAtLngt","TotalNo")
    dt_HHrec_flats_study[,(selcols):= lapply(.SD, as.numeric), .SDcols = selcols]
    
    # Multiplier = SubFactor if DataType ??? ("S","R") ====
    #! NOT DONE
    # Data Type can be reported as 
    # R (raw data reported, sorted catch might be sub-sampled), 
    # S (bulk unsorted catch was sub-sampled)
    # C (catch reported as CPUE).
    # d.	Count = HLNoAtLngt * Multiplier
    table(round(dt_HHrec_flats_study$SubFactor,1))
    dt_HHrec_flats_study[,HLNoAtLngt:=HLNoAtLngt*SubFactor]
  } # end loop for HL data
  
  study_spName <- paste0(gsub(" ","_",species_list[Valid_Aphia %in% study_species]$scientificname),collapse="%%")
  if(tosave)saveRDS(dt_HHrec_flats_study,file=paste0(inputdirflats,"/",study_spName,recordtypeLoop,"_",
                                                     paste0(unique(dt_HHrec_flats_study$Survey),collapse="+"),
                                                     "_in",paste0(range(study_period),collapse = "till"),"_in",
                                                     gsub("\\.","",paste0(study_area,collapse = "")),
                                                     ".RDS"))

} # End of loop where HL or CA are merged to HH data


