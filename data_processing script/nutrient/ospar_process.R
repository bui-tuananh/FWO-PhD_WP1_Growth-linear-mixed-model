# source: ospar_process_Charlotte


## setup ----
library(sf)
library(mapview)
library(data.table)
library(nngeo)
library(climwin)
library(statmod)
library(lme4)
library(tweedie)
library(nloptr)
library(units)
library(sf)
library(tidyverse)
library(tmap)

tmap_mode("view") #set tmap view

#source("C:/Users/cvanmoorleghem/OneDrive - ILVO/Development/RStudio/SEAwise T5.4/seawise_t5-4_juvsol7afg/R/Utilities/slidingwin2.R") # function from package climwin adapted to provide tweedie AIC instead of AICc
#source("C:/Users/cvanmoorleghem/OneDrive - ILVO/Development/RStudio/SEAwise T5.4/seawise_t5-4_juvsol7afg/R/Utilities/basewin3.R") # function from package climwin adapted to provide tweedie AIC instead of AICc
# slidingwin2 calculates the AIC for tweedie glm's with the tweedie:AICtweedie function

## note ----
# main rivers according to copernicus website (https://data.marine.copernicus.eu/viewer/expert?view=viewer&crs=epsg%3A4326&t=757605600000&z=-300&center=-3.2067157561987027%2C53.090325139036935&zoom=14.300000000000004&layers=W3siaWQiOiJjMCIsImxheWVySWQiOiJJTlNJVFVfR0xPX1BIWV9UU19ESVNDUkVURV9NWV8wMTNfMDAxL2NtZW1zX29icy1pbnNfZ2xvX3BoeS10ZW1wLXNhbF9teV9jb3JhX2lyci9fX0RFRkFVTFRfXyIsInpJbmRleCI6MTB9LHsiaWQiOiJjMSIsImxheWVySWQiOiJJTlNJVFVfR0xPX1BIWV9UU19ESVNDUkVURV9NWV8wMTNfMDAxL2NtZW1zX29icy1pbnNfZ2xvX3BoeS10ZW1wLXNhbF9teV9lYXN5Y29yYV9pcnIvX19ERUZBVUxUX18iLCJ6SW5kZXgiOjIwfSx7ImlkIjoiYzIiLCJsYXllcklkIjoiTldTSEVMRl9NVUxUSVlFQVJfUEhZXzAwNF8wMDkvY21lbXNfbW9kX253c19waHktc3NzX215XzdrbS0yRF9QVDFILWkvc28iLCJ6SW5kZXgiOjMwLCJsb2dTY2FsZSI6ZmFsc2V9XQ%3D%3D) de belangrijkste geselecteerd:
#   Seine v - 7d
#   Severn v - 7fg
#   Trent
#   Thames v - 4bc
#   Scheldt v - 4bc
#   Meuse v - 4bc
#   Rhine v - 4bc
#   Weser v - 4bc
#   Elbe v - 4bc
#   Ems v - 4bc
#   Skagerrak/Kattegat
#   Haven, Welland, Nene, Great Ouse /v/v/v
#   Mersey/Dee v/v
#   Kent/Leven v/v
#   Esk/Eden v/v
#   Humber v - 4bc
#   IJselmeren zijn er ook nog, maar grootste invloed komt toch precies van Rijn, Maas & Schelde, dus niet gebruiken?

# rivers that were mainly important for the Bay of Biscay, which is for now excluded from the dataset
#   Garonne v
#   Dordogne v
#     These two together form the Gironde
#   Loire v

# misschien valt er nog literatuur te vinden zoals Lacroix et al 2004, maar dan voor andere maritime regio's?

# v : in Van Leeuwen et al 2023 dataset

## get the location of all rivers ----
riverdataList <- list.files(path = "./Env/Eutrophication_OSPAR ICG-EMO/CurrentState1/",
                            recursive = TRUE,
                            pattern = "\\_1940_2022.dat$",
                            full.names = TRUE)

river_data_map <- data.frame()
for(r in riverdataList[1:length(riverdataList)]) { print(paste0("processing data for the river ",gsub(".*//(.+)_1940_2022.dat.*", "\\1",r)))
  x <- gsub(".*//(.+)_1940_2022.dat.*", "\\1",r)
  y <- sub(") ","",sub(".*(N,E)", "",readLines(r, n=1)))
  z <- read.table(r, header=TRUE, skip=6)
  z <- z[2:nrow(z),]
  #z <- z[which(z$year >= 1980),]
  z <- z %>% dplyr::filter(year == "2022", m == "1", d == "1") #get 1 value only to get coord
  z$reg_river <- x
  z$coord <- y
  river_data_map <- rbind(river_data_map,z)
}

river_data_map$Region <- sub("/.*","",river_data_map$reg_river)
river_data_map$River <- sub(".*/","",river_data_map$reg_river)
river_data_map$long <- as.numeric(sub(".*, ","",river_data_map$coord))
river_data_map$lat <- as.numeric(sub(", .*","",river_data_map$coord))
river_data_map <- st_as_sf(river_data_map, coords = c("long","lat"),crs = st_crs(4326))
river_data_map$long <- as.numeric(sub(".*, ","",river_data_map$coord))
river_data_map$lat <- as.numeric(sub(", .*","",river_data_map$coord))
river_data_map$coord <- NULL
river_data_map$reg_river <- NULL
# river_data <- river_data[,c("Region","River","lat","long")]

# keep only region, river, lat, long, and geometry
river_data_map <- river_data_map %>% select(Region:lat)
write_sf(river_data_map, file.path("./Env/Eutrophication_OSPAR ICG-EMO", "ospar_map.gpkg"))

# plot map
dir_gis <- "D:/OneDrive - UGent/data/Admin"
ices_div <- read_sf(file.path(dir_gis, "ices_areas_sub_group_4326_new.gpkg"))

river_data_map <- read_sf(file.path("./Env/Eutrophication_OSPAR ICG-EMO", "ospar_map.gpkg"))

ggplot() +
  geom_sf(data = ices_div) +
  geom_sf(data = river_data_map, aes(color = Region)) 

## get the data of selected rivers ----
# remove region Baltic and region PT
riverdataList_df <- tibble(id = 1:length(riverdataList),
                           riverdataList = riverdataList,
                           reg_river = gsub(".*//(.+)_1940_2022.dat.*", "\\1",riverdataList),
                           region = sub("/.*","", reg_river),
                           river = sub(".*/","", reg_river)) %>%
  filter(!region %in% c("Baltic", "PT"))
riverdataList_sub <- riverdataList_df$riverdataList

river_data <- data.frame()
for(r in riverdataList_sub[1:length(riverdataList_sub)]) { print(paste0("processing data for the river ",gsub(".*//(.+)_1940_2022.dat.*", "\\1",r)))
  x <- gsub(".*//(.+)_1940_2022.dat.*", "\\1",r)
  y <- sub(") ","",sub(".*(N,E)", "",readLines(r, n=1)))
  z <- read.table(r, header=TRUE, skip=6)
  z <- z[2:nrow(z),]
  #z <- z[which(z$year >= 1980),]
  z$reg_river <- x
  z$coord <- y
  river_data <- rbind(river_data,z)
}

river_data$Region <- sub("/.*","",river_data$reg_river)
river_data$River <- sub(".*/","",river_data$reg_river)
river_data$long <- as.numeric(sub(".*, ","",river_data$coord))
river_data$lat <- as.numeric(sub(", .*","",river_data$coord))
river_data <- st_as_sf(river_data, coords = c("long","lat"),crs = st_crs(4326))
river_data$long <- as.numeric(sub(".*, ","",river_data$coord))
river_data$lat <- as.numeric(sub(", .*","",river_data$coord))
river_data$coord <- NULL
river_data$reg_river <- NULL
# river_data <- river_data[,c("Region","River","lat","long")]

# set character vectors as integer and numeric where appropriate
river_data$year <- as.integer(river_data$year)
river_data$m <- as.integer(river_data$m)
river_data$d <- as.integer(river_data$d)
river_data$Q <- as.numeric(river_data$Q)
river_data$TN <- as.numeric(river_data$TN)
river_data$NO3 <- as.numeric(river_data$NO3)
river_data$NH4 <- as.numeric(river_data$NH4)
river_data$DIN <- as.numeric(river_data$DIN)
river_data$TP <- as.numeric(river_data$TP)
river_data$PO4 <- as.numeric(river_data$PO4)
river_data$Silicate <- as.numeric(river_data$Silicate)
river_data$TALK <- as.numeric(river_data$TALK)
river_data$DIC <- as.numeric(river_data$DIC)
river_data$DOC <- as.numeric(river_data$DOC)
river_data$SPM <- as.numeric(river_data$SPM)
river_data$NO2 <- as.numeric(river_data$NO2)
river_data$Fe <- as.numeric(river_data$Fe)
river_data$POC <- as.numeric(river_data$POC)
river_data$TOC <- as.numeric(river_data$TOC)

write_rds(river_data, file.path("./Env/Eutrophication_OSPAR ICG-EMO", "ospar_subset_1940-2022.rds"))

## get data in studied areas - 4bc, 7a, 8ab ----
#### data
## ospar
ospar <- as.data.frame(read_rds(file.path("./Env/Eutrophication_OSPAR ICG-EMO", "ospar_subset_1940-2022.rds")))

## ospar_map and ices_div
ospar_map <- read_sf(file.path("./Env/Eutrophication_OSPAR ICG-EMO", "ospar_map.gpkg"))

dir_gis <- "./Admin"
ices_div <- read_sf(file.path(dir_gis, "ices_area_sub_4abc_4326.gpkg"))

tm_shape(ices_div) +
  tm_borders() +
  tm_shape(ospar_map) +
  tm_dots()

### buffer ospar_map then intersect with ices_div ----
# to get list of rivers going to the study area
ospar_map_buffer <- st_buffer(ospar_map, 5000) #unit m
ospar_map_buffer_ices <- st_intersection(ospar_map_buffer, ices_div)
ospar_map_buffer_ices <- ospar_map_buffer_ices %>%
  mutate(id = seq(1, nrow(ospar_map_buffer_ices)))

## check if any duplicate - same rivers having different ices_div (due to too big buffer area)
ospar_dup <- ospar_map_buffer_ices %>% filter(duplicated(River) == T)
ospar_dup <- ospar_map_buffer_ices %>% filter(River %in% ospar_dup$River)
# no duplication

tm_shape(ospar_dup) +
  tm_borders() +
  tm_shape(ospar_map %>% filter(River %in% c("LUD", "UGIE"))) +
  tm_dots() +
  tm_shape(ices_div) +
  tm_borders()
 
## duplication
## UK UGIE: 4a 4b -> remove 4b (id 80)
## UK LUD: 4b 4c -> remove 4c (id 24)
ospar_map_buffer_ices <- ospar_map_buffer_ices %>% filter(!id %in% c(24, 80), )

## remove River Mulkear, Shannon (inland) due to error in ices_div file
ospar_map_buffer_ices <- ospar_map_buffer_ices %>% filter(!River %in% c("Mulkear", "Shannon") )

## remove River Ry (id 49) - seems like running to Kattegat (not North Sea)
ospar_map_buffer_ices <- ospar_map_buffer_ices %>% filter(!id %in% c(49) )

## plot
tm_shape(ospar_map_buffer_ices) +
  tm_polygons("Area_27") +
  tm_dots() +
  tm_shape(ices_div) +
  tm_borders()

## save ospar_map_buffer_ices
st_write(ospar_map_buffer_ices, file.path("./Env/Eutrophication_OSPAR ICG-EMO", "ospar_map_buffer5km_ices.gpkg"))

### extract ospar data in studied area ----
ospar_map_buffer_ices <- read_sf(file.path("./Env/Eutrophication_OSPAR ICG-EMO", "ospar_map_buffer5km_ices.gpkg"))

ospar_map_buffer_ices_sub <- as.data.frame(ospar_map_buffer_ices) %>%  
  mutate(IcesArea = Area_27) %>%
  select(River, IcesArea, Area_Full)

ospar_ices <- ospar %>% 
  filter(River %in% ospar_map_buffer_ices$River) %>%
  left_join(ospar_map_buffer_ices_sub)

## save ospar_ices
write_rds(ospar_ices, file.path("./Env/Eutrophication_OSPAR ICG-EMO", "ospar_subset_1940-2022_ices.rds"))

### summarize ospar data by ices_area and by year ----
ospar_ices <- read_rds(file.path("./Env/Eutrophication_OSPAR ICG-EMO", "ospar_subset_1940-2022_ices.rds"))

#### overview of river contribution ----
ospar_ices_river <- ospar_ices %>% 
  mutate(TN = if_else(TN < 0, NA, TN),
         TP = if_else(TP < 0, NA, TP)) %>%
  group_by(IcesArea, Region, River, year) %>%
  summarize(TN = sum(TN),
            TP = sum(TP)) 

river_cont <- ospar_ices_river %>% 
  mutate(is_TN = if_else(is.na(TN) == F, "yes", "no"),
         is_TP = if_else(is.na(TP) == F, "yes", "no")) %>%
  group_by(IcesArea) %>%
  mutate(sum_TN_per_area = sum(TN, na.rm = T),
         sum_TP_per_area = sum(TP, na.rm = T)) %>%
  group_by(Region, River) %>%
  mutate(sum_TN_per_river = sum(TN, na.rm = T),
         sum_TP_per_river = sum(TP, na.rm = T),
         prop_TN = sum_TN_per_river/sum_TN_per_area*100,
         prop_TP = sum_TP_per_river/sum_TP_per_area*100) %>%
  select(IcesArea, Region, River, is_TN, is_TP, 
         sum_TN_per_river, sum_TP_per_river, sum_TN_per_area, sum_TP_per_area, 
         prop_TN, prop_TP) %>%
  unique()

river_cont <- river_cont %>% 
  arrange(IcesArea, desc(prop_TN)) %>% 
  group_by(IcesArea) %>%
  mutate(cumsum_TN = cumsum(prop_TN),
         cumsum_TP = cumsum(prop_TP)) 


## save file
write_rds(river_cont, file.path("./Env/Eutrophication_OSPAR ICG-EMO", "ospar_subset_1940-2022_ices_river-contribution.rds"))

#### check issue of climatology data of major rivers in each IcesArea -----
# major rivers: account more than 80% of total TN and TP in each IcesArea
# 4bc: "Rhine", "Elbe", "Meuse", "HUMBER", "Weser", "Lake_IJssel_West", "THAMES", "Lake_IJssel_East", "Scheldt", "Ems", 
# "OUSE_AT_KINGS_LYNN", "North_Sea_Canal", 
# 7a: "Avoca", "MERSEY", "RIBBLE", "WEAVER", "DEE_AT_CHESTER", "EDEN_AT_CARLISLE", 
# "TEIFI", "DOUGLAS", "Tolka", "LUNE", "Barrow", "CLWYD", "Nore", "YSTWYTH", "DERWENT"
# 8ab: "Loire", "Garonne", "Vilaine", "Dordogne"

ospar_ices_test_river <- ospar_ices %>% 
  filter(River %in% c("Rhine", "Elbe", "Meuse", "HUMBER", "Weser", "Lake_IJssel_West", "THAMES", "Lake_IJssel_East", "Scheldt", "Ems", 
                      "OUSE_AT_KINGS_LYNN", "North_Sea_Canal",
                      "Avoca", "MERSEY", "RIBBLE", "WEAVER", "DEE_AT_CHESTER", "EDEN_AT_CARLISLE",
                      "TEIFI", "DOUGLAS", "Tolka", "LUNE", "Barrow", "CLWYD", "Nore", "YSTWYTH", "DERWENT",
                      "Loire", "Garonne", "Vilaine", "Dordogne")) %>%
  mutate(TN = if_else(TN < 0, NA, TN),
         TP = if_else(TP < 0, NA, TP)) %>%
  group_by(IcesArea, River, year) %>%
  summarize(TN = sum(TN, na.rm = T),
            TP = sum(TP, na.rm = T)) 

# Check TN_flag and TP_flag by river
river = c("HUMBER") # change river name and check data of each river
test <- ospar_ices %>% 
  filter(River %in% c(river) ) %>%
  filter(year >= 1978, year <= 2017) %>%
  select(River, year, m, d, TN, TN_flag, TP, TP_flag)
unique(test$TN_flag)
unique(test$TP_flag)
View(test %>% filter(TN_flag == "C"))

## 4bc - 1978 - 2017 - most major rivers having D: data flag for T and P
ggplot(data = ospar_ices_test_river %>% filter(IcesArea %in% c("4abc")), 
       aes(x = year, y = TN, color = River)) +
  geom_line() +
  facet_wrap(~ IcesArea) 

ggplot(data = ospar_ices_test_river %>% filter(IcesArea %in% c("4abc")), 
       aes(x = year, y = TP, color = River)) +
  geom_line() +
  facet_wrap(~ IcesArea) 

## 7a - mix of C and other flag -> maybe should not consider
# TN
ggplot(data = ospar_ices_test_river %>% filter(IcesArea %in% c("7a"),
                                               River == "Avoca"), 
       aes(x = year, y = TN, color = River)) +
  geom_line() +
  facet_wrap(~ IcesArea) 

ggplot(data = ospar_ices_test_river %>% filter(IcesArea %in% c("7a"),
                                               River != "Avoca"), 
       aes(x = year, y = TN, color = River)) +
  geom_line() +
  facet_wrap(~ IcesArea) 

# TP
ggplot(data = ospar_ices_test_river %>% filter(IcesArea %in% c("7a"),
                                               River == "DEE_AT_CHESTER"), 
       aes(x = year, y = TP, color = River)) +
  geom_line() +
  facet_wrap(~ IcesArea) 

ggplot(data = ospar_ices_test_river %>% filter(IcesArea %in% c("7a"),
                                               River != "DEE_AT_CHESTER"), 
       aes(x = year, y = TP, color = River)) +
  geom_line() +
  facet_wrap(~ IcesArea) 

# 8ab - mix of C and other flag -> maybe should not consider
ggplot(data = ospar_ices_test_river %>% filter(IcesArea %in% c("8ab")), 
       aes(x = year, y = TN, color = River)) +
  geom_line() +
  facet_wrap(~ IcesArea) 

ggplot(data = ospar_ices_test_river %>% filter(IcesArea %in% c("8ab")), 
       aes(x = year, y = TP, color = River)) +
  geom_line() +
  facet_wrap(~ IcesArea) 

TN_flag <- unique(ospar_ices$TN_flag)
ospar_ices_test <- ospar_ices %>% 
  filter(!grepl("C", TN_flag) )

#### summary ----

## 7a, 8ab data seem to have very high uncertainty as almost all major rivers include TN_flag and TP_flag C (climatology)

## 4bc data seem to have low uncertainty as almost all major rivers include TN_flag and TP_flag D
## except for OUSE_AT_KINGS_LYNN, HUMBER and THAMES
# years included in subset data for 4abc: 1978-2017 (data before or after this period are almost all C flag)
# rivers included in subset data for 4abc:
# "Rhine", "Elbe", "Meuse", "Weser", "Lake_IJssel_West", "Lake_IJssel_East", "Scheldt", "Ems", "North_Sea_Canal"
river_4abc <- c("Rhine", "Elbe", "Meuse", "Weser", "Lake_IJssel_West", "Lake_IJssel_East", "Scheldt", "Ems", "North_Sea_Canal")

# total contribution of select rivers: 73.1% TN, 64% TP out of total TN in 4abc estimated from the ospar data 
river_cont %>% filter(River %in% river_4abc) %>%
  summarize(prop_TN = sum(prop_TN),
            prop_TP = sum(prop_TP))

# keep Region, River
ospar_ices_4abc <- ospar_ices %>% 
  filter(River %in% river_4abc,
         year >= 1978, year <= 2017) %>%
  mutate(TN = if_else(TN < 0, NA, TN),
         TP = if_else(TP < 0, NA, TP)) %>%
  group_by(IcesArea, Region, River, year) %>%
  summarize(TN = sum(TN, na.rm = T),
            TP = sum(TP, na.rm = T)) 

ggplot(data = ospar_ices_4abc %>% filter(River == "Rhine"), 
       aes(x = year, y = TN)) +
  geom_line() 

ggplot(data = ospar_ices_4abc %>% filter(River == "Rhine"), 
       aes(x = year, y = TP)) +
  geom_line() 

## save file
write_rds(ospar_ices_4abc, file.path("./Env/Eutrophication_OSPAR ICG-EMO", "ospar_subset_1978-2017_ices_4abc.rds"))

## plot - selected rivers
ospar_map_4abc <- ospar_map %>% filter(River %in% ospar_ices_4abc$River)

tm_shape(ospar_map) +
  tm_dots() +
  tm_shape(ospar_map_4abc) +
  tm_dots(col = "red") 

#### overall year trend (all data flag) ----
# overall trend
ospar_ices_year <- ospar_ices %>% 
  mutate(TN = if_else(TN < 0, NA, TN),
         TP = if_else(TP < 0, NA, TP)) %>%
  group_by(IcesArea, year) %>%
  summarize(TN = sum(TN, na.rm = T),
            TP = sum(TP, na.rm = T)) 

## save file
write_rds(ospar_ices_year, file.path("./Env/Eutrophication_OSPAR ICG-EMO", "ospar_subset_1940-2022_ices_year.rds"))

## TN
ggplot(data = ospar_ices_year %>% filter(IcesArea %in% c("4bc")), 
       aes(x = year, y = TN, color = IcesArea)) +
  geom_line() +
  xlim(1975, 2022)

ggplot(data = ospar_ices_year %>% filter(IcesArea %in% c("7a")), 
       aes(x = year, y = TN, color = IcesArea)) +
  geom_line() +
  xlim(1975, 2022)

ggplot(data = ospar_ices_year %>% filter(IcesArea %in% c("8ab")), 
       aes(x = year, y = TN, color = IcesArea)) +
  geom_line() +
  xlim(1975, 2022)

## TP
ggplot(data = ospar_ices_year %>% filter(IcesArea %in% c("4bc")), 
       aes(x = year, y = TP, color = IcesArea)) +
  geom_line() +
  xlim(1975, 2022)

ggplot(data = ospar_ices_year %>% filter(IcesArea %in% c("7a")), 
       aes(x = year, y = TP, color = IcesArea)) +
  geom_line() +
  xlim(1975, 2022)

ggplot(data = ospar_ices_year %>% filter(IcesArea %in% c("8ab")), 
       aes(x = year, y = TP, color = IcesArea)) +
  geom_line() +
  xlim(1975, 2022)


