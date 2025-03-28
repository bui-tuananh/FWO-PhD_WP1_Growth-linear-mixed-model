# 1. SETUP ----
library(tidyverse)
library(sf)
library(tmap)

# 2. LOAD DATA ----

## ICES statistical rectangles and division ----
## stat_rec
dir_stat_rec <- "./data/ices/ICES_statistical-rectangle"
ices_rec <- read_sf(file.path(dir_stat_rec, "ICES_Statistical_Rectangles_Eco.shp"))

## ices_div 
# note: area is the whole 4abc, but keep naming 4bc to be consistent with other processed data
dir_admin <- "./data/admin"
ices_div <- read_sf(file.path(dir_admin, "ices_areas_sub_group_4abc_4326_new.gpkg")) #ices_areas_sub_group_4abc_4326_new.gpkg #ices_areas_sub_group_4326_new.gpkg
ices_div <- ices_div %>% 
  filter(Area_27 %in% c("4abc", "7a", "8ab"))

## DATRAS data ----
dir_datras <- "./data/ices"
hh <- read_rds(file.path(dir_datras, "sf_HHflats_1985till2022_in3a204a4b4c7a7d7e7f7g7h7j28a8b.RDS"))
hl <- read_rds(file.path(dir_datras, "sf_Solea_soleaHL_withAbs_BTS+BTS-VIII+DYFS+SNS_in1985till2021_in3a204a4b4c7a7d7e7f7g7h7j28a8b.RDS"))

# 3. PROCESS DATA ----
## 3.1. ICES stat_rec within ices_div ----
## check crs
st_crs(ices_rec) == st_crs(ices_div) # TRUE

## intersect ices_rec ices_div
ices_rec <- st_intersection(ices_rec, ices_div)

## ices_rec at borders between nearby ices_div can be duplicated -> check
# view to see ices
tmap_mode("view")

tm_shape(ices_rec %>% select(ICESNAME)) +
  tm_polygons() +
  tm_shape(ices_div) +
  tm_borders(lwd = 2)

# check duplication
ices_rec_dup <- as_tibble(ices_rec) %>% 
  group_by(ICESNAME) %>%
  summarize(n = n()) %>%
  filter(n > 1) %>% 
  inner_join(ices_rec %>% select(ICESNAME, Area_27))

write_rds(ices_rec, file.path(dir_stat_rec, "ices_rec_4abc7a8ab.rds")) #error in writing sf

# 3.2. sol distribution area (DATRAS survey data) ----

# survey locations hh
hh <- hh %>% 
  select(StatRec) %>%
  unique()

tm_shape(ices_rec %>% select(ICESNAME)) +
  tm_polygons() +
  tm_shape(hh) +
  tm_dots()

# survey with catch of sol - hl
# note: only quarter 3,4, as quarter 1 only a few points in 7a (mainly 7efg)

# survey locations
hl_sum <- hl %>% 
  filter(!is.na(HLNoAtLngt),
         Quarter %in% c(3,4),
         Area_27 %in% c("4.a", "4.b", "4.c", "7.a", "8.a", "8.b")) %>%
  select(StatRec) %>%
  unique()

# survey ices rec
hl_loc <- as_tibble(hl) %>% 
  filter(!is.na(HLNoAtLngt),
         Quarter %in% c(3,4)) %>%
  select(StatRec) %>%
  unique()

# join with ices_rec 
hl_loc <- hl_loc %>% 
  inner_join(ices_rec, by = join_by(StatRec == ICESNAME)) %>%
  st_as_sf() %>%
  select(StatRec, Area_27) 

# dissolve and calculate area
hl_loc <- hl_loc %>%
  group_by(Area_27) %>%
  summarize(Area_27 = first(Area_27))

hl_loc <- hl_loc %>%
  mutate(area_km2 = as.numeric(st_area(hl_loc)/(1000^2)))

# visualisation
tm_shape(hl_loc) +
  tm_borders() +
  tm_shape(hl_sum) +
  tm_dots(alpha = 0.1)

# save file
write_sf(hl_loc, file.path(dir_datras, "hl_loc_4abc7a8ab.gpkg"))
