#BEFORE
#      3_report files

#AFTER
#      report


# 1. SETUP ----------------------------------------------------------------------
# PACKAGES
library(tidyverse)  # process data frame
library(lme4)       # fit model
library(AICcmodavg) # calculate AIC
library(effects)    # display effect
library(MuMIn)      # compare models
library(sf)         # process geospatial data
library(reactable)  # interactive table
library(RColorBrewer)# color for visualization
library(broom.mixed)# summarize model results
library(writexl)    # write excel
library(patchwork)  # arrange plot
library(sjPlot)     # summarize table
library(knitr)      # print sjPlot table
library(magrittr)   # print sjPlot table
#library(arm)       # to calculate se - do not load to not override select function
library(stars)
library(tmap)
library(FSA)        # age reading precision
library(lmerTest)   # lmerTest

# DEFAULT THEME
theme_set(theme_classic()) 
theme_update(panel.border = element_rect(colour = "black", fill=NA))

# FUNCTION
#centring function
c. <- function (x) {(x - mean(x))} 
# scale function
s. <- function (x) {(x - mean(x))/sd(x)} 
# function to transform all variable value to +
fun_plus <- function(x) {if_else(is.na(x) == F, "+", "")}

# rsquared function
rsquared.glmm=function(modlist) {
  do.call(rbind,lapply(modlist,function(i) {
    if(inherits(i,"merMod") | class(i)=="merLmerTest") {
      VarF=var(as.vector(fixef(i) %*% t(i@pp$X))) 
      VarRand=colSums(do.call(rbind,lapply(VarCorr(i),function(j) j[1])))
      VarResid=attr(VarCorr(i),"sc")^2
      Rm=VarF/(VarF+VarRand+VarResid)
      Rc=(VarF+VarRand)/(VarF+VarRand+VarResid)
      Rsquared.mat=data.frame(Class=class(i),Marginal=Rm,Conditional=Rc,
                              AIC=AIC(update(i,REML=F))) } 
    else { print("Function requires models of class lm, lme, mer, or    merMod") 
    } } ) ) }

# WORKING DIRECTORY 
dir_output <- "./output_revised"
dir_report <- "./report"

# POPULATION NAME
# create df_pop
df_pop <- tibble(IcesAreaGroup = c("4bc", "7a", "8ab"),
                 pop = c("4bc", "7a", "8ab"),
                 pop.name = factor(c("North Sea", "Irish Sea", "Bay of Biscay"),
                                   levels = c("North Sea", "Irish Sea", "Bay of Biscay")))

# COLOR SCALE
col_scale <- c("ISIMIP" = "#F8766D", "ORAS5" = "#00BFC4", "NEMO-MEDUSA" = "#7CAE00")

# HTML TABLE CSS STYLE
# css style
css_list = list(
  css.firsttablerow =  "font-weight:bold; font-style:normal; border:1px solid;",
  css.depvarhead = "font-style:normal; font-weight:bold; text-align:center; padding-top:.8em; border:1px solid;",
  css.firsttablecol = "text-align:left; border:1px solid;",
  css.centeralign = "text-align:center; border:1px solid;",
  css.randomparts = "font-weight:bold; text-align:left; padding-top:.8em; border:1px solid;",
  css.summary = "border:1px solid;",
  css.summarydata = "text-align:left;"
)

# 2. LOAD DATA ------------------------------------------------------------

## Otolith data ---------------------------------------------------

dir_otl <- "./data"
otl <- read_rds(file.path(dir_otl, "otl_full.rds"))

## check 0 increment
incremet_0 <- otl %>% filter(AnnulusDiameterIncrement.um == 0) 

## check number of increment/year
data_sum <- otl %>% 
  group_by(GrowingYear, IcesAreaGroup) %>% 
  summarize(n_increment = n())
# >= 10 increment/year for each population
data_sum_10 <- otl %>% 
  group_by(GrowingYear, IcesAreaGroup) %>% 
  summarize(n_increment = n()) %>%
  filter(n_increment >= 10) %>%
  mutate(pop.year = paste0(IcesAreaGroup, ":", GrowingYear))

## create ssb index by dividing ssb by ices area
dir_gis <- "./data/admin"
ices_4bc <- as.data.frame(st_read(file.path(dir_gis, "ices_areas_sub_group_4326_new.gpkg")))
ices_4bc <- ices_4bc %>% 
  filter(Area_27 %in% c("4bc", "7a", "8ab")) %>% 
  rename(IcesAreaGroup.area.km2 = Area_km2,
         IcesAreaGroup = Area_27) %>%
  select(IcesAreaGroup, IcesAreaGroup.area.km2) 

## process otl data 
# 1. remove 1 fish 8ab cohort < 1985
# 2. remove 1 wur fish with very small age 1 increment (< 200nm)
# 3. 0 increment, add a small value 0.05 to avoid log issue (smallest increment from ILVO is 0.06)
# 4. add area of ices area

data_otl <- otl %>% 
  filter(FishID != "SOL_52_G1_Z.121_26-03-1990_1602") %>%
  filter(FishID != "sol_fab_0575") %>%
  mutate(AnnulusDiameterIncrement.um = if_else(AnnulusDiameterIncrement.um == 0, 0.05, AnnulusDiameterIncrement.um)) %>%
  left_join(ices_4bc, by = "IcesAreaGroup") 

# 5. rename variables to formulate models easier
data_otl <- data_otl %>% mutate(fishid = FishID,
                                increment = AnnulusDiameterIncrement.um,
                                log.increment = log(AnnulusDiameterIncrement.um),
                                age = Age,
                                log.age = log(Age),
                                aac = AgeAtCapture,
                                log.aac = log(AgeAtCapture), 
                                method = OtolithProcessingMethod,
                                datasource = DataSource,
                                pop = IcesAreaGroup,
                                pop.area.km2 = IcesAreaGroup.area.km2,
                                year = GrowingYear,
                                cohort = Cohort,
                                fyear = factor(GrowingYear),
                                fcohort = factor(Cohort),
                                pop.year = paste0(pop, ":", year),
                                pop.cohort = paste0(pop, ":", cohort)
                                # predictors
                                #f = FishingMortality,
                                #ssb.i = SpawningStockBiomass.1000t*1000/IcesAreaGroup.area.km2, #unit: ton/km2
)

# 6. keep only years with >= 10 increments/year for each pop  
data_otl <- data_otl %>% filter(pop.year %in% data_sum_10$pop.year)

## Temperature data ----
## load and process data
dir_temp <- "./data/temp"

# isimip
isimip <- read_rds(file.path(dir_temp, "isimip_sbt_datras_hist_ssp585.rds")) %>%
  mutate(date = Date,
         year = Year, 
         pop = if_else(IcesArea == "4abc", "4bc", IcesArea),
         temp = as.numeric(isimip_sbt),
         source = "isimip") %>%
  select(date:source) 

# oras5
oras <- read_rds(file.path(dir_temp, "oras5_datras.rds")) %>%
  mutate(date = Date,
         year = Year, 
         pop = if_else(IcesArea == "4abc", "4bc", IcesArea),
         temp = oras_sbt,
         source = "oras5") %>%
  select(date:source)

# nemo-medusa
nm <- read_rds(file.path(dir_temp, "nemomedusa_datras.rds")) %>%
  mutate(date = Date,
         year = as.numeric(Year), 
         pop = if_else(IcesArea == "4abc", "4bc", IcesArea),
         temp = nemomedusa_sbt,
         source = "nemo-medusa") %>%
  select(date:source)

## merge all temp data
# merge and summarize mean(temp) 
data_temp <- bind_rows(isimip, oras, nm) %>%
  filter(pop %in% c("4bc", "7a", "8ab"),
         year >= min(data_otl$year), year <= max(data_otl$year)) %>%
  group_by(source, pop, year) %>%
  summarize(temp = mean(temp))

# get centered temp (c.temp) and average temperature (ave.temp) for each pop
data_temp <- data_temp %>%
  group_by(source, pop) %>%
  mutate(c.temp = c.(temp),
         ave.temp = mean(temp))

# add source name
df_temp_name <- tibble(source = c("oras5", "isimip", "nemo-medusa"),
                       source_name = factor(c("ORAS5", "ISIMIP", "NEMO-MEDUSA"), 
                                            levels = c("ORAS5", "ISIMIP", "NEMO-MEDUSA")))
data_temp <- data_temp %>% left_join(df_temp_name)

# add pop name
data_temp <- data_temp %>% left_join(df_pop)

## Fishing data ----
# Fishing mortality, Spawning Stock Biomass, Recruitment from ICES Stock Assessment

dir_ices <- "./data/ices"

# sole distribution area - survey datras
datras <- read_sf(file.path(dir_ices, "hl_loc_4abc7a8ab.gpkg"))
datras <- as_tibble(datras) %>%
  mutate(pop = if_else(Area_27 == "4abc", "4bc", Area_27)) %>%
  select(pop, area_km2)

# sole stock assessment
data_sol <- read_rds(file.path(dir_ices, "stock-assessment_2023.rds")) %>%
  rename(ssb = SSB,
         f = `F`) %>%
  select(pop, year, ssb, f, recruitment) %>%
  left_join(datras) %>%
  mutate(ssb.i = ssb/area_km2,
         recruitment.i = recruitment/area_km2) %>%
  select(-f) #remove fbar from stock assessment (different fbar age range across populations)

# add fbar from f-at-age 3-7 (consistent fbar age range across all populations)
data_sol_fbar <- read_rds(file.path(dir_ices, "sol_fbar_age3-7_stock-assessment_2023.rds"))
data_sol <- data_sol %>% left_join(data_sol_fbar)

## Nutrient data ----
dir_nu <- "./data/nutrient"
data_nu <- read_rds(file.path(dir_nu, "ospar_subset_1978-2017_ices_4abc.rds"))
# summarize all river by year 
data_nu <- data_nu %>% 
  group_by(IcesArea, year) %>% 
  summarize(TN = sum(TN)/1000,
            TP = sum(TP)/1000) %>%
  mutate(IcesArea = "4bc") %>%
  rename(pop = IcesArea)
# unit: 1000tN/year, 1000tP/year

## Model data ----

#### BEST INTRINSIC STRUCTURE: 
## no scaled
m3 <- read_rds(file.path(dir_output, "intrinsic.model_best.rds"))
## scaled
m3_s <- read_rds(file.path(dir_output, "intrinsic.model_best_scaled.rds"))

#### BEST EXTRINSIC STRUCTURE: 
## no scaled
m4_isimip <- read_rds(file.path(dir_output, "extrinsic.model_best_isimip.rds"))
m4_oras5 <- read_rds(file.path(dir_output, "extrinsic.model_best_oras5.rds"))

## scaled
m4_isimip_s <- read_rds(file.path(dir_output, "extrinsic.model_best_scaled_isimip.rds"))
m4_oras5_s <- read_rds(file.path(dir_output, "extrinsic.model_best_scaled_oras5.rds"))

#### BEST EXTRINSIC STRUCTURE WITH NUTRIENT DATA:
## no scaled
m4_nu_isimip <- read_rds(file.path(dir_output, "extrinsic.model.nu_best_isimip.rds"))
m4_nu_oras5 <- read_rds(file.path(dir_output, "extrinsic.model.nu_best_oras5.rds"))

## scaled
m4_nu_isimip_s <- read_rds(file.path(dir_output, "extrinsic.model.nu_best_scaled_isimip.rds"))
m4_nu_oras5_s <- read_rds(file.path(dir_output, "extrinsic.model.nu_best_scaled_oras5.rds"))

#### BEST EXTRINSIC EXTENDED STRUCTURE: 
## no scaled
m5_isimip <- read_rds(file.path(dir_output, "extrinsic.model.ext_best_isimip.rds"))
m5_oras5 <- read_rds(file.path(dir_output, "extrinsic.model.ext_best_oras5.rds"))

## scaled
m5_isimip_s <- read_rds(file.path(dir_output, "extrinsic.model.ext_best_scaled_isimip.rds"))
m5_oras5_s <- read_rds(file.path(dir_output, "extrinsic.model.ext_best_scaled_oras5.rds"))

#### BEST EXTRINSIC EXTENDED STRUCTURE WITH NUTRIENT DATA: 
## no scaled
m5_nu_isimip <- read_rds(file.path(dir_output, "extrinsic.model.ext.nu_best_isimip.rds"))
m5_nu_oras5 <- read_rds(file.path(dir_output, "extrinsic.model.ext.nu_best_oras5.rds"))

## scaled
m5_nu_isimip_s <- read_rds(file.path(dir_output, "extrinsic.model.ext.nu_best_scaled_isimip.rds"))
m5_nu_oras5_s <- read_rds(file.path(dir_output, "extrinsic.model.ext.nu_best_scaled_oras5.rds"))

# 3. FIGURE ----

## fig S1-5 - otolith sampling - Appendix S1. otolith (smartdots) ----
## fig S6 - otolith sampling - age at capture distribution -----
# setup
otl <- otl %>% 
  filter(FishID != "SOL_52_G1_Z.121_26-03-1990_1602") %>%
  mutate(pop.name = if_else(IcesAreaGroup == "4bc", "North Sea",
                            if_else(IcesAreaGroup == "7a", "Irish Sea", 
                                    "Bay of Biscay"))) %>%
  mutate(pop.name = factor(pop.name, levels = c("North Sea", "Irish Sea", "Bay of Biscay")))

otl_sum <- otl %>% group_by(pop.name, FishID, AgeAtCapture) %>% summarize(n = n()) %>%
  group_by(pop.name, AgeAtCapture) %>% summarize(n_fish = n())

# plot
ggplot(data = otl_sum, aes(x = AgeAtCapture, y = n_fish)) + 
  geom_bar(stat = "identity") +
  facet_wrap(~ pop.name) +
  labs(x = "Age at capture (years)",
       y = "Number of otoliths") +
  theme(axis.title = element_text(size = 9))

# save plot
ggsave(last_plot(), file = file.path(dir_report, "figS6_age at capture distribution.pdf"),
       device = cairo_pdf,
       width = 19, height = 6.3,
       units = "cm")
ggsave(last_plot(), file = file.path(dir_report, "figS6_age at capture distribution.png"),
       width = 19, height = 6.3,
       units = "cm",
       dpi = 1000)

## fig S7 - otolith sampling - sampling size ----
## sampling size vs SamplingYear/Cohort vs pop
otl_sum <- otl %>% group_by(pop.name, SamplingYear, Cohort, AgeAtCapture) %>% summarize(n_fish = n())

### SamplingYear
ggplot(data = otl_sum, aes(x = SamplingYear, y = AgeAtCapture, size = n_fish)) + 
  geom_point(alpha = 0.2) + 
  facet_wrap(~ pop.name) +
  labs(x = "Collection year",
       y = "Age at capture (years)",
       size = "Number of otoliths") +
  scale_size(breaks = c(10, 50, 100, 250)) +
  theme(legend.position = "bottom",
        legend.justification = c(0.5, 0.5)) +
  theme(axis.title = element_text(size = 9),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9))

# save plot
ggsave(last_plot(), file = file.path(dir_report, "figS7_sampling size.pdf"),
       device = cairo_pdf,
       width = 19, height = 6.3,
       units = "cm")
ggsave(last_plot(), file = file.path(dir_report, "figS7_sampling size.png"),
       width = 19, height = 6.3,
       units = "cm",
       dpi = 1000)

## fig S8 - otolith sampling - otolith size vs length ----
# setup
otl_sub <- data_otl %>% 
  group_by(FishID) %>%
  summarize(width = max(OtolithWidth.um),
            fish.length = max(Length.mm)) %>%
  filter(width >= 1500)

lm <- lm(fish.length ~ width, data = otl_sub)
summary(lm)
# Relationship between total fish length and otolith width (measured along the measurement axis on the otolith section) (Fish length ~ 51.42 + 0.076 * Otolith width, adjusted R2 = 0.53, p-value: < 0.001).
cor.test(otl_sub$fish.length, otl_sub$width)

ggplot(data = otl_sub, aes(x = width, y = fish.length)) +
  geom_point(alpha = 0.3) +
  geom_smooth(color = "black", method = "lm") +
  theme_classic() +
  labs(x = "Otolith diameter (μm)",
       y = "Fish length (mm)") +
  annotate("text", 
           x = 5500, y = 240, 
           label = expression("" ~ R^2 ~ " = 0.53, p < 0.001, n = 2152"), 
           size = 9/.pt) +
  theme(axis.title = element_text(size = 9))

# note: 
nrow(otl_sub) #2152 
# 2 fish with otolith width very small width (< 1500 um) likely due to error were not included
# R correlation: 0.73
# R2: 0.53

# save plot
ggsave(last_plot(), file = file.path(dir_report, "figS8_otolith vs fish size.pdf"),
       device = cairo_pdf,
       width = 19, height = 12.6,
       units = "cm")
ggsave(last_plot(), file = file.path(dir_report, "figS8_otolith vs fish size.png"),
       width = 19, height = 12.6,
       units = "cm",
       dpi = 1000)

## fig S9 - temperature - CTD by year (Appendix S4) ----
#### setup 
dir_gis <- "./data/admin"

# continents
continents <- read_sf(file.path(dir_gis, "esri_continent.gpkg"))

# countries
countries <- read_sf(file.path(dir_gis, "esri_countries.gpkg"))

# ices area 
ices_area <- read_sf(file.path(dir_gis, "ices_areas_sub_group_4abc_4326_new.gpkg")) 
ices_area <- ices_area %>% 
  filter(Area_27 %in% c("4abc", "7a", "8ab")) %>% 
  mutate(pop.name = if_else(Area_27 == "4abc", "North Sea",
                            if_else(Area_27 == "7a", "Irish Sea", 
                                    "Bay of Biscay"))) %>%
  st_simplify(dTolerance = 2000)

# ctd 
ctd <- read_sf(file.path(dir_temp, "ices_ctd_sea-bottom-temperature.gpkg"))
ctd <- ctd %>% 
  mutate(IcesArea = Area_27) %>%
  filter(IcesArea %in% c("7a", "8ab", "4abc")) %>%
  mutate(ctd = temp_degC, 
         date_match = as.Date(paste0("01", "-", Month, "-", Year),format = "%d-%m-%Y"))

#### plot 
st_bbox(ices_area)

list_plot <- list()

p1 <- ggplot() +
  geom_sf(data = countries, fill = "white", color = "grey",linewidth = 0.5) + 
  geom_sf(data = continents, fill = NA, color = "grey", linewidth = 1) + 
  geom_sf(data = ices_area, fill = NA, color = "black", linewidth = 1) +
  geom_sf(data = ctd %>% filter(Year == 1980), alpha = 0.3, size = 1) +
  coord_sf(xlim = c(-9, 9.5), ylim = c(43, 62.5), expand = FALSE) +
  facet_grid(~ Year)

p2 <- ggplot() +
  geom_sf(data = countries, fill = "white", color = "grey",linewidth = 0.5) + 
  geom_sf(data = continents, fill = NA, color = "grey", linewidth = 1) + 
  geom_sf(data = ices_area, fill = NA, color = "black", linewidth = 1) +
  geom_sf(data = ctd %>% filter(Year == 2000), alpha = 0.3, size = 1) +
  coord_sf(xlim = c(-9, 9.5), ylim = c(43, 62.5), expand = FALSE) +
  facet_grid(~ Year)

p3 <- ggplot() +
  geom_sf(data = countries, fill = "white", color = "grey",linewidth = 0.5) + 
  geom_sf(data = continents, fill = NA, color = "grey", linewidth = 1) + 
  geom_sf(data = ices_area, fill = NA, color = "black", linewidth = 1) +
  geom_sf(data = ctd %>% filter(Year == 2020), alpha = 0.3, size = 1) +
  coord_sf(xlim = c(-9, 9.5), ylim = c(43, 62.5), expand = FALSE) +
  facet_grid(~ Year)

p1 + p2 + p3

# save plot
ggsave(last_plot(), file = file.path(dir_report, "figS9_ctd by year.pdf"),
       device = cairo_pdf,
       width = 19, height = 12.6,
       units = "cm")
ggsave(last_plot(), file = file.path(dir_report, "figS9_ctd by year.png"),
       width = 19, height = 12.6,
       units = "cm",
       dpi = 1000)

## fig S10 - temperature - CTD vs temp (Appendix S4) ----
#### setup
# read temp data
oras_ras <- read_stars(file.path(dir_temp, "oras5.tif"))

ctd_sub <- ctd %>% 
  filter(Year == 2000, Month == 1, Day %in% c(23, 24)) %>%
  mutate(label = if_else(Day == 23, 
                         "CTD 23/01/2000 vs ORAS5 01/2000", 
                         "CTD 24/01/2000 vs ORAS5 01/2000"))

#### plot
bb <- st_bbox(c(xmin = 0, xmax = 6, ymax = 51, ymin = 55), crs = st_crs(4326))
breaks = c(6, 7, 8)

p1 <- ggplot() +
  geom_stars(data = oras_ras[,,,517][bb]) +
  geom_sf(data = continents, fill = "grey") + 
  geom_sf(data = countries , fill = NA, color = "black",linewidth = 0.5) + 
  geom_sf(data = ctd_sub %>% filter(Year == 2000, Month == 1, Day == 23)) +
  coord_sf(xlim = c(0, 6), ylim = c(51, 55), expand = FALSE) +
  scale_fill_distiller(palette = "RdYlBu", 
                       breaks = breaks, 
                       na.value = "white",
                       guide = guide_colorbar(title = "ORAS5 Temperature", 
                                              title.position = "top",
                                              title.hjust = 0.5,
                                              barwidth = 10)) +
  theme_bw() +
  labs(x = NULL,
       y = NULL) +
  theme(legend.position = "bottom") +
  facet_grid(~ label)

p2 <- ggplot() +
  geom_stars(data = oras_ras[,,,517][bb]) +
  geom_sf(data = continents, fill = "grey") + 
  geom_sf(data = countries, fill = NA, color = "black",linewidth = 0.5) + 
  geom_sf(data = ctd_sub %>% filter(Year == 2000, Month == 1, Day == 24)) +
  coord_sf(xlim = c(0, 6), ylim = c(51, 55), expand = FALSE) +
  scale_fill_distiller(palette = "RdYlBu", 
                       breaks = breaks, 
                       na.value = "white",
                       guide = guide_colorbar(title = "ORAS5 Temperature", 
                                              title.position = "top",
                                              title.hjust = 0.5,
                                              barwidth = 10)) +
  theme_bw() +
  labs(x = NULL,
       y = NULL) +
  theme(legend.position = "bottom") +
  facet_grid(~ label)
  
(p1 + p2) + 
  plot_layout(axis_titles = "collect", guides = "collect") & 
  theme(legend.position = "bottom")

# save plot
ggsave(last_plot(), file = file.path(dir_report, "figS10_ctd vs temp.pdf"),
       device = cairo_pdf,
       width = 19, height = 12.6,
       units = "cm")
ggsave(last_plot(), file = file.path(dir_report, "figS10_ctd vs temp.png"),
       width = 19, height = 12.6,
       units = "cm",
       dpi = 1000)

## fig S11 - distribution area (Appendix S5) ----
#### setup
dir_datras <- "./data/ices"

hh <- read_rds(file.path(dir_datras, "sf_HHflats_1985till2022_in3a204a4b4c7a7d7e7f7g7h7j28a8b.RDS"))
hl <- read_rds(file.path(dir_datras, "sf_Solea_soleaHL_withAbs_BTS+BTS-VIII+DYFS+SNS_in1985till2021_in3a204a4b4c7a7d7e7f7g7h7j28a8b.RDS"))

hh_loc <- hh %>% 
  select(Survey, Area_27) %>%
  filter(Area_27 %in% c("4.a", "4.b", "4.c", "7.a", "8.a", "8.b")) %>%
  unique()

hl_loc <- hl %>% 
  filter(!is.na(HLNoAtLngt),
         #Quarter %in% c(1,3,4),
         Area_27 %in% c("4.a", "4.b", "4.c", "7.a", "8.a", "8.b")) %>%
  select(Survey, Area_27) %>%
  unique()

# gis
dir_gis <- "./data/admin"
# continents
continents <- read_sf(file.path(dir_gis, "esri_continent.gpkg"))
# countries
countries <- read_sf(file.path(dir_gis, "esri_countries.gpkg"))
# ices area 
ices_area <- read_sf(file.path(dir_gis, "ices_areas_sub_group_4abc_4326_new.gpkg")) 
ices_area <- ices_area %>% 
  filter(Area_27 %in% c("4abc", "7a", "8ab")) %>% 
  mutate(pop.name = if_else(Area_27 == "4abc", "North Sea",
                            if_else(Area_27 == "7a", "Irish Sea", 
                                    "Bay of Biscay"))) %>%
  st_simplify(dTolerance = 2000)

#### plot
p1 <- ggplot() +
  geom_sf(data = countries, fill = "white", color = "grey",linewidth = 0.5) + 
  geom_sf(data = continents, fill = NA, color = "grey", linewidth = 1) + 
  geom_sf(data = ices_area, fill = NA, color = "black", linewidth = 1) +
  geom_sf(data = hh_loc, aes(color = Survey)) +
  coord_sf(xlim = c(-9, 9.5), ylim = c(43, 62.5), expand = FALSE) +
  theme_bw()

p2 <- ggplot() +
  geom_sf(data = countries, fill = "white", color = "grey",linewidth = 0.5) + 
  geom_sf(data = continents, fill = NA, color = "grey", linewidth = 1) + 
  geom_sf(data = ices_area, fill = NA, color = "black", linewidth = 1) +
  geom_sf(data = hl_loc, aes(color = Survey)) +
  coord_sf(xlim = c(-9, 9.5), ylim = c(43, 62.5), expand = FALSE) +
  theme_bw()

(p1 + p2) + 
  plot_annotation(tag_levels = 'A') +
  plot_layout(guides = "collect")  &
  theme(legend.position = "bottom",
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9))

# save plot
ggsave(last_plot(), file = file.path(dir_report, "figS11_datras survey.pdf"),
       device = cairo_pdf,
       width = 19, height = 12.6,
       units = "cm")
ggsave(last_plot(), file = file.path(dir_report, "figS11_datras survey.png"),
       width = 19, height = 12.6,
       units = "cm",
       dpi = 1000)

## fig S12 - ORAS5 vs ISIMIP trend (Appendix S6) ----
list_pop <- c("4bc", "7a", "8ab")
list_pop.name <- c("North Sea", "Irish Sea", "Bay of Biscay")

df <- tibble()
for (i in 1:3) {
  lm_oras5 <- lm(temp ~ year, data = filter(data_temp, source == "oras5", pop == list_pop[i]))
  lm_isimip <- lm(temp ~ year, data = filter(data_temp, source == "isimip", pop == list_pop[i]))
  
  df_oras5 <- tibble(pop = list_pop[i],
                     est = lm_oras5$coefficients["year"],
                     est_lower = confint(lm_oras5)["year",1],
                     est_upper = confint(lm_oras5)["year",2],
                     source = "ORAS5")
  
  df_isimip <- tibble(pop = list_pop[i],
                      est = lm_isimip$coefficients["year"],
                      est_lower = confint(lm_isimip)["year",1],
                      est_upper = confint(lm_isimip)["year",2],
                      source = "ISIMIP")
  
  df <- bind_rows(df, df_oras5, df_isimip)
}
df <- df %>%
  left_join(df_pop)

ggplot(data = df) +
  geom_point(aes(x = est, y = source)) +
  geom_linerange(aes(x = est, y = source, xmin = est_lower, xmax = est_upper)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_grid(~ pop.name) +
  labs(x = "Estimate of linear year effect",
       y = NULL)

#### save plot
ggsave(last_plot(), file = file.path(dir_report, "figS12_oras5 vs isimip_trend.pdf"),
       device = cairo_pdf,
       width = 9, height = 6.3,
       units = "cm")
ggsave(last_plot(), file = file.path(dir_report, "figS12_oras5 vs isimip_trend.png"),
       width = 9, height = 6.3,
       units = "cm",
       dpi = 1000)

## fig S13 - temperature - temp ----
#### plot
ggplot(data = data_temp %>%
         filter(source %in% c("isimip", "oras5")) %>%
         left_join(df_pop), 
       aes(x = year, y = temp, color = pop.name)) +
  geom_line() +
  facet_grid(source_name ~ .) +
  labs(x = "Year",
       y = "Temperature (°C)",
       color = "Population") +
  theme(legend.position = "bottom") +
  scale_color_brewer(palette = "Dark2", direction = -1) +
  theme(axis.title = element_text(size = 9))

#### save plot
ggsave(last_plot(), file = file.path(dir_report, "figS13_temperature.pdf"),
       device = cairo_pdf,
       width = 19, height = 12.6,
       units = "cm")
ggsave(last_plot(), file = file.path(dir_report, "figS13_temperature.png"),
       width = 19, height = 12.6,
       units = "cm",
       dpi = 1000)

## fig S14 - temperature - c.temp ----
#### plot
ggplot(data = data_temp %>% 
         filter(source %in% c("isimip", "oras5")) %>%
         left_join(df_pop), 
       aes(x = year, y = c.temp, color = source_name)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(pop.name ~ .) +
  labs(x = "Year",
       y = expression('T'['population-anomaly'] * ' (°C)'),
       color = "Dataset") +
  theme(legend.position = "bottom") +
  scale_color_manual(values = col_scale) +
  theme(axis.title = element_text(size = 9),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9))

#### save plot
ggsave(last_plot(), file = file.path(dir_report, "figS14_temperature pop anomaly.pdf"),
       device = cairo_pdf,
       width = 19, height = 12.6,
       units = "cm")
ggsave(last_plot(), file = file.path(dir_report, "figS14_temperature pop anomaly.png"),
       width = 19, height = 12.6,
       units = "cm",
       dpi = 1000)


## fig S15 - variable - spawning stock biomass, recruitment ---- 
data_sol <- data_sol %>% left_join(df_pop)

#### plot
p1 <- ggplot(data = data_sol, aes(x = year, y = ssb, color = pop.name)) +
  geom_line() +
  labs(x = "Year",
       y = "Spawning stock biomass (tonne)",
       color = "Population") +
  theme(legend.position = "none",
        legend.justification = c(0, 1),
        legend.title = element_blank()) +
  scale_color_brewer(palette = "Dark2", direction = -1) +
  xlim(1958, 2019)

p2 <- ggplot(data = data_sol, aes(x = year, y = recruitment, color = pop.name)) +
  geom_line() +
  labs(x = "Year",
       y = "Recruitment (thoudsand)",
       color = "Population") +
  theme(legend.position = c(0.6, 1.05),
        legend.justification = c(0.1, 1),
        legend.background = element_rect(fill = "transparent"),
        legend.title = element_blank()) +
  scale_color_brewer(palette = "Dark2", direction = -1) +
  xlim(1958, 2019)

p1 + p2 + 
  plot_annotation(tag_levels = 'A') +
  theme(axis.title = element_text(size = 9))

#### save plot
ggsave(last_plot(), file = file.path(dir_report, "figS15_raw ssb rec.pdf"),
       device = cairo_pdf,
       width = 19, height = 6.3,
       units = "cm")
ggsave(last_plot(), file = file.path(dir_report, "figS15_raw ssb rec.png"),
       width = 19, height = 6.3,
       units = "cm",
       dpi = 1000)

## fig S16 - variable - TN and TP ----
#### plot
p1 <- ggplot() + 
  geom_line(data = data_nu, aes(x = year, y = TN)) +
  labs(x = "Year",
       y = "Total nitrogen (kilotonne/year)")  +
  xlim(1958, 2019)

p2 <- ggplot() + 
  geom_line(data = data_nu, aes(x = year, y = TP)) +
  labs(x = "Year",
       y = "Total phosphorus (kilotonne/year)")  +
  xlim(1958, 2019)

p1 + p2 + 
  plot_annotation(tag_levels = 'A') +
  theme(axis.title = element_text(size = 9))

#### save plot
ggsave(last_plot(), file = file.path(dir_report, "figS16_nutrient.pdf"),
       device = cairo_pdf,
       width = 19, height = 6.3,
       units = "cm")
ggsave(last_plot(), file = file.path(dir_report, "figS16_nutrient.png"),
       width = 19, height = 6.3,
       units = "cm",
       dpi = 1000)

## fig S17,18,19 - model diagnostics ----
### intrinsic model ----
model <- m3
data <- model@frame
data <- data %>% 
  mutate(res = resid(model),
         fit = predict(model)) %>%
  mutate(col = "Intrinsic model")

ggplot(data = data, aes(x = fit, y = res)) +
  geom_point(alpha = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Fitted values",
       y = "Residuals") +
  facet_grid(~ col)
p1 <- last_plot()

ggplot(data = data, aes(sample = res)) +
  geom_qq(alpha = 0.1) +
  geom_qq_line(linetype = "dashed") +
  labs(x = "Theoretical quantiles",
       y = "Sample quantiles") +
  facet_grid(~ col)
p2 <- last_plot()

### extrinsic model ----
#### oras5 ----
model <- m4_oras5
data <- model@frame
data <- data %>% 
  mutate(res = resid(model),
         fit = predict(model)) %>%
  mutate(col = "Population-level extrinsic model (ORAS5)")

ggplot(data = data, aes(x = fit, y = res)) +
  geom_point(alpha = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Fitted values",
       y = "Residuals") +
  facet_grid(~ col)
p3 <- last_plot()

ggplot(data = data, aes(sample = res)) +
  geom_qq(alpha = 0.1) +
  geom_qq_line(linetype = "dashed") +
  labs(x = "Theoretical quantiles",
       y = "Sample quantiles") +
  facet_grid(~ col)
p4 <- last_plot()

#### isimip ----
model <- m4_isimip
data <- model@frame
data <- data %>% 
  mutate(res = resid(model),
         fit = predict(model)) %>%
  mutate(col = "Population-level extrinsic model (ISIMIP)")

ggplot(data = data, aes(x = fit, y = res)) +
  geom_point(alpha = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Fitted values",
       y = "Residuals") +
  facet_grid(~ col)
p5 <- last_plot()

ggplot(data = data, aes(sample = res)) +
  geom_qq(alpha = 0.1) +
  geom_qq_line(linetype = "dashed") +
  labs(x = "Theoretical quantiles",
       y = "Sample quantiles") +
  facet_grid(~ col)
p6 <- last_plot()

### extrinsic extended model ----
#### oras5 ----
model <- m5_oras5
data <- model@frame
data <- data %>% 
  mutate(res = resid(model),
         fit = predict(model)) %>%
  mutate(col = "Individual-level extrinsic model (ORAS5)")

ggplot(data = data, aes(x = fit, y = res)) +
  geom_point(alpha = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Fitted values",
       y = "Residuals") +
  facet_grid(~ col)
p7 <- last_plot()

ggplot(data = data, aes(sample = res)) +
  geom_qq(alpha = 0.1) +
  geom_qq_line(linetype = "dashed") +
  labs(x = "Theoretical quantiles",
       y = "Sample quantiles") +
  facet_grid(~ col)
p8 <- last_plot()

#### isimip ----
model <- m5_isimip
data <- model@frame
data <- data %>% 
  mutate(res = resid(model),
         fit = predict(model)) %>%
  mutate(col = "Individual-level extrinsic model (ISIMIP)")

ggplot(data = data, aes(x = fit, y = res)) +
  geom_point(alpha = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Fitted values",
       y = "Residuals") +
  facet_grid(~ col)
p9 <- last_plot()

ggplot(data = data, aes(sample = res)) +
  geom_qq(alpha = 0.1) +
  geom_qq_line(linetype = "dashed") +
  labs(x = "Theoretical quantiles",
       y = "Sample quantiles") +
  facet_grid(~ col)
p10 <- last_plot()

### merge and save plot ----
#### intrinsic ----
(p1 | p2) 

ggsave(last_plot(), file = file.path(dir_report, "figS17_model diagnostics_intrinsic model.pdf"),
       device = cairo_pdf,
       width = 19, height = 6.3,
       units = "cm")
ggsave(last_plot(), file = file.path(dir_report, "figS17_model diagnostics_intrinsic model.png"),
       width = 19, height = 6.3,
       units = "cm",
       dpi = 1000)

#### extrinsic ----
(p3 | p4)/(p5 | p6) 
ggsave(last_plot(), file = file.path(dir_report, "figS18_model diagnostics_extrinsic model.pdf"),
       device = cairo_pdf,
       width = 19, height = 12.6,
       units = "cm")
ggsave(last_plot(), file = file.path(dir_report, "figS18_model diagnostics_extrinsic model.png"),
       width = 19, height = 12.6,
       units = "cm",
       dpi = 1000)

#### extrinsic extended ----
(p7 | p8)/(p9 | p10) 
ggsave(last_plot(), file = file.path(dir_report, "figS19_model diagnostics_extrinsic extended model.pdf"),
       device = cairo_pdf,
       width = 19, height = 12.6,
       units = "cm")
ggsave(last_plot(), file = file.path(dir_report, "figS19_model diagnostics_extrinsic extended model.png"),
       width = 19, height = 12.6,
       units = "cm",
       dpi = 1000)

## fig S20 - ave.temp, c.temp, c.temp_within, c.temp_between ----

list_plot <- list()

### ave.temp ----
#### setup
# note: fit within range [9,12.8] as temperatures has different range 

m4_oras5_temp <- read_rds(file.path(dir_output, "extrinsic.model_best_oras5_temp.rds"))

list_model <- list(
  "isimip" = m4_isimip,
  "oras5" = m4_oras5_temp)

pred <- tibble()
for (i in 1:length(list_model)) {
  model <- list_model[[i]]
  data <- model@frame
  temp_source <- names(list_model[i])
  
  pred_temp <- as.data.frame(Effect(c("ave.temp"), 
                                    model, 
                                    xlevels = list(ave.temp = seq(9, 12.8, 0.1)))) %>%
    mutate(source = temp_source)
  
  pred <- bind_rows(pred, pred_temp)
}

pred <- pred %>% 
  left_join(df_temp_name)

#### plot
yrange <- exp(c(min(pred$lower), max(pred$upper)))
  
ggplot(data = pred %>% filter(source_name == "ORAS5")) + 
  geom_line(aes(x = ave.temp, 
                y = exp(fit)),
            linewidth = 1,
            color = "grey") +
  geom_ribbon(aes(x = ave.temp,
                  ymin = exp(lower),
                  ymax = exp(upper)), 
              alpha = 0.1) +
  facet_grid( ~ source_name) +
  labs(x = expression('T'['population-average'] * ' (°C)'),
       y = "Predicted increment growth (μm)",
       color = "Dataset",
       fill = "Dataset") +
  ylim(yrange)
p1 <- last_plot()

ggplot(data = pred %>% filter(source_name == "ISIMIP")) + 
  geom_line(aes(x = ave.temp, 
                y = exp(fit)),
            linewidth = 1) +
  geom_ribbon(aes(x = ave.temp,
                  ymin = exp(lower),
                  ymax = exp(upper)), 
              alpha = 0.1) +
  facet_grid( ~ source_name) +
  labs(x = expression('T'['population-average'] * ' (°C)'),
       y = "Predicted increment growth (μm)",
       color = "Dataset",
       fill = "Dataset") +
  ylim(yrange)
p2 <- last_plot()

### c.temp ----
#### setup
# note: fit within range [-1,1] as temperatures has different range 

list_model <- list(
  "isimip" = m4_isimip,
  "oras5" = m4_oras5_temp
)

pred <- tibble()
for (i in 1:length(list_model)) {
  model <- list_model[[i]]
  data <- model@frame
  temp_source <- names(list_model[i])
  
  pred_temp <- as.data.frame(Effect(c("c.temp"), 
                                    model, 
                                    xlevels = list(c.temp = seq(-1, 1, 0.1)))) %>%
    mutate(source = temp_source)
  
  pred <- bind_rows(pred, pred_temp)
}

pred <- pred %>% 
  left_join(df_temp_name)

#### plot
yrange <- exp(c(min(pred$lower), max(pred$upper)))

ggplot(data = pred %>% filter(source_name == "ORAS5")) + 
  geom_line(aes(x = c.temp, 
                y = exp(fit)),
            linewidth = 1,
            color = "grey") +
  geom_ribbon(aes(x = c.temp,
                  ymin = exp(lower),
                  ymax = exp(upper)), 
              alpha = 0.1) +
  labs(x = expression('T'['population-anomaly'] * ' (°C)'),
       y = "Predicted increment growth (μm)") +
  ylim(yrange)
p3 <- last_plot()

ggplot(data = pred %>% filter(source_name == "ISIMIP")) + 
  geom_line(aes(x = c.temp, 
                y = exp(fit)),
            linewidth = 1) +
  geom_ribbon(aes(x = c.temp,
                  ymin = exp(lower),
                  ymax = exp(upper)), 
              alpha = 0.1) +
  labs(x = expression('T'['population-anomaly'] * ' (°C)'),
       y = "Predicted increment growth (μm)") +
  ylim(yrange)
p4 <- last_plot()

### c.temp_between ----
#### setup
# note: fit within range [-1,1] as temperatures has different range 
# note: oras5 has smaller range [-0.8,0.9]

list_model <- list("isimip" = m5_isimip,
                   "oras5" = m5_oras5)

pred <- tibble()
for (i in 1:length(list_model)) {
  model <- list_model[[i]]
  data <- model@frame
  temp_source <- names(list_model[i])
  
  pred_temp <- as.data.frame(Effect(c("log.age", "c.temp_between"), 
                                    model, 
                                    xlevels = list(log.age = unique(data$log.age),
                                                   c.temp_between = seq(-1, 1, 0.1)))) %>%
    mutate(age = round(exp(log.age),1)) %>%
    mutate(source = temp_source)
  
  pred <- bind_rows(pred, pred_temp)
}

pred <- pred %>% 
  left_join(df_temp_name)

#### plot
yrange <- exp(c(min(pred$lower), max(pred$upper)))

ggplot(data = pred %>% filter(age < 6, source_name == "ORAS5")) + 
  geom_line(aes(x = c.temp_between, 
                y = exp(fit),
                linetype = factor(age)),
            linewidth = 1,
            alpha = 0.3) +
  geom_ribbon(aes(x = c.temp_between,
                  ymin = exp(lower),
                  ymax = exp(upper),
                  group = factor(age)), 
              alpha = 0.1) +
  labs(x = expression('T'['individual-average'] * ' (°C)'),
       y = "Predicted increment growth (μm)",
       linetype = "Age") +
  scale_linetype_manual(values = seq(2,6)) +
  guides(linetype = guide_legend(override.aes = list(alpha = 1)))
p5 <- last_plot()

ggplot(data = pred %>% filter(age < 6, source_name == "ISIMIP")) + 
  geom_line(aes(x = c.temp_between, 
                y = exp(fit),
                linetype = factor(age)),
            linewidth = 1) +
  geom_ribbon(aes(x = c.temp_between,
                  ymin = exp(lower),
                  ymax = exp(upper),
                  group = factor(age)), 
              alpha = 0.1) +
  labs(x = expression('T'['individual-average'] * ' (°C)'),
       y = "Predicted increment growth (μm)",
       linetype = "Age") +
  scale_linetype_manual(values = seq(2,6)) 
p6 <- last_plot()

### c.temp_within ----
#### setup
# note: fit within range [-1,1] as temperatures has different range 

list_model <- list("isimip" = m5_isimip,
                   "oras5" = m5_oras5)

pred <- tibble()
for (i in 1:length(list_model)) {
  model <- list_model[[i]]
  data <- model@frame
  temp_source <- names(list_model[i])
  
  pred_temp <- as.data.frame(Effect(c("log.age", "c.temp_within"), 
                                    model, 
                                    xlevels = list(log.age = unique(data$log.age),
                                                   c.temp_within = seq(-1, 1, 0.1)))) %>%
    mutate(age = round(exp(log.age),1)) %>%
    mutate(source = temp_source)
  
  pred <- bind_rows(pred, pred_temp)
}

pred <- pred %>% 
  left_join(df_temp_name)

#### plot
yrange <- exp(c(min(pred$lower), max(pred$upper)))

ggplot(data = pred %>% filter(age < 6, source_name == "ORAS5")) + 
  geom_line(aes(x = c.temp_within, 
                y = exp(fit),
                linetype = factor(age)),
            linewidth = 1) +
  geom_ribbon(aes(x = c.temp_within,
                  ymin = exp(lower),
                  ymax = exp(upper),
                  group = factor(age)), 
              alpha = 0.1) +
  labs(x = expression('T'['individual-anomaly'] * ' (°C)'),
       y = "Predicted increment growth (μm)",
       linetype = "Age") +
  ylim(yrange) +
  scale_linetype_manual(values = seq(2,6)) 
p7 <- last_plot()

ggplot(data = pred %>% filter(age < 6, source_name == "ISIMIP")) + 
  geom_line(aes(x = c.temp_within, 
                y = exp(fit),
                linetype = factor(age)),
            linewidth = 1) +
  geom_ribbon(aes(x = c.temp_within,
                  ymin = exp(lower),
                  ymax = exp(upper),
                  group = factor(age)), 
              alpha = 0.1) +
  labs(x = expression('T'['individual-anomaly'] * ' (°C)'),
       y = "Predicted increment growth (μm)",
       linetype = "Age") +
  ylim(yrange) +
  scale_linetype_manual(values = seq(2,6)) 
p8 <- last_plot()

### merge and save plot ----
# merge plot
(p1 + p2 + 
   p3 + p4 +
   p5 + p6 +
   p7 + p8) +
  plot_layout(nrow = 4, axis_titles = "collect_y", guides = "collect") +
  plot_annotation(tag_levels = 'A')  &  
  theme(legend.position = "bottom", 
        plot.tag.position  = c(0.92, 0.82), # c(0.95, 0.85)
        axis.title = element_text(size = 9),
        legend.key.size = unit(10, "pt"),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9))


(p0 + p0 +
    p1 + p2 +
    p3 + p4 +
    p5 + p6) +
  plot_layout(nrow = 4, axis_titles = "collect", guides = "collect") +
  plot_annotation(tag_levels = 'A') &  
  theme(legend.position = "bottom", 
        plot.tag.position  = c(0.92, 0.82), # c(0.95, 0.85)
        axis.title = element_text(size = 9),
        legend.key.size = unit(10, "pt"),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9))

# save file
ggsave(last_plot(), file = file.path(dir_report, "figS20_temp effects.pdf"),
       device = cairo_pdf,
       width = 19, height = 28,
       units = "cm")
ggsave(last_plot(), file = file.path(dir_report, "figS20_temp effects.png"),
       width = 19, height = 28,
       units = "cm",
       dpi = 1000)

## fig S21 - density effect with and without standardisation ----
#### save file
df_all <- read_rds(file.path(dir_output, "robustness_density without standardisation_extrinsic.model.rds"))

df_all <- df_all %>%
  mutate(source_name = factor(source_name, levels = c("ORAS5", "ISIMIP")))

#### plot
ggplot(data = df_all, 
       aes(x = estimate, 
           y = data_name)) +
  geom_point() +
  geom_linerange(aes(xmin = estimate - 1.96*std.error,
                     xmax = estimate + 1.96*std.error)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_grid(source_name ~ term_name) +
  labs(x = "Parameter estimate",
       y = NULL) +
  scale_x_continuous(breaks = seq(0,2.5,0.5)) +
  theme(axis.title = element_text(size = 9))

#### save plot
ggsave(last_plot(), file = file.path(dir_report, "figS21_density effect vs standardisation.pdf"),
       device = cairo_pdf,
       width = 19, height = 12.6,
       units = "cm")
ggsave(last_plot(), file = file.path(dir_report, "figS21_density effect vs standardisation.png"),
       width = 19, height = 12.6,
       units = "cm", 
       dpi = 1000)

## fig S22 - extrinsic effects ----
list_plot <- list()

### ssb.i ----
list_model <- list("isimip" = m5_isimip,
                   "oras5" = m5_oras5
)
# 0.02-0.46 - nemo-medusa different range ssb.i
pred <- tibble()
for (i in 1:length(list_model)) {
  model <- list_model[[i]]
  data <- model@frame
  temp_source <- names(list_model[i])
  
  pred_temp <- as.data.frame(Effect(c("ssb.i"), 
                                    model, 
                                    xlevels = list(ssb.i = seq(0.02, 0.46, 0.01)))) %>%
    mutate(source = temp_source)
  
  pred <- bind_rows(pred, pred_temp)
}

pred <- pred %>% 
  left_join(df_temp_name)

#### plot
col_scale <- c("ISIMIP" = "#F8766D", "ORAS5" = "#00BFC4", "NEMO-MEDUSA" = "#7CAE00")

ggplot(data = pred) + 
  geom_line(aes(x = ssb.i, 
                y = exp(fit),
                color = source_name),
            linewidth = 1) +
  geom_ribbon(aes(x = ssb.i,
                  ymin = exp(lower),
                  ymax = exp(upper),
                  fill = source_name), 
              alpha = 0.1) +
  labs(x = "Spawning stock biomass (tonne/km²)",
       y = "Predicted increment growth (μm)",
       color = "Dataset",
       fill = "Dataset") +
  scale_color_manual(values = col_scale) +
  scale_fill_manual(values = col_scale) 

p1 <- last_plot()

### recruitment.i ----
list_model <- list("isimip" = m5_isimip,
                   "oras5" = m5_oras5)
pred <- tibble()
for (i in 1:length(list_model)) {
  model <- list_model[[i]]
  data <- model@frame
  temp_source <- names(list_model[i])
  
  pred_temp <- as.data.frame(Effect(c("recruitment.i"), 
                                    model, 
                                    xlevels = list(recruitment.i = seq(0.01, 2.53, 0.01)))) %>%
    mutate(source = temp_source)
  
  pred <- bind_rows(pred, pred_temp)
}

pred <- pred %>% 
  left_join(df_temp_name)

#### plot
col_scale <- c("ISIMIP" = "#F8766D", "ORAS5" = "#00BFC4", "NEMO-MEDUSA" = "#7CAE00")

ggplot(data = pred) + 
  geom_line(aes(x = recruitment.i, 
                y = exp(fit),
                color = source_name),
            linewidth = 1) +
  geom_ribbon(aes(x = recruitment.i,
                  ymin = exp(lower),
                  ymax = exp(upper),
                  fill = source_name), 
              alpha = 0.1) +
  labs(x = "Recruitment (thoudsand/km²)",
       y = "Predicted increment growth (μm)",
       color = "Dataset",
       fill = "Dataset") +
  scale_color_manual(values = col_scale) +
  scale_fill_manual(values = col_scale)

p2 <- last_plot()

### nutrient ----

list_model <- list("isimip" = m5_nu_isimip,
                   "oras5" = m5_nu_oras5)

pred <- tibble()
for (i in 1:length(list_model)) {
  model <- list_model[[i]]
  data <- model@frame
  temp_source <- names(list_model[i])
  
  pred_temp <- as.data.frame(Effect(c("TP"), 
                                    model, 
                                    xlevels = list(TP = unique(data$TP)))) %>%
    mutate(source = temp_source)
  
  pred <- bind_rows(pred, pred_temp)
}

pred <- pred %>% 
  left_join(df_temp_name)

#### plot
col_scale <- c("ISIMIP" = "#F8766D", "ORAS5" = "#00BFC4", "NEMO-MEDUSA" = "#7CAE00")

ggplot(data = pred) + 
  geom_line(aes(x = TP, 
                y = exp(fit),
                color = source_name),
            linewidth = 1) +
  geom_ribbon(aes(x = TP,
                  ymin = exp(lower),
                  ymax = exp(upper),
                  fill = source_name), 
              alpha = 0.1) +
  labs(x = "Total phosphorus (kilotonne/year)",
       y = "Predicted increment growth (μm)",
       color = "Dataset",
       fill = "Dataset") +
  scale_color_manual(values = col_scale) +
  scale_fill_manual(values = col_scale)

p3 <- last_plot()

### merge and save plot ----
# merge plot
(p1 + p2 + p3) +
  plot_layout(nrow = 3, axis_titles = "collect", guides = "collect") +
  plot_annotation(tag_levels = 'A') &  
  theme(legend.position = "bottom", 
        plot.tag.position  = c(0.95, 0.9),
        axis.title = element_text(size = 9),
        legend.key.size = unit(10, "pt"),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9)
  )

# save file
ggsave(last_plot(), file = file.path(dir_report, "figS22_extrinsic effects.pdf"),
       device = cairo_pdf,
       width =  9, height = 18,
       units = "cm")
ggsave(last_plot(), file = file.path(dir_report, "figS22_extrinsic effects.png"),
       width =  9, height = 18,
       units = "cm",
       dpi = 1000)

## fig S23 - interaction c.temp_within vs env ----

list_plot <- list()

### ssb.i ----
# range: p25, 50, 75
df_age <- tibble(age = c(1, 3, 5),
                 age_name = c("Age 1", "Age 3", "Age 5"))

# ssb - isimip
model <- m5_isimip
data <- model@frame

df_ssb <- tibble(p_name = c("P25", "P75"),
                 ssb.i = round(quantile(data$ssb.i, probs = c(0.25, 0.75)),3)) 

pred <- as.data.frame(Effect(c("log.age", "c.temp_within", "ssb.i"), 
                             model, 
                             xlevels = list(log.age = unique(data$log.age),
                                            c.temp_within = seq(-1, 1, 0.1),
                                            ssb.i = df_ssb$ssb.i))) %>%
  mutate(age = round(exp(log.age),1)) %>%
  mutate(source = "isimip")

pred <- pred %>% 
  left_join(df_age) %>%
  left_join(df_temp_name) %>%
  left_join(df_ssb)

#### plot
# add a point of ORAS5, NEMO-MEDUSA to pred to merge legend 
pred_add <- pred %>% 
  filter(age %in% c(1,3,5), 
         c.temp_within == 0) %>%
  mutate(source_name = c("ORAS5"))

pred <- bind_rows(pred, pred_add) %>%
  mutate(source_name = factor(source_name, levels = c("ISIMIP", "ORAS5", "NEMO-MEDUSA"))) %>%
  mutate(var = "Spawning stock biomass")

col_scale <- c("ISIMIP" = "#F8766D", "ORAS5" = "#00BFC4", "NEMO-MEDUSA" = "#7CAE00")


list_plot <- list()
list_age <- c(1,3,5)
for(i in 1:3) {
  
  ggplot(data = pred %>% filter(age %in% list_age[i])) + 
    geom_line(aes(x = c.temp_within, 
                  y = exp(fit),
                  color = source_name,
                  linetype = p_name),
              linewidth = 1) +
    geom_ribbon(aes(x = c.temp_within,
                    ymin = exp(lower),
                    ymax = exp(upper),
                    fill = source_name,
                    linetype = p_name),
                alpha = 0.1) +
    facet_grid(. ~ age_name)  +
    labs(x = expression('T'['individual-anomaly'] * ' (°C)'),
         y = "Predicted increment growth (μm)",
         color = "Dataset",
         fill = "Dataset",
         linetype = "Variable value") +
    scale_color_manual(values = col_scale) +
    scale_fill_manual(values = col_scale) +
    guides(
      colour = guide_legend(order = 1),
      fill = guide_legend(order = 1),
      linedtype = guide_legend(order = 2)
    )
  
  # add facet variable name at the last panel
  if(i == 3) {
    ggplot(data = pred %>% filter(age %in% list_age[i])) + 
      geom_line(aes(x = c.temp_within, 
                    y = exp(fit),
                    color = source_name,
                    linetype = p_name),
                linewidth = 1) +
      geom_ribbon(aes(x = c.temp_within,
                      ymin = exp(lower),
                      ymax = exp(upper),
                      fill = source_name,
                      linetype = p_name),
                  alpha = 0.1) +
      facet_grid(var ~ age_name)  +
      labs(x = expression('T'['individual-anomaly'] * ' (°C)'),
           y = "Predicted increment growth (μm)",
           color = "Dataset",
           fill = "Dataset",
           linetype = "Variable value") +
      scale_color_manual(values = col_scale) +
      scale_fill_manual(values = col_scale) +
      guides(
        colour = guide_legend(order = 1),
        fill = guide_legend(order = 1),
        linedtype = guide_legend(order = 2)
      )
  }
  
  list_plot[[i]] <- last_plot()
}

p1 <- list_plot[[1]]
p2 <- list_plot[[2]]
p3 <- list_plot[[3]]

### recuitment.i ----
# recruitment - oras5, nemo-medusa
# range: p25, 50, 75
df_age <- tibble(age = c(1,3,5),
                 age_name = c("Age 1", "Age3", "Age 5"))

list_model <- list("oras5" = m5_oras5)

pred <- tibble()
for (i in 1:length(list_model)) {
  model <- list_model[[i]]
  data <- model@frame
  temp_source <- names(list_model[i])
  
  df_recruitment <- tibble(p_name = c("P25", "P75"),
                           recruitment.i = round(quantile(data$recruitment.i, probs = c(0.25, 0.75)),3)) 
  
  pred_temp <- as.data.frame(Effect(c("log.age", "c.temp_within", "recruitment.i"), 
                                    model, 
                                    xlevels = list(log.age = unique(data$log.age),
                                                   c.temp_within = seq(-1, 1, 0.1),
                                                   recruitment.i = df_recruitment$recruitment.i))) %>%
    mutate(age = round(exp(log.age),1)) %>%
    mutate(source = temp_source) %>%
    left_join(df_recruitment)
  
  pred <- bind_rows(pred, pred_temp)
}

pred <- pred %>% 
  left_join(df_age) %>%
  left_join(df_temp_name)

#### plot
# add a point of ISIMIP, ORAS5 to pred to merge legend 
# add a point of ORAS5, NEMO-MEDUSA to pred to merge legend 
pred_add <- pred %>% 
  filter(age %in% c(1,3,5), 
         c.temp_within == 0) %>%
  mutate(source_name = c("ISIMIP"))

pred <- bind_rows(pred, pred_add) %>%
  mutate(source_name = factor(source_name, levels = c("ISIMIP", "ORAS5", "NEMO-MEDUSA"))) %>%
  mutate(var = "Recruitment")

col_scale <- c("ISIMIP" = "#F8766D", "ORAS5" = "#00BFC4", "NEMO-MEDUSA" = "#7CAE00")


list_plot <- list()
list_age <- c(1,3,5)
for(i in 1:3) {
  
  ggplot(data = pred %>% filter(age %in% list_age[i])) + 
    geom_line(aes(x = c.temp_within, 
                  y = exp(fit),
                  color = source_name,
                  linetype = p_name),
              linewidth = 1) +
    geom_ribbon(aes(x = c.temp_within,
                    ymin = exp(lower),
                    ymax = exp(upper),
                    fill = source_name,
                    linetype = p_name),
                alpha = 0.1) +
    facet_grid(. ~ age_name)  +
    labs(x = expression('T'['individual-anomaly'] * ' (°C)'),
         y = "Predicted increment growth (μm)",
         color = "Dataset",
         fill = "Dataset",
         linetype = "Variable value") +
    scale_color_manual(values = col_scale) +
    scale_fill_manual(values = col_scale) +
    guides(
      colour = guide_legend(order = 1),
      fill = guide_legend(order = 1),
      linedtype = guide_legend(order = 2)
    )
  
  # add facet variable name at the last panel
  if(i == 3) {
    ggplot(data = pred %>% filter(age %in% list_age[i])) + 
      geom_line(aes(x = c.temp_within, 
                    y = exp(fit),
                    color = source_name,
                    linetype = p_name),
                linewidth = 1) +
      geom_ribbon(aes(x = c.temp_within,
                      ymin = exp(lower),
                      ymax = exp(upper),
                      fill = source_name,
                      linetype = p_name),
                  alpha = 0.1) +
      facet_grid(var ~ age_name)  +
      labs(x = expression('T'['individual-anomaly'] * ' (°C)'),
           y = "Predicted increment growth (μm)",
           color = "Dataset",
           fill = "Dataset",
           linetype = "Variable value") +
      scale_color_manual(values = col_scale) +
      scale_fill_manual(values = col_scale) +
      guides(
        colour = guide_legend(order = 1),
        fill = guide_legend(order = 1),
        linedtype = guide_legend(order = 2)
      )
  }
  
  list_plot[[i]] <- last_plot()
}

p4 <- list_plot[[1]]
p5 <- list_plot[[2]]
p6 <- list_plot[[3]]

### f ----
# f - oras5
# range: p25, 50, 75
df_age <- tibble(age = c(1,3,5),
                 age_name = c("Age 1", "Age 3", "Age 5"))

model <- m5_oras5
data <- model@frame

df_f <- tibble(p_name = c("P25", "P75"),
               f = round(quantile(data$f, probs = c(0.25, 0.75)),3)) 

pred <- as.data.frame(Effect(c("log.age", "c.temp_within", "f"), 
                             model, 
                             xlevels = list(log.age = unique(data$log.age),
                                            c.temp_within = seq(-1, 1, 0.1),
                                            f = df_f$f))) %>%
  mutate(age = round(exp(log.age),1)) %>%
  mutate(source = "oras5")

pred <- pred %>% 
  left_join(df_age) %>%
  left_join(df_temp_name) %>%
  left_join(df_f)

#### plot
# add a point of ISIMIP, NEMO-MEDUSA to pred to merge legend 
pred_add <- pred %>% 
  filter(age %in% c(1,3,5), 
         c.temp_within == 0) %>%
  mutate(source_name = c("ISIMIP"))

pred <- bind_rows(pred, pred_add) %>%
  mutate(source_name = factor(source_name, levels = c("ISIMIP", "ORAS5", "NEMO-MEDUSA"))) %>%
  mutate(var = "Fishing mortality")

col_scale <- c("ISIMIP" = "#F8766D", "ORAS5" = "#00BFC4", "NEMO-MEDUSA" = "#7CAE00")

list_plot <- list()
list_age <- c(1,3,5)
for(i in 1:3) {
  
  ggplot(data = pred %>% filter(age %in% list_age[i])) + 
    geom_line(aes(x = c.temp_within, 
                  y = exp(fit),
                  color = source_name,
                  linetype = p_name),
              linewidth = 1) +
    geom_ribbon(aes(x = c.temp_within,
                    ymin = exp(lower),
                    ymax = exp(upper),
                    fill = source_name,
                    linetype = p_name),
                alpha = 0.1) +
    facet_grid(. ~ age_name)  +
    labs(x = expression('T'['individual-anomaly'] * ' (°C)'),
         y = "Predicted increment growth (μm)",
         color = "Dataset",
         fill = "Dataset",
         linetype = "Variable value") +
    scale_color_manual(values = col_scale) +
    scale_fill_manual(values = col_scale) +
    guides(
      colour = guide_legend(order = 1),
      fill = guide_legend(order = 1),
      linedtype = guide_legend(order = 2)
    )
  
  # add facet variable name at the last panel
  if(i == 3) {
    ggplot(data = pred %>% filter(age %in% list_age[i])) + 
      geom_line(aes(x = c.temp_within, 
                    y = exp(fit),
                    color = source_name,
                    linetype = p_name),
                linewidth = 1) +
      geom_ribbon(aes(x = c.temp_within,
                      ymin = exp(lower),
                      ymax = exp(upper),
                      fill = source_name,
                      linetype = p_name),
                  alpha = 0.1) +
      facet_grid(var ~ age_name)  +
      labs(x = expression('T'['individual-anomaly'] * ' (°C)'),
           y = "Predicted increment growth (μm)",
           color = "Dataset",
           fill = "Dataset",
           linetype = "Variable value") +
      scale_color_manual(values = col_scale) +
      scale_fill_manual(values = col_scale) +
      guides(
        colour = guide_legend(order = 1),
        fill = guide_legend(order = 1),
        linedtype = guide_legend(order = 2)
      )
  }
  
  list_plot[[i]] <- last_plot()
}

p7 <- list_plot[[1]]
p8 <- list_plot[[2]]
p9 <- list_plot[[3]]

### c.temp_between ----
# f - oras5
# range: p25, 50, 75
df_age <- tibble(age = c(1, 3, 5),
                 age_name = c("Age 1", "Age 3", "Age 5"))

list_model <- list("oras5" = m5_oras5)

pred <- tibble()
for (i in 1:length(list_model)) {
  model <- list_model[[i]]
  data <- model@frame
  temp_source <- names(list_model[i])
  
  df_c.temp_between <- tibble(p_name = c("P25", "P75"),
                              c.temp_between = round(quantile(data$c.temp_between, probs = c(0.25, 0.75)),3)) 
  
  pred_temp <- as.data.frame(Effect(c("log.age", "c.temp_within", "c.temp_between"), 
                                    model, 
                                    xlevels = list(log.age = unique(data$log.age),
                                                   c.temp_within = seq(-1, 1, 0.1),
                                                   c.temp_between = df_c.temp_between$c.temp_between))) %>%
    mutate(age = round(exp(log.age),1)) %>%
    mutate(source = temp_source) %>%
    left_join(df_c.temp_between)
  
  pred <- bind_rows(pred, pred_temp)
}

pred <- pred %>% 
  left_join(df_age) %>%
  left_join(df_temp_name) 
pred <- pred %>% 
  mutate(var = "T_ind_average")

#### plot
# add a point of ISIMIP, NEMO-MEDUSA to pred to merge legend 
pred_add <- pred %>% 
  filter(age %in% c(1,3,5), 
         c.temp_within == 0) %>%
  mutate(source_name = c("ISIMIP"))

pred <- bind_rows(pred, pred_add) %>%
  mutate(source_name = factor(source_name, levels = c("ISIMIP", "ORAS5", "NEMO-MEDUSA"))) %>%
  mutate(var = "T_ind_average")

col_scale <- c("ISIMIP" = "#F8766D", "ORAS5" = "#00BFC4", "NEMO-MEDUSA" = "#7CAE00")

list_plot <- list()
list_age <- c(1,3,5)
for(i in 1:3) {
  
  ggplot(data = pred %>% filter(age %in% list_age[i])) + 
    geom_line(aes(x = c.temp_within, 
                  y = exp(fit),
                  color = source_name,
                  linetype = p_name),
              linewidth = 1) +
    geom_ribbon(aes(x = c.temp_within,
                    ymin = exp(lower),
                    ymax = exp(upper),
                    fill = source_name,
                    linetype = p_name),
                alpha = 0.1) +
    facet_grid(. ~ age_name)  +
    labs(x = expression('T'['individual-anomaly'] * ' (°C)'),
         y = "Predicted increment growth (μm)",
         color = "Dataset",
         fill = "Dataset",
         linetype = "Variable value") +
    scale_color_manual(values = col_scale) +
    scale_fill_manual(values = col_scale) +
    guides(
      colour = guide_legend(order = 1),
      fill = guide_legend(order = 1),
      linedtype = guide_legend(order = 2)
    )
  
  # add facet variable name at the last panel
  if(i == 3) {
    # add label as expression for geom_grid label
    pred <- pred %>% 
      mutate(age_name2 = factor("age5", labels = expression('Age' * ' 5'))) %>%
      mutate(var2 = factor("T_ind_average",
                           labels = expression('T'['individual-average']))) 
    
    ggplot(data = pred %>% filter(age %in% list_age[i])) + 
      geom_line(aes(x = c.temp_within,
                    y = exp(fit),
                    color = source_name,
                    linetype = p_name),
                linewidth = 1) +
      geom_ribbon(aes(x = c.temp_within,
                      ymin = exp(lower),
                      ymax = exp(upper),
                      fill = source_name,
                      linetype = p_name),
                  alpha = 0.1) +
      facet_grid(var2 ~ age_name2, labeller = label_parsed) +
      labs(x = expression('T'['individual-anomaly'] * ' (°C)'),
           y = "Predicted increment growth (μm)",
           color = "Dataset",
           fill = "Dataset",
           linetype = "Variable value") +
      scale_color_manual(values = col_scale) +
      scale_fill_manual(values = col_scale) +
      guides(
        colour = guide_legend(order = 1),
        fill = guide_legend(order = 1),
        linedtype = guide_legend(order = 2)
      )
  }
  
  list_plot[[i]] <- last_plot()
}

p10 <- list_plot[[1]]
p11 <- list_plot[[2]]
p12 <- list_plot[[3]]

### merge and save plot ----
# merge plot
(p1 + p2 + p3 +
   p4 + p5 + p6 +
   p7 + p8 + p9 +
   p10 + p11 + p12) +
  plot_layout(nrow = 4, axis_titles = "collect", guides = "collect") +
  plot_annotation(tag_levels = 'A') &  
  theme(legend.position = "bottom", 
        legend.box = "vertical",
        axis.title = element_text(size = 9),
        plot.tag.position  = c(0.81, 0.82),
        #plot.tag = element_text(size = 9),
        legend.key.size = unit(15, "pt"),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        strip.text.y = element_text(size = 7))

# save file
ggsave(last_plot(), file = file.path(dir_report, "figS23_temp vs env.pdf"),
       device = cairo_pdf,
       width = 19, height = 25,
       units = "cm")
ggsave(last_plot(), file = file.path(dir_report, "figS23_temp vs env.png"),
       width = 19, height = 25,
       units = "cm",
       dpi = 1000)

## fig S24 - ratio variance ----
#### data 
list_source <- c("isimip", 
  "oras5" 
  #"nemo-medusa"
)

list_pair <- tibble(pop_pair = c("4bc/8ab", "7a/8ab", "4bc/7a"),
                    pop_pair_name = factor(c("North Sea/Bay of Biscay", 
                                             "Irish Sea/Bay of Biscay", 
                                             "North Sea/Irish Sea"),
                                           levels = c("North Sea/Bay of Biscay", 
                                                      "Irish Sea/Bay of Biscay", 
                                                      "North Sea/Irish Sea")))

var_ratio <- tibble()
for (i in 1:length(list_source)) {
  var_ratio_temp <- read_rds(file.path(dir_output, paste0("ind.plastic_var.ratio_", list_source[i], ".rds")))
  
  var_ratio_temp <- var_ratio_temp %>%
    group_by(pop_pair, age_range, n_pop1, n_pop2) %>%
    summarise(across(where(is.numeric), ~ mean(.x))) %>% 
    mutate(across(var_ratio:ci_upper, ~ round(., 2))) %>%
    mutate(p.value = round(p.value, 3))
  
  var_ratio_temp <- var_ratio_temp %>% 
    mutate(source = list_source[i]) %>%
    mutate(sig_level = if_else(p.value < 0.05, "significant", "non-significant")) %>%
    mutate(subset = paste(substr(age_range, 5, 5),
                          substr(age_range, 13, 14),
                          sep = "-"))
  
  var_ratio <- bind_rows(var_ratio, var_ratio_temp)
  
}

var_ratio <- var_ratio %>% 
  left_join(df_temp_name) %>%
  left_join(list_pair) 

#### plot
ggplot() +
  geom_point(data = var_ratio, 
             aes(x = subset,
                 y = var_ratio),
             position = position_dodge(width=0.5)) +
  geom_linerange(data = var_ratio, 
                 aes(x = subset, 
                     y = var_ratio, 
                     ymin = ci_lower, 
                     ymax = ci_upper),
                 position = position_dodge(width=0.5)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  facet_grid(source_name ~ pop_pair_name) +
  scale_color_brewer(palette = "Paired") +
  labs(x = "Increment measurement range",
       y = "Variance ratio") +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 7),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 9))

#### save plot
ggsave(last_plot(), file = file.path(dir_report, "figS24_variance ratio.pdf"),
       device = cairo_pdf,
       width = 18, height = 12.6,
       units = "cm")
ggsave(last_plot(), file = file.path(dir_report, "figS24_variance ratio.png"),
       width = 18, height = 12.6,
       units = "cm",
       dpi = 1000)

## fig S25 - variance individual plasticity ----
#### setup
cor_isimip <- read_rds(file.path(dir_output, paste0("ind.plastic_cor.test_", "isimip", ".rds"))) %>%
  mutate(source = "isimip")
cor_oras5 <- read_rds(file.path(dir_output, paste0("ind.plastic_cor.test_", "oras5", ".rds"))) %>%
  mutate(source = "oras5")
#cor_nm <- read_rds(file.path(dir_output, paste0("ind.plastic_cor.test_", "nemo-medusa", ".rds"))) %>%
#  mutate(source = "nemo-medusa")

list_var <- tibble(var2 = c("c.temp", "ssb.i", "recruitment.i", "f"),
                   var_name = factor(c("Temperature",
                                       "Spawning Stock Biomass",
                                       "Recruitment",
                                       "Fishing mortality"),
                                     levels = c("Temperature",
                                                "Spawning Stock Biomass",
                                                "Recruitment",
                                                "Fishing mortality")))

list_subset <- tibble(subset = c("all ages", "age 6+"),
                      subset_name = c("full range",
                                      "6+"))

cor_all <- bind_rows(cor_isimip, 
                     cor_oras5 
                     #cor_nm
                     )
cor_all <- cor_all %>% 
  left_join(df_temp_name) %>%
  left_join(list_var) %>%
  left_join(list_subset) %>%
  left_join(df_pop)

#### plot
ggplot() +
  geom_point(data = cor_all %>% filter(sig_level == "significant"), 
             aes(x = pop.name, y = cor, shape = subset_name),
             position=position_dodge(width=0.5)) +
  geom_linerange(data = cor_all %>% filter(sig_level == "significant"), 
                 aes(x = pop.name, y = cor, ymin = conf.low, ymax = conf.high, 
                     group = subset_name),
                 position=position_dodge(width=0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(source_name ~ var_name) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom") +
  scale_shape_manual(values = c(1,16)) +
  #scale_color_brewer(palette = "Paired") +
  labs(x = NULL,
       y = "Pearson correlation",
       shape = "Increment measurement range") +
  theme(strip.text = element_text(size = 9),
        #axis.text = element_text(size = 7),
        axis.title = element_text(size = 9),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9))

#### save plot
ggsave(last_plot(), file = file.path(dir_report, "figS25_variance individual plasticity.pdf"),
       device = cairo_pdf,
       width = 19, height = 12.6,
       units = "cm")
ggsave(last_plot(), file = file.path(dir_report, "figS25_variance individual plasticity.png"),
       width = 19, height = 12.6,
       units = "cm", 
       dpi = 1000)

# 4. TABLE ----
## table S1 - otolith sampling - age reading precision ----

#### setup
otl_reage <- read_rds(file.path(dir_otl, "otl_ilvo_reage.rds"))
otl_ilvo <- read_rds(file.path(dir_otl, "otl_ilvo.rds")) #data used for analysis
otl_ilvo <- otl_ilvo %>% filter(FishID != "SOL_52_G1_Z.121_26-03-1990_1602") #remove 1 fish 8ab cohort < 1985

otl_reage <- otl_reage %>% 
  filter(FishID %in% otl_ilvo$FishID) %>%
  filter(FishID != "SOL_52_G1_Z.121_26-03-1990_1602") %>%
  select(FishID, Reader, AgeAtCapture) %>% 
  unique() %>%
  pivot_wider(names_from = Reader, values_from = AgeAtCapture)

#### info
# number of otoliths re-age
sum(!is.na(otl_reage$imaertens))/n_distinct(otl_reage$FishID)*100 #99.5% otoliths re-age

# aging precision tbui vs imaertens
otl_tbui <- otl_reage %>% 
  select(-kdiaz) %>%
  filter(!is.na(tbui))

ap_tbui <- agePrecision(~ tbui + imaertens, data = otl_tbui)
#summary(ap_tbui, what="absolute difference")
summary(ap_tbui, what="precision")

# aging precision kdiaz vs imaertens
otl_kdiaz <- otl_reage %>% 
  select(-tbui) %>%
  filter(!is.na(kdiaz))

ap_kdiaz <- agePrecision(~ kdiaz + imaertens, data = otl_kdiaz)
#summary(ap_kdiaz,what="absolute difference")
summary(ap_kdiaz, what="precision")

## table S2 - temperature - CTD vs temp (monthly average) (Appendix S2) ---- 
# see si_sbt_in-situ ctd vs modeled temperature.Rmd

## table S3 - temperature - annual vs seasonal temp ----
#### setup
temp_season <- bind_rows(isimip, oras, nm) %>%
  filter(pop %in% c("4bc", "7a", "8ab"),
         year >= min(data_otl$year), year <= max(data_otl$year)) %>%
  mutate(month = month(date),
         season = if_else(month %in% c(12,1,2), "Winter (Dec-Feb)",
                          if_else(month %in% c(3,4,5), "Spring (Mar-May)",
                                  if_else(month %in% c(6,7,8), "Summer (Jun-Aug)", "Autumn (Sep-Nov)")))) %>%
  mutate(season = factor(season, levels = c("Summer (Jun-Aug)", "Autumn (Sep-Nov)", "Winter (Dec-Feb)", "Spring (Mar-May)")))

# test by source, pop, season
list_source <- unique(temp_season$source)
list_pop <- unique(temp_season$pop)
list_season <- unique(temp_season$season)

#### table
df <- tibble()
for(i in 1:3) {
  
  for (p in 1:3) {
    
    for (s in 1:4) {
      
      df_year <- temp_season %>%  
        filter(source == list_source[i], pop == list_pop[p]) %>%
        group_by(source, pop, year) %>%
        summarize(temp = mean(temp))
      
      df_season <- temp_season %>%  
        filter(source == list_source[i], pop == list_pop[p]) %>%
        filter(season == list_season[s]) %>%
        group_by(source, pop, year) %>%
        summarize(temp = mean(temp))
      
      lm <- lm(df_year$temp ~ df_season$temp)
      
      df_temp <- tibble(r_square = c(summary(lm)$adj.r.squared),
                        p_value = c(summary(lm)$coefficients[2,4]),
                        season = list_season[s],
                        pop = list_pop[p],
                        source = list_source[i]) %>%
        mutate(r_square = round(r_square, 2),
               p_value = round(p_value, 3))
      
      df <- bind_rows(df, df_temp)
    }
    
  }
  
}

df <- df %>% 
  arrange(season) %>%
  left_join(df_pop) %>%
  left_join(df_temp_name) %>%
  select(season, source_name, pop.name, r_square, p_value) %>%
  rename(Season = season,
         `Temperature dataset` = source_name,
         `Population` = pop.name,
         `R2` = `r_square`,
         `p-value` = p_value) 

#### save table
col.header <- c("Season", 
                "Temperature dataset", 
                "Population", 
                "R\u00b2",
                "p-value")

tab_df(df,
       CSS = css_list,
       col.header = col.header,
       file = file.path(dir_report, "tableS3_temp season.html"))

## table S4, S5 - ORAS5 vs ISIMIP mean and var  ----
data_temp_oras5 <- data_temp %>% filter(source == "oras5")
data_temp_isimip <- data_temp %>% filter(source == "isimip")

### mean ----
# t-test
test_4bc <- t.test(filter(data_temp_oras5, pop == "4bc")$temp, 
                   filter(data_temp_isimip, pop == "4bc")$temp)
test_7a <- t.test(filter(data_temp_oras5, pop == "7a")$temp, 
                  filter(data_temp_isimip, pop == "7a")$temp)
test_8ab <- t.test(filter(data_temp_oras5, pop == "8ab")$temp, 
                   filter(data_temp_isimip, pop == "8ab")$temp)
test_all <- c(test_4bc$p.value, test_7a$p.value, test_8ab$p.value)

# table
df_mean <- data_temp %>%
  filter(source %in% c("oras5", "isimip")) %>%
  group_by(source, pop.name) %>% 
  summarize(mean = sprintf("%.2f", mean(temp))) %>%
  pivot_wider(names_from = source, values_from = c(mean)) %>%
  select(pop.name, oras5, isimip) %>%
  rename(Population = pop.name,
         ORAS5 = oras5,
         ISIMIP = isimip) %>%
  mutate(`t-test p-value` = sprintf("%.2f", test_all))
df_mean

tab_df(df_mean,
       file = file.path(dir_report, "tableS4_oras5 vs isimip_mean.html"))

### variance ----
# variance-test
test_4bc <- var.test(filter(data_temp_oras5, pop == "4bc")$temp, 
                     filter(data_temp_isimip, pop == "4bc")$temp)
test_7a <- var.test(filter(data_temp_oras5, pop == "7a")$temp, 
                    filter(data_temp_isimip, pop == "7a")$temp)
test_8ab <- var.test(filter(data_temp_oras5, pop == "8ab")$temp, 
                     filter(data_temp_isimip, pop == "8ab")$temp)
test_all <- c(test_4bc$p.value, test_7a$p.value, test_8ab$p.value)

# table
df_var <- data_temp %>%
  filter(source %in% c("oras5", "isimip")) %>%
  group_by(source, pop.name) %>% 
  summarize(var = sprintf("%.2f", var(temp))) %>%
  pivot_wider(names_from = source, values_from = c(var)) %>%
  select(pop.name, oras5, isimip) %>%
  rename(Population = pop.name,
         ORAS5 = oras5,
         ISIMIP = isimip) %>%
  mutate(`variation ratio test p-value` = sprintf("%.2f", test_all))
df_var

tab_df(df_var,
       file = file.path(dir_report, "tableS4_oras5 vs isimip_var.html"))

## table S6 - model - random structure ----
#### setup
data <- data_otl

# refit all models
m1a <- lmer(log.increment ~ 1 + log.age + log.aac + log.age*pop + 
              (1 | fishid), data = data, REML = T)
m2a <- lmer(log.increment ~ 1 + log.age + log.aac + log.age*pop +  
              (1 | fishid) + (1 | pop.year), data = data, REML = T)
m2b <- lmer(log.increment ~ 1 + log.age + log.aac + log.age*pop + 
              (1 | fishid) + (1 + log.age | pop.year), data = data, REML = T)
m2e <- lmer(log.increment ~ 1 + log.age + log.aac + log.age*pop +  
              (1 | fishid) + (1 | pop.year) + (1 | pop.cohort), data = data, REML = T)
m2f <- lmer(log.increment ~ 1 + log.age + log.aac + log.age*pop +  
              (1 | fishid) + (1 + log.age | pop.year) + (1 | pop.cohort), data = data, REML = T)
m2g <- lmer(log.increment ~ 1 + log.age + log.aac + log.age*pop +  
              (1 | fishid) + (1 | pop.year) + (1 + log.age | pop.cohort), data = data, REML = T)

# Model comparison
models <- list(m1a, m2a, m2b, m2e, m2f, m2g)
Modnames <- c('m1a', 'm2a', 'm2b', 'm2e', 'm2f', 'm2g')
m_compare <- as.data.frame(aictab(cand.set = models, modnames = Modnames, sort = TRUE))

#### table
df <- tibble(Modnames = c("m1a", "m2a", "m2b", "m2e", "m2f", "m2g"),
             FishID = c(rep("+", 6)),
             `Population:Year` = c("", "+", "+", "+", "+", "+"),
             `Age|Population:Year` = c("", "", "+", "+", "+", ""),
             `Population:Cohort` = c("", "", "", "+", "+", "+"),
             `Age|Population:Cohort` = c(rep("", 5), "+")) %>%
  left_join(m_compare) %>%
  select(-K, -(ModelLik:Cum.Wt)) %>%
  select(- Modnames) %>%
  arrange(AICc) %>%
  mutate(AICc = round(AICc, 2),
         Delta_AICc = round(Delta_AICc, 2)) 

#### save table
col.header <- c("FishID", 
                "Population:Year", 
                "ln(Age)|Population:Year", 
                "Population:Cohort", 
                "ln(Age)|Population:Cohort", 
                "AICc",
                "\u0394AIC")

tab_df(df,
       CSS = css_list,
       col.header = col.header,
       file = file.path(dir_report, "tableS6_random structure.html"))

## table S7 - model - intrinsic structure ----
#### setup
m_compare <- read_rds(file.path(dir_output, "intrinsic.model_comparison.rds"))

df_aic <- as.data.frame(m_compare) %>% 
  select(AICc, delta, id, var_added) %>%
  round(2)

df_var <- as.data.frame(m_compare) %>% 
  select(1:5) %>%
  mutate_all(fun_plus) 

#### table
df <- bind_cols(df_var, df_aic) %>%
  select(`(Intercept)`, log.age, log.aac, pop, `log.age:pop`, AICc, delta, id, var_added)

#### save table
names(df)
col.header <- c("Intercept", 
                 "ln(Age)", 
                 "ln(Age at capture)", 
                 "Population",
                 "ln(Age) * Population", 
                 "AIC", 
                 "\u0394AIC", 
                "ModelID",
                "No. variable added")
tab_df(df,
       CSS = css_list,
       col.header = col.header,
       file = file.path(dir_report, "tableS7_intrinsic structure.html")) 

## table S8 - model - extrinsic structure ----
#### setup
list_model <- tibble(source = c("oras5", "isimip"),
                     id_best = c(65, 1)) # best model id see 2_analyze data.R
df <- tibble()
for(i in 1:2) {
  m_compare <- read_rds(file.path(dir_output, paste0("extrinsic.model_comparison_", list_model$source[i], ".rds")))
  
  df_aic <- as.data.frame(m_compare) %>% 
    select(AICc, delta) %>%
    round(2)
  
  df_info <- as.data.frame(m_compare) %>% 
    select(source, year_range, var_added, id)
  
  df_var <- as.data.frame(m_compare) %>% 
    select(-(df:id)) %>% 
    mutate_all(fun_plus) 
  
  if(list_model$id_best[i] == 1) {
    df_temp <- bind_cols(df_var, df_aic, df_info) %>%
      filter(id %in% c(1, 2) | var_added == min(var_added))
  } else {
    df_temp <- bind_cols(df_var, df_aic, df_info) %>%
      filter(id %in% c(1, list_model$id_best) | var_added == min(var_added))
  }
  
  df <- bind_rows(df, df_temp) 
}

df <- df %>%
  left_join(df_temp_name) %>%
  select(ssb.i, `log.age:ssb.i`, 
         recruitment.i, `log.age:recruitment.i`,
         f, 
         ave.temp, 
         c.temp, `c.temp:log.age`, `c.temp:ssb.i`, `c.temp:recruitment.i`, `c.temp:f`, `ave.temp:c.temp`,
         AICc, delta, id, var_added, source_name)

#### save table
names(df)
col.header <- c("SSB", 
                "SSB * ln(Age)", 
                "Rec", 
                "Rec * ln(Age)", 
                "F", 
                "Tpopulation-average",
                "Tpopulation-anomaly",
                "Tpopulation-anomaly * ln(Age)",
                "Tpopulation-anomaly * SSB",
                "Tpopulation-anomaly * Rec",
                "Tpopulation-anomaly * F",
                "Tpopulation-anomaly * Tpopulation-average",
                "AICc",
                "\u0394AIC",
                "ModelID",
                "No. variable added",
                "Temperature dataset")

tab_df(df,
       CSS = css_list,
       col.header = col.header,
       file = file.path(dir_report, "tableS8_extrinsic structure.html")) 

## table S9 - model - extrinsic extended - nutrient ----
list_model <- tibble(source = c("oras5", "isimip"))
df <- tibble()
for(i in 1:2) {
  m_compare <- read_rds(file.path(dir_output, paste0("extrinsic.model.nu_comparison_", list_model$source[i], ".rds")))
  
  df_aic <- as.data.frame(m_compare) %>% 
    select(AICc, delta) %>%
    round(2)
  
  df_var <- as.data.frame(m_compare) %>% 
    select(-(df:weight)) %>% 
    mutate_all(fun_plus) 
  
  df_temp <- bind_cols(df_var, df_aic) %>%
    mutate(source = list_model$source[i])
  
  df <- bind_rows(df, df_temp) 
}

df <- df %>%
  left_join(df_temp_name) %>%
  select(ssb.i, 
         recruitment.i, 
         c.temp, 
         `c.temp:ssb.i`,
         TN, TP,
         AICc, delta, source_name) %>%
  replace_na(list(c.temp = "", `c.temp:ssb.i` = ""))

#### save table
names(df)
col.header <- c("SSB", 
                "Rec", 
                "Tpopulation-anomaly",
                "Tpopulation-anomaly * SSB",
                "Total nitrogen",
                "Total phosphorus",
                "AICc",
                "\u0394AIC",
                "Temperature dataset")

tab_df(df,
       CSS = css_list,
       col.header = col.header,
       file = file.path(dir_report, "tableS9_extrinsic structure_nutrient.html")) 

## table S10 - model - extrinsic extended structure ----
#### setup
list_model <- tibble(source = c("oras5", "isimip"),
                     id_best = c(1, 1)) # best model id see 2_analyze data.R
df <- tibble()
for(i in 1:2) {
  m_compare <- read_rds(file.path(dir_output, paste0("extrinsic.model.ext_comparison_", list_model$source[i], ".rds")))
  
  df_aic <- as.data.frame(m_compare) %>% 
    select(AICc, delta) %>%
    round(2)
  
  df_info <- as.data.frame(m_compare) %>% 
    select(source, year_range, var_added, id)
  
  df_var <- as.data.frame(m_compare) %>% 
    select(-(df:id)) %>% 
    mutate_all(fun_plus) 
  
  if(list_model$id_best[i] == 1) {
    df_temp <- bind_cols(df_var, df_aic, df_info) %>%
      filter(id %in% c(1, 2) | var_added == min(var_added))
  } else {
    df_temp <- bind_cols(df_var, df_aic, df_info) %>%
      filter(id %in% c(1, list_model$id_best[i]) | var_added == min(var_added))
  }
  
  df <- bind_rows(df, df_temp) 
}

df <- df %>%
  left_join(df_temp_name) %>%
  select(ssb.i, 
         recruitment.i, 
         f, 
         c.temp_within, `c.temp_within:log.age`, `c.temp_within:ssb.i`, `c.temp_within:recruitment.i`, `c.temp_within:f`, `c.temp_between:c.temp_within`, `ave.temp:c.temp_within`,
         c.temp_between, `c.temp_between:log.age`,
         ave.temp,
         AICc, delta, id, var_added, source_name)

#### save table
names(df)
col.header <- c("SSB", 
                "Rec", 
                "F",
                "Tindividual-anomaly",
                "Tindividual-anomaly * ln(Age)",
                "Tindividual-anomaly * SSB",
                "Tindividual-anomaly * Rec",
                "Tindividual-anomaly * F",
                "Tindividual-anomaly * Tindividual-average",
                "Tindividual-anomaly * Tpopulation-average",
                "Tindividual-average",
                "Tindividual-average * ln(Age)",
                "Tpopulation-average",
                "AICc",
                "\u0394AIC",
                "ModelID",
                "No. variable added",
                "Temperature dataset")

tab_df(df,
       CSS = css_list,
       col.header = col.header,
       file = file.path(dir_report, "tableS10_extrinsic extended structure.html")) 

## table S11 - model - extrinsic extended structure - nutrient ----
list_model <- tibble(source = c("oras5", "isimip"))
df <- tibble()
for(i in 1:2) {
  m_compare <- read_rds(file.path(dir_output, paste0("extrinsic.model.ext.nu_comparison_", list_model$source[i], ".rds")))
  
  df_aic <- as.data.frame(m_compare) %>% 
    select(AICc, delta) %>%
    round(2)
  
  df_var <- as.data.frame(m_compare) %>% 
    select(-(df:weight)) %>% 
    mutate_all(fun_plus) 
  
  df_temp <- bind_cols(df_var, df_aic) %>%
    mutate(source = list_model$source[i])
  
  df <- bind_rows(df, df_temp) 
}

df <- df %>%
  left_join(df_temp_name) %>%
  select(ssb.i, 
         recruitment.i, 
         f,
         c.temp_within, `c.temp_within:log.age`, `c.temp_within:ssb.i`, `c.temp_within:recruitment.i`, `c.temp_within:f`, `c.temp_between:c.temp_within`, 
         c.temp_between, `c.temp_between:log.age`,
         TN, TP,
         AICc, delta, source_name) %>%
  replace_na(list(f = "",
                  `c.temp_between:log.age` = "", 
                  `c.temp_within:ssb.i` = "",
                  `c.temp_within:recruitment.i` = "",
                  `c.temp_within:f` = "",
                  `c.temp_between:c.temp_within` = ""
                  ))

#### save table
names(df)
col.header <- c("SSB", 
                "Rec", 
                "F",
                "Tindividual-anomaly",
                "Tindividual-anomaly * ln(Age)",
                "Tindividual-anomaly * SSB",
                "Tindividual-anomaly * Rec",
                "Tindividual-anomaly * F",
                "Tindividual-anomaly * Tindividual-average",
                "Tindividual-average",
                "Tindividual-average * ln(Age)",
                "Total nitrogen",
                "Total phosphorus",
                "AICc",
                "\u0394AIC",
                "Temperature dataset")

tab_df(df,
       CSS = css_list,
       col.header = col.header,
       file = file.path(dir_report, "tableS11_extrinsic structure_nutrient.html")) 

## table S12 - model summary - extrinsic structure ----
pred_labels <- c(
  "Intercept",
  "ln(Age)",
  "ln(Age at capture)",
  "Spawning Stock Biomass",
  "Recruitment",
  "Tpopulation-anomaly",
  "Tpopulation-average",
  "Tpopulation-anomaly * Spawning Stock Biomass"

)

# summary table
tab_model(m4_oras5, m4_isimip, show.se = TRUE, 
                 show.ci = NULL, show.p = TRUE,
          p.threshold = c(0.05),
          p.style = c("stars"),
                 dv.labels = c("Population-level extrinsic model (ORAS5)", 
                               "Population-level extrinsic model (ISIMIP)"),
                 pred.labels = pred_labels,
          order.terms = c(1,2,3,4,5,6,8,7),
                 string.pred = "Fixed Effects",
                 string.est = "Estimate (SE)",
                 collapse.se = TRUE,
                 show.icc = FALSE,
                 digits.rsq = 2,
                 CSS = css_list,
                 file = file.path(dir_report, "tableS12_model_summary.html")
)


## table S13 - model summary - extrinsic extended structure ----
pred_labels <- c(
  "Intercept",
  "ln(Age)",
  "ln(Age at capture)",
  "Tindividual-anomaly",
  "Tindividual-average",
  "Fishing mortality",
  "Spawning Stock Biomass",
  "Recruitment",
  "Tindividual-anomaly * ln(Age)",
  "Tindividual-anomaly * Tindividual-average",
  "Tindividual-anomaly * Fishing mortality",
  "Tindividual-anomaly * Recruitment",
  "Tpopulation-average",
  "Tindividual-average * ln(Age)",
  "Tindividual-anomaly * Spawning Stock Biomass"
)

# summary table
tab_model(m5_oras5, m5_isimip, show.se = TRUE, 
          show.ci = NULL, show.p = TRUE,
          p.threshold = c(0.05),
          p.style = c("stars"),
          dv.labels = c("Individual-level extrinsic model (ORAS5)", 
                        "Individual-level extrinsic model (ISIMIP)"),
          pred.labels = pred_labels,
          order.terms = c(1,2,3,7,8,6,4,5,9,15,12,11,10,14,13),
          string.pred = "Fixed Effects",
          string.est = "Estimate (SE)",
          collapse.se = TRUE,
          show.icc = FALSE,
          digits.rsq = 2,
          CSS = css_list,
          file = file.path(dir_report, "tableS13_model_summary_ext.html")
)


# 5. SI text ----
## preparation method: sectioned vs sectioned/stained -----
data <- data_otl

# sectioned data is only available in 8ab
data_sub <- data %>% filter(pop == "8ab")
data_sub <- data_sub %>% filter(method %in% c('sectioned', 'sectioned/stained'))

# year range
# View(data_sub %>% group_by(method) %>%
#        summarize(year_range = range(year)))
# sectioned: 2002-2016
# sectioned/stained: 1997-2018

# year range with more than 10 increments
# View(data_sub %>% group_by(method, year) %>%
#        summarize(n = n()))
# sectioned: 2006-2014
# sectioned/stained: 2003-2018

# overlapping period: 2006-2014
data_sub2 <- data_sub %>% filter(year >= 2006, year <= 2014)
table(data_sub2$method)

# model
m3 <- lmer(log.increment ~ 1 + log.age + log.aac + method + 
             (1 | fishid) + (1 + log.age | pop.year) + (1 | pop.cohort), data = data_sub2, 
           REML = T, na.action = "na.fail")
summary(m3)
# note: no difference between preparation method "sectioned" and "sectioned/stained" reading for the same period (1990-1995) and the same preparation method (sectioned/stained)

## broken and sectioned same period ----
data <- data_otl %>% 
  mutate(method2 = if_else(method == "broken/burned", "broken", "sectioned"))

data %>% group_by(pop, method2) %>% summarize(year_range = range(year))

### 4bc ----
data_sub <- data %>%
  filter(pop == "4bc", 
         datasource == "ILVO") 
data_sub %>%
  group_by(pop, method2) %>% summarize(year_range = range(year)) 
# broken: 1960-1998
# sectioned: 1988-2019

# year range with more than 10 increments
# View(data_sub %>% 
#        group_by(method2, year) %>%
#        summarize(n = n()))
# broken: 1969-1998
# sectioned: 1990-2019

# overlapping period: 1990-1998
data_sub2 <- data_sub %>% filter(year >= 1990, year <= 1998)
table(data_sub2$method2)
# model
m3 <- lmer(log.increment ~ 1 + log.age + log.aac + method2 + 
             (1 | fishid) + (1 + log.age | pop.year) + (1 | pop.cohort), data = data_sub, 
           REML = T, na.action = "na.fail")
summary(m3)

### 7a ----
data_sub <- data %>%
  filter(pop == "7a", 
         datasource == "ILVO") 
data_sub %>%
  group_by(pop, method2) %>% summarize(year_range = range(year)) 
# broken: 1971-1998
# sectioned: 1984-2019

# year range with more than 10 increments
# View(data_sub %>% 
#        group_by(method2, year) %>%
#        summarize(n = n()))
# broken: 1971-1997
# sectioned: 1993-2019

# overlapping period: 1993-1997
data_sub2 <- data_sub %>% filter(year >= 1993, year <= 1997)
table(data_sub2$method2)
# model
m3 <- lmer(log.increment ~ 1 + log.age + log.aac + method2 + 
             (1 | fishid) + (1 + log.age | pop.year) + (1 | pop.cohort), data = data_sub, 
           REML = T, na.action = "na.fail")
summary(m3)

### 8ab ----
data_sub <- data %>%
  filter(pop == "8ab", 
         datasource == "ILVO") 
data_sub %>%
  group_by(pop, method2) %>% summarize(year_range = range(year)) 
# broken: 1990-2003
# sectioned: 1997-2018

# year range with more than 10 increments
View(data_sub %>% 
       group_by(method2, year) %>%
       summarize(n = n()))
# broken: 1971-2000
# sectioned: 2002-2018

# no overlapping period

## reading institute (datasource): WUR vs ILVO -----
data <- data_otl

# different reading institutes (datasource) are only available in 4bc
data_sub <- data %>% filter(pop == "4bc")
data_sub <- data_sub %>% filter(datasource %in% c('WUR', 'ILVO'))

# year range
# View(data_sub %>% group_by(datasource, method) %>%
#        summarize(year_range = range(year)))

# year range with more than 10 increments
# View(data_sub %>% 
#        filter(method == "sectioned/stained") %>%
#        group_by(datasource, year) %>%
#        summarize(n = n()))
# ILVO: 1990-2019
# WUR: 1958-1995

# overlapping period: 1990-1995
data_sub2 <- data_sub %>% filter(method == "sectioned/stained", year >= 1990, year <= 1995)
table(data_sub2$datasource)

# model
m3 <- lmer(log.increment ~ 1 + log.age + log.aac + datasource + 
             (1 | fishid) + (1 + log.age | pop.year) + (1 | pop.cohort), data = data_sub2, 
           REML = T, na.action = "na.fail")
# singular model with (1 + log.age | pop.year) -> refit with (1 | log.year)
m3 <- lmer(log.increment ~ 1 + log.age + log.aac + datasource + 
             (1 | fishid) + (1 | pop.year) + (1 | pop.cohort), data = data_sub2, 
           REML = T, na.action = "na.fail")
summary(m3)

# note: no difference between WUR and ILVO reading for the same period (1990-1995) and the same preparation method (sectioned/stained)

# 6. Other ----
## Fishing mortality ----
### new (same age group 3-7 across pops) vs old (diff age groups across pops) Fbar 
data_sol_raw <- read_rds(file.path(dir_ices, "stock-assessment_2023.rds")) %>%
  rename(ssb = SSB,
         f = `F`) %>%
  left_join(df_pop)

data_sol_fbar <- read_rds(file.path(dir_ices, "sol_fbar_age3-7_stock-assessment_2023.rds")) %>%
  left_join(df_pop)

ggplot() +
  geom_line(data = data_sol_raw, aes(x = year, y = f, linetype = "old Fbar")) +
  geom_line(data = data_sol_fbar, aes(x = year, y = f, linetype = "new Fbar")) +
  facet_grid(pop.name ~ .) +
  labs(x = "Year",
       y = "Fishing mortality",
       linetype = NULL) +
  theme(legend.position = "bottom")

### fbar 3-7 vs age at capture 3-20+
library(tidyverse)
library(readxl)
dir_ices <- "./data/ices"

# f-at-age range vs stock assessment Fbar: 
## 4: 1-10 (Fbar 2-6)
## 7a: 2-8 (Fbar 4-7)
## 8ab: 2-8 (Fbar 3-6)
# -> new Fbar 3-7

f4 <- read_xlsx(file.path(dir_ices, "sol_4_f-at-age.xlsx")) 
f7a <- read_xlsx(file.path(dir_ices, "sol_7a_f-at-age.xlsx")) 
f8ab <- read_xlsx(file.path(dir_ices, "sol_8ab_f-at-age.xlsx")) 

f4 <- f4 %>% 
  pivot_longer(cols = age1:age10, names_to = "age", values_to = "f") %>%
  mutate(age = as.integer(substr(age, 4, 5))) %>%
  mutate(pop = "4bc")

f7a <- f7a %>% 
  pivot_longer(cols = age2:age8, names_to = "age", values_to = "f") %>%
  mutate(age = as.integer(substr(age, 4, 4))) %>%
  mutate(pop = "7a")

f8ab <- f8ab %>% 
  pivot_longer(cols = age2:age8, names_to = "age", values_to = "f") %>%
  mutate(age = as.integer(substr(age, 4, 4))) %>%
  mutate(pop = "8ab")

# combine all f together
data_sol_f <- bind_rows(f4, f7a, f8ab) %>%
  left_join(df_pop)

ggplot(data = data_sol_f, aes(x = age, y = f)) +
  geom_rect(aes(xmin = 3, xmax = 7, ymin = -Inf, ymax = Inf), fill = "grey", alpha = 0.01) +
  geom_point(alpha = 0.1) +
  facet_grid(pop.name ~ .) +
  scale_x_continuous(breaks = seq(1,10)) +
  labs(x = "Age",
       y = "Fishing mortality at age")
p1 <- last_plot()

data <- data_otl %>%
  left_join(df_pop)
ggplot(data = data, aes(x = AgeAtCapture, y = Length.mm)) +
  geom_point(alpha = 0.1) +
  #geom_hline(yintercept = 240, linetype = "dashed") +
  facet_grid(pop.name ~ .) +
  labs(x = "Age at capture",
       y = "Fish length (mm)")
p2 <- last_plot()

p1 | p2

## Random effect age specific ----
# year random intercept + log.age slope (1 + log.age | pop.year) 

#### data 
data <- data_otl
df_age <- unique(select(data, age, log.age)) %>% arrange(age)

pred_year <- ranef(m3)$`pop.year` 
pred_year_se <- as.data.frame(arm::se.ranef(m3)$`pop.year`)

pred_year <- pred_year %>%
  mutate(pop = sub(":.*", "", rownames(pred_year)),
         year = as.numeric(sub(".*:", "", rownames(pred_year))),
         intercept = `(Intercept)`,
         intercept_se = pred_year_se$`(Intercept)`,
         log.age_se =  pred_year_se$log.age) %>%
  left_join(df_pop)

# random effect by age
pred_year2 <- pred_year %>%
  mutate(`1` = intercept + log.age*df_age$log.age[1],
         `2` = intercept + log.age*df_age$log.age[2],
         `3` = intercept + log.age*df_age$log.age[3],
         `4` = intercept + log.age*df_age$log.age[4],
         `5` = intercept + log.age*df_age$log.age[5],
         `6` = intercept + log.age*df_age$log.age[6],
         `7` = intercept + log.age*df_age$log.age[7],
         `8` = intercept + log.age*df_age$log.age[8],
         `9` = intercept + log.age*df_age$log.age[9],
         `10` = intercept + log.age*df_age$log.age[10],
         `11` = intercept + log.age*df_age$log.age[11],
         `12` = intercept + log.age*df_age$log.age[12],
         `13` = intercept + log.age*df_age$log.age[13],
         `14` = intercept + log.age*df_age$log.age[14],
         `15` = intercept + log.age*df_age$log.age[15],
         `22` = intercept + log.age*df_age$log.age[22]) %>%
  select(pop.name, year,  `1`:`22`) %>%
  pivot_longer(cols = `1`:`22`, names_to = "age", values_to = "ranef") %>%
  mutate(age = as.numeric(age))

ggplot() +
  geom_line(data = pred_year2, aes(x = log(age), y = ranef, group = year)) +
  geom_point(data = pred_year2 %>% filter(age %in% c(1)), 
             aes(x = log(age), y = ranef, group = year, color = "age 1"), alpha = 0.5) +
  geom_point(data = pred_year2 %>% filter(age %in% c(5)), 
             aes(x = log(age), y = ranef, group = year, color = "age 5"), alpha = 0.5) +
  facet_grid(. ~ pop.name) +
  labs(x = "log(Age)",
       y = "Year random effect",
       color = NULL)


age1 <- pred_year2 %>% filter(age == 1) %>% mutate(age_group = "Age 1")
age5 <- pred_year2 %>% filter(age == 5) %>% mutate(age_group = "Age 5")
age_all <- bind_rows(age1, age5)
ggplot(data = age_all, aes(x = year, y = ranef)) +
  geom_line() + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(pop.name ~ age_group) +
  labs(x = "Year",
       y = "Year random effect")


## average age
pred_year3 <- pred_year %>%
  mutate(`1` = 7.8 + intercept + (log.age - 1.23)*df_age$log.age[1],
         `2` = 7.8 + intercept + (log.age - 1.23)*df_age$log.age[2],
         `3` = 7.8 + intercept + (log.age - 1.23)*df_age$log.age[3],
         `4` = 7.8 + intercept + (log.age - 1.23)*df_age$log.age[4],
         `5` = 7.8 + intercept + (log.age - 1.23)*df_age$log.age[5],
         `6` = 7.8 + intercept + (log.age - 1.23)*df_age$log.age[6],
         `22` = 7.8 + intercept + (log.age - 1.23)*df_age$log.age[22]
         ) %>%
  select(pop.name, year,  `1`:`22`) %>%
  pivot_longer(cols = `1`:`22`, names_to = "age", values_to = "fitted") %>%
  mutate(age = as.numeric(age))

age_effect <- tibble(log.age = df_age$log.age,
                     age = df_age$age,
                     fitted = 7.8 - 1.23*log.age)

ggplot() +
  geom_line(data = pred_year3, aes(x = log(age), y = fitted, group = year), color = "grey") +
  geom_line(data = age_effect, aes(x = log(age), y = fitted)) +
  facet_grid(~ pop.name) +
  labs(x = "log(Age (years))",
       y = "log(Predicted increment growth (μm))")




