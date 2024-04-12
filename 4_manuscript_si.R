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
dir_output <- "./output"
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
df_temp_name <- tibble(source = c("isimip", "oras5", "nemo-medusa"),
                       source_name = factor(c("ISIMIP", "ORAS5", "NEMO-MEDUSA"), 
                                            levels = c("ISIMIP", "ORAS5", "NEMO-MEDUSA")))
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
         recruitment.i = recruitment/area_km2)

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
m4_nm <- read_rds(file.path(dir_output, "extrinsic.model_best_nemo-medusa.rds"))

## scaled
m4_isimip_s <- read_rds(file.path(dir_output, "extrinsic.model_best_scaled_isimip.rds"))
m4_oras5_s <- read_rds(file.path(dir_output, "extrinsic.model_best_scaled_oras5.rds"))
m4_nm_s <- read_rds(file.path(dir_output, "extrinsic.model_best_scaled_nemo-medusa.rds"))

#### BEST EXTRINSIC STRUCTURE WITH NUTRIENT DATA:
## no scaled
m4_nu_isimip <- read_rds(file.path(dir_output, "extrinsic.model.nu_best_isimip.rds"))
m4_nu_oras5 <- read_rds(file.path(dir_output, "extrinsic.model.nu_best_oras5.rds"))
m4_nu_nm <- read_rds(file.path(dir_output, "extrinsic.model.nu_best_nemo-medusa.rds"))

## scaled
m4_nu_isimip_s <- read_rds(file.path(dir_output, "extrinsic.model.nu_best_scaled_isimip.rds"))
m4_nu_oras5_s <- read_rds(file.path(dir_output, "extrinsic.model.nu_best_scaled_oras5.rds"))
m4_nu_nm_s <- read_rds(file.path(dir_output, "extrinsic.model.nu_best_scaled_nemo-medusa.rds"))

#### BEST EXTRINSIC EXTENDED STRUCTURE: 
## no scaled
m5_isimip <- read_rds(file.path(dir_output, "extrinsic.model.ext_best_isimip.rds"))
m5_oras5 <- read_rds(file.path(dir_output, "extrinsic.model.ext_best_oras5.rds"))
m5_nm <- read_rds(file.path(dir_output, "extrinsic.model.ext_best_nemo-medusa.rds"))

## scaled
m5_isimip_s <- read_rds(file.path(dir_output, "extrinsic.model.ext_best_scaled_isimip.rds"))
m5_oras5_s <- read_rds(file.path(dir_output, "extrinsic.model.ext_best_scaled_oras5.rds"))
m5_nm_s <- read_rds(file.path(dir_output, "extrinsic.model.ext_best_scaled_nemo-medusa.rds"))

#### BEST EXTRINSIC EXTENDED STRUCTURE WITH NUTRIENT DATA: 
## no scaled
m5_nu_isimip <- read_rds(file.path(dir_output, "extrinsic.model.ext.nu_best_isimip.rds"))
m5_nu_oras5 <- read_rds(file.path(dir_output, "extrinsic.model.ext.nu_best_oras5.rds"))
m5_nu_nm <- read_rds(file.path(dir_output, "extrinsic.model.ext.nu_best_nemo-medusa.rds"))

## scaled
m5_nu_isimip_s <- read_rds(file.path(dir_output, "extrinsic.model.ext.nu_best_scaled_isimip.rds"))
m5_nu_oras5_s <- read_rds(file.path(dir_output, "extrinsic.model.ext.nu_best_scaled_oras5.rds"))
m5_nu_nm_s <- read_rds(file.path(dir_output, "extrinsic.model.ext.nu_best_scaled_nemo-medusa.rds"))

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
       y = "Number of otoliths")

# save plot
ggsave(last_plot(), file = file.path(dir_report, "figS6_age at capture distribution.tiff"),
       width = 17, height = 5.75,
       units = "cm",
       scaling = 0.8,
       dpi = 600)
ggsave(last_plot(), file = file.path(dir_report, "figS6_age at capture distribution.eps"),
       width = 17, height = 5.75,
       units = "cm",
       dpi = 600, 
       device = "eps")

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
        legend.justification = c(0.5, 0.5))

# save plot
ggsave(last_plot(), file = file.path(dir_report, "figS7_sampling size.tiff"),
       width = 17, height = 5.75,
       units = "cm",
       scaling = 0.8,
       dpi = 600)
ggsave(last_plot(), file = file.path(dir_report, "figS7_sampling size.eps"),
       width = 17, height = 5.75,
       units = "cm",
       dpi = 600, 
       device = "eps")

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
  labs(x = "Otolith width (μm)",
       y = "Fish length (mm)") +
  annotate("text", 
           x = 5500, y = 240, 
           label = expression("" ~ R^2 ~ " = 0.53, p < 0.001, n = 2152"), 
           size = 3) 

# note: 
nrow(otl_sub) #2152 
# 2 fish with otolith width very small width (< 1500 um) likely due to error were not included
# R correlation: 0.73
# R2: 0.53

# save plot
ggsave(last_plot(), file = file.path(dir_report, "figS8_otolith vs fish size.tiff"),
       width = 17, height = 11.5,
       units = "cm",
       scaling = 0.8,
       dpi = 600)
ggsave(last_plot(), file = file.path(dir_report, "figS8_otolith vs fish size.eps"),
       width = 17, height = 11.5,
       units = "cm",
       dpi = 600, 
       device = "eps")

## fig S9 - temperature - CTD by year (Appendix S2) ----
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
ggsave(last_plot(), file = file.path(dir_report, "figS9_ctd by year.tiff"),
       width = 17, height = 11.5,
       units = "cm",
       scaling = 0.8,
       dpi = 600)
ggsave(last_plot(), file = file.path(dir_report, "figS9_ctd by year.eps"),
       width = 17, height = 11.5,
       units = "cm",
       dpi = 600, 
       device = "eps")

## fig S10 - temperature - CTD vs temp (Appendix S2) ----
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
ggsave(last_plot(), file = file.path(dir_report, "figS10_ctd vs temp.tiff"),
       width = 17, height = 11.5,
       units = "cm",
       scaling = 0.8,
       dpi = 600)
ggsave(last_plot(), file = file.path(dir_report, "figS10_ctd vs temp.eps"),
       width = 17, height = 11.5,
       units = "cm",
       dpi = 600, 
       device = "eps")

## fig S11 - temperature - temp ----
#### plot
ggplot(data = data_temp %>% left_join(df_pop), 
       aes(x = year, y = temp, color = pop.name)) +
  geom_line() +
  facet_grid(source_name ~ .) +
  labs(x = "Year",
       y = "Temperature (°C)",
       color = "Population") +
  theme(legend.position = c(0.01, 0.2),
        legend.justification = c(0, 1),
        axis.title.x = element_blank(),
        legend.title = element_blank()) +
  scale_color_brewer(palette = "Dark2", direction = -1) 

#### save plot
ggsave(last_plot(), file = file.path(dir_report, "figS11_temperature.tiff"),
       width = 17, height = 11.5,
       units = "cm",
       scaling = 0.8,
       dpi = 600)
ggsave(last_plot(), file = file.path(dir_report, "figS11_temperature.eps"),
       width = 17, height = 11.5,
       units = "cm",
       dpi = 600, 
       device = "eps")

## fig S12 - temperature - c.temp ----
#### plot
ggplot(data = data_temp %>% left_join(df_pop), 
       aes(x = year, y = c.temp, color = source_name)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(pop.name ~ .) +
  labs(x = "Year",
       y = expression('T'['population-anomaly'] * ' (°C)'),
       color = "Dataset") +
  theme(legend.position = "bottom") +
  scale_color_manual(values = col_scale)

#### save plot
ggsave(last_plot(), file = file.path(dir_report, "figS12_temperature pop anomaly.tiff"),
       width = 17, height = 11.5,
       units = "cm",
       scaling = 0.8,
       dpi = 600)
ggsave(last_plot(), file = file.path(dir_report, "figS12_temperature pop anomaly.eps"),
       width = 17, height = 11.5,
       units = "cm",
       dpi = 600, 
       device = "eps")

## fig S13 - distribution area (Appendix S3) ----
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
  theme(legend.position = "bottom")

# save plot
ggsave(last_plot(), file = file.path(dir_report, "figS13_datras survey.tiff"),
       width = 17, height = 11.5,
       units = "cm",
       scaling = 0.8,
       dpi = 600)
ggsave(last_plot(), file = file.path(dir_report, "figS13_datras survey.eps"),
       width = 17, height = 11.5,
       units = "cm",
       dpi = 600, 
       device = "eps")

## fig S14 - variable - spawning stock biomass, recruitment ---- 
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
  scale_color_brewer(palette = "Dark2", direction = -1)  

p2 <- ggplot(data = data_sol, aes(x = year, y = recruitment, color = pop.name)) +
  geom_line() +
  labs(x = "Year",
       y = "Recruitment (thoudsand)",
       color = "Population") +
  theme(legend.position = c(0.6, 1.05),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = "transparent"),
        legend.title = element_blank()) +
  scale_color_brewer(palette = "Dark2", direction = -1)  

p1 + p2 + 
  plot_annotation(tag_levels = 'A') 

#### save plot
ggsave(last_plot(), file = file.path(dir_report, "figS14_raw ssb rec.tiff"),
       width = 17, height = 5.75,
       units = "cm",
       scaling = 0.8,
       dpi = 600)
ggsave(last_plot(), file = file.path(dir_report, "figS14_raw ssb rec.eps"),
       width = 17, height = 5.75,
       units = "cm",
       dpi = 600, 
       device = "eps")

## fig S15 - variable - TN and TP ----
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
  plot_annotation(tag_levels = 'A') 

#### save plot
ggsave(last_plot(), file = file.path(dir_report, "figS15_nutrient.tiff"),
       width = 17, height = 5.75,
       units = "cm",
       scaling = 0.8,
       dpi = 600)
ggsave(last_plot(), file = file.path(dir_report, "figS15_nutrient.eps"),
       width = 17, height = 5.75,
       units = "cm",
       dpi = 600, 
       device = "eps")

## fig S16 - intrinsic - reading institute (datasource) ----
#### setup
pred <- as.data.frame(Effect(c("log.age", "datasource"), 
                             m3, 
                             xlevels = list(log.age = unique(data$log.age)))) %>%
  mutate(age = round(exp(log.age),0)) %>%
  arrange(age)

#### plot
ggplot(data = pred) + 
  geom_line(aes(x = age, 
                y = exp(fit),
                linetype = datasource),
            linewidth = 0.5) +
  geom_linerange(aes(x = age,
                     ymin = exp(lower),
                     ymax = exp(upper),
                     linetype = datasource)) +
  labs(x = "Age (years)",
       y = "Predicted increment growth (μm)",
       color = "Reading institute",
       fill = "Reading institute",
       linetype = "Reading institute") +
  theme(legend.position = "bottom")

#### save plot
ggsave(last_plot(), file = file.path(dir_report, "figS16_age x reading institute.tiff"),
       width = 8.5, height = 6.5,
       units = "cm",
       scaling = 0.8,
       dpi = 600)
ggsave(last_plot(), file = file.path(dir_report, "figS16_age x reading institute.eps"),
       width = 8.5, height = 6.5,
       units = "cm",
       dpi = 600, 
       device = "eps")

## fig S17 - intrinsic - pop ----
#### setup
pred_pop <- ranef(m3)$pop
pred_pop_se <- as.data.frame(arm::se.ranef(m3)$pop)

pred_pop <- pred_pop %>%
  mutate(pop = rownames(pred_pop),
         intercept = `(Intercept)`,
         se = pred_pop_se$`(Intercept)`
  ) %>%
  left_join(df_pop)

#### plot
ggplot(data = pred_pop) + 
  geom_point(aes(x = pop.name, 
                 y = intercept)) +
  geom_linerange(aes(x = pop.name,
                    ymin = intercept + se,
                    ymax = intercept - se)) +
  labs(x = "Population",
       y = "Population random effect")

#### save plot
ggsave(last_plot(), file = file.path(dir_report, "figS17_pop.tiff"),
       width = 8.5, height = 5.75,
       units = "cm",
       scaling = 0.8,
       dpi = 600)
ggsave(last_plot(), file = file.path(dir_report, "figS17_pop.eps"),
       width = 8.5, height = 5.75,
       units = "cm",
       dpi = 600, 
       device = "eps")

## fig S18 - comparing extrinsic effect - without nutrient (scaled model) ----
### extrinsic oras5 ----
#### setup
# df_age
data_age <- data_otl %>%
  mutate(s.log.age = s.(log.age))
df_age <- select(data_age, age, s.log.age) %>% unique() %>% arrange(age) %>% round(2)

# data and model
data <- m4_oras5_s@frame
model <- m4_oras5_s

# temperature
pred_temp <- as.data.frame(Effect(c("s.log.age", "s.c.temp"),
                                  model,
                                  xlevels = list(s.log.age = df_age$s.log.age,
                                                 s.c.temp = seq(0, 1, 0.1)))) %>%
  mutate(var = s.c.temp) %>%
  left_join(df_age) %>%
  mutate(fit_trans = exp(fit),
         upper_trans = exp(upper),
         lower_trans = exp(lower)) %>%
  select(age, var, fit_trans, upper_trans, lower_trans) %>%
  mutate(var_name = "Temperature",
         age_name = paste0("Age ", age))

# ssb
pred_ssb <- as.data.frame(Effect(c("s.log.age", "s.ssb.i"),
                                 model,
                                 xlevels = list(s.log.age = df_age$s.log.age,
                                                s.ssb.i = seq(0, 1, 0.1)))) %>%
  mutate(var = s.ssb.i) %>%
  left_join(df_age) %>%
  mutate(fit_trans = exp(fit),
         upper_trans = exp(upper),
         lower_trans = exp(lower)) %>%
  select(age, var, fit_trans, upper_trans, lower_trans) %>%
  mutate(var_name = "Spawning stock biomass",
         age_name = paste0("Age ", age))


# recruitment
pred_rec <- as.data.frame(Effect(c("s.log.age", "s.recruitment.i"),
                                 model,
                                 xlevels = list(s.log.age = df_age$s.log.age,
                                                s.recruitment.i = seq(0, 1, 0.1)))) %>%
  mutate(var = s.recruitment.i) %>%
  left_join(df_age) %>%
  mutate(fit_trans = exp(fit),
         upper_trans = exp(upper),
         lower_trans = exp(lower)) %>%
  select(age, var, fit_trans, upper_trans, lower_trans) %>%
  mutate(var_name = "Recruitment",
         age_name = paste0("Age ", age))

# all extrinsic variables
pred <- bind_rows(pred_temp, pred_ssb, pred_rec)

#### summary
pred_sum <- pred %>% 
  filter(var %in% c(0,1)) %>%
  select(var_name, age, var, fit_trans) %>%
  group_by(var_name, age) %>%
  pivot_wider(names_from = var, values_from = fit_trans) %>%
  mutate(diff = (`1`-`0`)/`0`*100) %>%
  mutate(across(where(is.numeric), round, 2)) %>%
  arrange(age)

pred_sum <- pred_sum %>%
  mutate(var_name = factor(var_name, levels = c("Temperature", 
                                                "Spawning stock biomass",
                                                "Recruitment")),
         sign = if_else(diff < 0, "Negative", "Positive"),
         source_name = "ORAS5")

#### plot
scale_shape = c("Temperature" = 15, 
                "Spawning stock biomass" = 2, 
                "Recruitment" = 5,
                "Total phosphorus" = 6)

unique(pred_sum$var_name)

ggplot(data = pred_sum) +
  geom_point(aes(x = age, y = abs(diff), shape = var_name, color = sign)) +
  scale_shape_manual(values = scale_shape,
                     labels = expression('Temperature'['population-anomaly'],
                                         'Spawning Stock Biomass', 
                                         'Recruitment')) +
  theme(legend.text.align = 0) +
  labs(x = "Age",
       y = "Growth change (%)",
       shape = "Predictor",
       color = "Effect direction") +
  facet_grid(. ~ source_name)

p1 <- last_plot()   

### extrinsic nemo-medusa ----
#### setup
# df_age
data_age <- data_otl %>%
  mutate(s.log.age = s.(log.age))
df_age <- select(data_age, age, s.log.age) %>% unique() %>% arrange(age) %>% round(2)

# data and model
data <- m4_nm_s@frame
model <- m4_nm_s

# temperature
pred_temp <- as.data.frame(Effect(c("s.log.age", "s.c.temp"),
                                  model,
                                  xlevels = list(s.log.age = df_age$s.log.age,
                                                 s.c.temp = seq(0, 1, 0.1)))) %>%
  mutate(var = s.c.temp) %>%
  left_join(df_age) %>%
  mutate(fit_trans = exp(fit),
         upper_trans = exp(upper),
         lower_trans = exp(lower)) %>%
  select(age, var, fit_trans, upper_trans, lower_trans) %>%
  mutate(var_name = "Temperature",
         age_name = paste0("Age ", age))

# ssb
pred_ssb <- as.data.frame(Effect(c("s.log.age", "s.ssb.i"),
                                 model,
                                 xlevels = list(s.log.age = df_age$s.log.age,
                                                s.ssb.i = seq(0, 1, 0.1)))) %>%
  mutate(var = s.ssb.i) %>%
  left_join(df_age) %>%
  mutate(fit_trans = exp(fit),
         upper_trans = exp(upper),
         lower_trans = exp(lower)) %>%
  select(age, var, fit_trans, upper_trans, lower_trans) %>%
  mutate(var_name = "Spawning stock biomass",
         age_name = paste0("Age ", age))


# recruitment
pred_rec <- as.data.frame(Effect(c("s.log.age", "s.recruitment.i"),
                                 model,
                                 xlevels = list(s.log.age = df_age$s.log.age,
                                                s.recruitment.i = seq(0, 1, 0.1)))) %>%
  mutate(var = s.recruitment.i) %>%
  left_join(df_age) %>%
  mutate(fit_trans = exp(fit),
         upper_trans = exp(upper),
         lower_trans = exp(lower)) %>%
  select(age, var, fit_trans, upper_trans, lower_trans) %>%
  mutate(var_name = "Recruitment",
         age_name = paste0("Age ", age))

# all extrinsic variables
pred <- bind_rows(pred_temp, pred_ssb, pred_rec)

#### summary
pred_sum <- pred %>% 
  filter(var %in% c(0,1)) %>%
  select(var_name, age, var, fit_trans) %>%
  group_by(var_name, age) %>%
  pivot_wider(names_from = var, values_from = fit_trans) %>%
  mutate(diff = (`1`-`0`)/`0`*100) %>%
  mutate(across(where(is.numeric), round, 2)) %>%
  arrange(age)

pred_sum <- pred_sum %>%
  mutate(var_name = factor(var_name, levels = c("Temperature", 
                                                "Spawning stock biomass",
                                                "Recruitment")),
         sign = if_else(diff < 0, "Negative", "Positive"),
         source_name = "NEMO-MEDUSA")

#### plot
scale_shape = c("Temperature" = 15, 
                "Spawning stock biomass" = 2, 
                "Recruitment" = 5,
                "Total phosphorus" = 6)

unique(pred_sum$var_name)

ggplot(data = pred_sum) +
  geom_point(aes(x = age, y = abs(diff), shape = var_name, color = sign)) +
  scale_shape_manual(values = scale_shape,
                     labels = expression('Temperature'['population-anomaly'],
                                         'Spawning Stock Biomass', 
                                         'Recruitment')) +
  theme(legend.text.align = 0) +
  labs(x = "Age",
       y = "Growth change (%)",
       shape = "Predictor",
       color = "Effect direction") +
  facet_grid(. ~ source_name)

p2 <- last_plot()

### merge and save plot ----
p1 + p2 +
  plot_annotation(tag_levels = 'A') +
  plot_layout(guides = "collect")
  

ggsave(last_plot(), file = file.path(dir_report, "figS18_extrinsic comparison without nutrient.tiff"),
       width = 17, height = 5.75,
       units = "cm",  scaling = 0.8,
       dpi = 600)
ggsave(last_plot(), file = file.path(dir_report, "figS18_extrinsic comparison without nutrient.eps"),
       width = 17, height = 5.75,
       units = "cm", 
       dpi = 600, 
       device = "eps")

## fig S19 - comparing extrinsic effect - with nutrient (scaled model) ----
### extrinsic oras5 ----
#### setup
# model
model <- m4_nu_oras5_s

# data and model
data <- data_otl %>% 
  left_join(filter(data_temp, source == "oras5"),
            by = join_by(pop, year)) %>%
  filter(is.na(c.temp) == F) %>%
  left_join(data_sol, by = join_by(pop, year)) %>%
  left_join(data_nu, by = join_by(pop, year))
data <- data %>% filter(!is.na(TN) & !is.na(TP))
data <- data %>%
  mutate(s.log.age = s.(log.age),
         s.log.aac = s.(log.aac),
         s.c.temp = s.(c.temp),
         s.ssb.i = s.(ssb.i),
         s.recruitment.i = s.(recruitment.i),
         s.TP = s.(TP))

# df_age
df_age <- select(data, age, s.log.age) %>% unique() %>% arrange(age) %>% round(2)

# temperature
pred_temp <- as.data.frame(Effect(c("s.log.age", "s.c.temp"),
                                  model,
                                  xlevels = list(s.log.age = df_age$s.log.age,
                                                 s.c.temp = seq(0, 1, 0.1)))) %>%
  mutate(var = s.c.temp) %>%
  left_join(df_age) %>%
  mutate(fit_trans = exp(fit),
         upper_trans = exp(upper),
         lower_trans = exp(lower)) %>%
  select(age, var, fit_trans, upper_trans, lower_trans) %>%
  mutate(var_name = "Temperature",
         age_name = paste0("Age ", age))

# ssb
pred_ssb <- as.data.frame(Effect(c("s.log.age", "s.ssb.i"),
                                 model,
                                 xlevels = list(s.log.age = df_age$s.log.age,
                                                s.ssb.i = seq(0, 1, 0.1)))) %>%
  mutate(var = s.ssb.i) %>%
  left_join(df_age) %>%
  mutate(fit_trans = exp(fit),
         upper_trans = exp(upper),
         lower_trans = exp(lower)) %>%
  select(age, var, fit_trans, upper_trans, lower_trans) %>%
  mutate(var_name = "Spawning stock biomass",
         age_name = paste0("Age ", age))


# recruitment
pred_rec <- as.data.frame(Effect(c("s.log.age", "s.recruitment.i"),
                                 model,
                                 xlevels = list(s.log.age = df_age$s.log.age,
                                                s.recruitment.i = seq(0, 1, 0.1)))) %>%
  mutate(var = s.recruitment.i) %>%
  left_join(df_age) %>%
  mutate(fit_trans = exp(fit),
         upper_trans = exp(upper),
         lower_trans = exp(lower)) %>%
  select(age, var, fit_trans, upper_trans, lower_trans) %>%
  mutate(var_name = "Recruitment",
         age_name = paste0("Age ", age))

# TP
pred_tp <- as.data.frame(Effect(c("s.log.age", "s.TP"),
                                model,
                                xlevels = list(s.log.age = df_age$s.log.age,
                                               s.TP = seq(0, 1, 0.1)))) %>%
  mutate(var = s.TP) %>%
  left_join(df_age) %>%
  mutate(fit_trans = exp(fit),
         upper_trans = exp(upper),
         lower_trans = exp(lower)) %>%
  select(age, var, fit_trans, upper_trans, lower_trans) %>%
  mutate(var_name = "Total phosphorus",
         age_name = paste0("Age ", age))

# all extrinsic variables
pred <- bind_rows(pred_temp, pred_ssb, pred_rec, pred_tp)

#### summary
pred_sum <- pred %>% 
  filter(var %in% c(0,1)) %>%
  select(var_name, age, var, fit_trans) %>%
  group_by(var_name, age) %>%
  pivot_wider(names_from = var, values_from = fit_trans) %>%
  mutate(diff = (`1`-`0`)/`0`*100) %>%
  mutate(across(where(is.numeric), round, 2)) %>%
  arrange(age)

pred_sum <- pred_sum %>%
  mutate(var_name = factor(var_name, levels = c("Temperature", 
                                                "Spawning stock biomass",
                                                "Recruitment",
                                                "Total phosphorus")),
         sign = if_else(diff < 0, "Negative", "Positive"),
         source_name = "ORAS5")

#### plot
scale_shape = c("Temperature" = 15, 
                "Spawning stock biomass" = 2, 
                "Recruitment" = 5,
                "Total phosphorus" = 6)

ggplot(data = pred_sum) +
  geom_point(aes(x = age, y = abs(diff), shape = var_name, color = sign)) +
  scale_shape_manual(values = scale_shape,
                     labels = expression('Temperature'['population-anomaly'],
                                         'Spawning stock biomass', 
                                         'Recruitment',
                                         'Total phosphorus')) +
  theme(legend.text.align = 0) +
  labs(x = "Age",
       y = "Growth change (%)",
       shape = "Predictor",
       color = "Effect direction") +
  facet_grid(. ~ source_name)

p1 <- last_plot()

### extrinsic nemo-medusa ----
#### setup
# model
model <- m4_nu_nm_s

# data and model
data <- data_otl %>% 
  left_join(filter(data_temp, source == "nemo-medusa"),
            by = join_by(pop, year)) %>%
  filter(is.na(c.temp) == F) %>%
  left_join(data_sol, by = join_by(pop, year)) %>%
  left_join(data_nu, by = join_by(pop, year))
data <- data %>% filter(!is.na(TN) & !is.na(TP))
data <- data %>%
  mutate(s.log.age = s.(log.age),
         s.log.aac = s.(log.aac),
         s.c.temp = s.(c.temp),
         s.ssb.i = s.(ssb.i),
         s.recruitment.i = s.(recruitment.i),
         s.TP = s.(TP))

# df_age
df_age <- select(data, age, s.log.age) %>% unique() %>% arrange(age) %>% round(2)

# temperature
pred_temp <- as.data.frame(Effect(c("s.log.age", "s.c.temp"),
                                  model,
                                  xlevels = list(s.log.age = df_age$s.log.age,
                                                 s.c.temp = seq(0, 1, 0.1)))) %>%
  mutate(var = s.c.temp) %>%
  left_join(df_age) %>%
  mutate(fit_trans = exp(fit),
         upper_trans = exp(upper),
         lower_trans = exp(lower)) %>%
  select(age, var, fit_trans, upper_trans, lower_trans) %>%
  mutate(var_name = "Temperature",
         age_name = paste0("Age ", age))

# ssb
pred_ssb <- as.data.frame(Effect(c("s.log.age", "s.ssb.i"),
                                 model,
                                 xlevels = list(s.log.age = df_age$s.log.age,
                                                s.ssb.i = seq(0, 1, 0.1)))) %>%
  mutate(var = s.ssb.i) %>%
  left_join(df_age) %>%
  mutate(fit_trans = exp(fit),
         upper_trans = exp(upper),
         lower_trans = exp(lower)) %>%
  select(age, var, fit_trans, upper_trans, lower_trans) %>%
  mutate(var_name = "Spawning stock biomass",
         age_name = paste0("Age ", age))

# recruitment
pred_rec <- as.data.frame(Effect(c("s.log.age", "s.recruitment.i"),
                                 model,
                                 xlevels = list(s.log.age = df_age$s.log.age,
                                                s.recruitment.i = seq(0, 1, 0.1)))) %>%
  mutate(var = s.recruitment.i) %>%
  left_join(df_age) %>%
  mutate(fit_trans = exp(fit),
         upper_trans = exp(upper),
         lower_trans = exp(lower)) %>%
  select(age, var, fit_trans, upper_trans, lower_trans) %>%
  mutate(var_name = "Recruitment",
         age_name = paste0("Age ", age))


# all extrinsic variables
pred <- bind_rows(pred_temp, pred_ssb, pred_rec)

#### summary
pred_sum <- pred %>% 
  filter(var %in% c(0,1)) %>%
  select(var_name, age, var, fit_trans) %>%
  group_by(var_name, age) %>%
  pivot_wider(names_from = var, values_from = fit_trans) %>%
  mutate(diff = (`1`-`0`)/`0`*100) %>%
  mutate(across(where(is.numeric), round, 2)) %>%
  arrange(age)

# summary all ages
pred_sum <- pred_sum %>%
  mutate(var_name = factor(var_name, levels = c("Temperature", 
                                                "Spawning stock biomass",
                                                "Recruitment",
                                                "Total phosphorus")),
         sign = if_else(diff < 0, "Negative", "Positive"),
         source_name = "NEMO-MEDUSA")

#### plot
scale_shape = c("Temperature" = 15, 
                "Spawning stock biomass" = 2, 
                "Recruitment" = 5,
                "Total phosphorus" = 6)

ggplot(data = pred_sum) +
  geom_point(aes(x = age, y = abs(diff), shape = var_name, color = sign)) +
  scale_shape_manual(values = scale_shape,
                     labels = expression('Temperature'['population-anomaly'],
                                         'Spawning stock biomass', 
                                         'Recruitment',
                                         'Total phosphorus')) +
  theme(legend.text.align = 0) +
  theme(legend.position = "none") +
  labs(x = "Age",
       y = "Growth change (%)",
       shape = NULL,
       color = NULL) +
  facet_grid(. ~ source_name)

p2 <- last_plot()

### merge and save plot ----
p1 + p2 +
  plot_annotation(tag_levels = 'A') + 
  plot_layout(guides = "collect")

ggsave(last_plot(), file = file.path(dir_report, "figS19_extrinsic comparison with nutrient.tiff"),
       width = 17, height = 5.75,
       units = "cm",  scaling = 0.8,
       dpi = 600)
ggsave(last_plot(), file = file.path(dir_report, "figS19_extrinsic comparison with nutrient.eps"),
       width = 17, height = 5.75,
       units = "cm", 
       dpi = 600, 
       device = "eps")

## fig S20 - density effect with and without standardisation ----
#### save file
df_all <- read_rds(file.path(dir_output, "robustness_density without standardisation.rds"))

#### plot
ggplot(data = df_all, 
       aes(x = estimate, 
           y = data_name)) +
  geom_point() +
  geom_linerange(aes(xmin = estimate - std.error,
                     xmax = estimate + std.error)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_grid(source_name ~ term_name) +
  labs(x = "Parameter estimate",
       y = NULL) +
  scale_x_continuous(breaks = seq(0,2.5,0.5))

#### save plot
ggsave(last_plot(), file = file.path(dir_report, "figS20_density effect vs standardisation.tiff"),
       width = 17, height = 11.5,
       units = "cm",  scaling = 0.8,
       dpi = 600)
ggsave(last_plot(), file = file.path(dir_report, "figS20_density effect vs standardisation.eps"),
       width = 17, height = 11.5,
       units = "cm", 
       dpi = 600, 
       device = "eps")

## fig S21 - variance individual plasticity ----
#### setup
cor_isimip <- read_rds(file.path(dir_output, paste0("ind.plastic_cor.test_", "isimip", ".rds"))) %>%
  mutate(source = "isimip")
cor_oras5 <- read_rds(file.path(dir_output, paste0("ind.plastic_cor.test_", "oras5", ".rds"))) %>%
  mutate(source = "oras5")
cor_nm <- read_rds(file.path(dir_output, paste0("ind.plastic_cor.test_", "nemo-medusa", ".rds"))) %>%
  mutate(source = "nemo-medusa")

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

cor_all <- bind_rows(cor_isimip, cor_oras5, cor_nm)
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
       shape = "Increment measurement range") 

#### save plot
ggsave(last_plot(), file = file.path(dir_report, "figS21_variance individual plasticity.tiff"),
       width = 17, height = 11.5,
       units = "cm",  scaling = 0.8,
       dpi = 600)
ggsave(last_plot(), file = file.path(dir_report, "figS21_variance individual plasticity.eps"),
       width = 17, height = 11.5,
       units = "cm", 
       dpi = 600, 
       device = "eps")

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

## table S4 - model - random structure ----
#### setup
# refit all models
m1a <- lmer(log.increment ~ 1 + log.age*method + log.age*datasource + log.aac*method + log.aac*datasource + 
              (1 | fishid), data = data, REML = T)
m2a <- lmer(log.increment ~ 1 + log.age*method + log.age*datasource + log.aac*method + log.aac*datasource + 
              (1 | fishid) + (1 | pop) + (1 | pop.year), data = data, REML = T,
            control = lmerControl(optimizer ="Nelder_Mead"))
m2b <- lmer(log.increment ~ 1 + log.age*method + log.age*datasource + log.aac*method + log.aac*datasource + 
             (1 | fishid) + (1 | pop) + (1 + log.age | pop.year), data = data, REML = T)
m2e <- lmer(log.increment ~ 1 + log.age*method + log.age*datasource + log.aac*method + log.aac*datasource + 
              (1 | fishid) + (1 | pop) + (1 | pop.cohort), data = data, REML = T)
m2k <- lmer(log.increment ~ 1 + log.age*method + log.age*datasource + log.aac*method + log.aac*datasource + 
              (1 | fishid) + (1 | pop) + (1 | pop.year) + (1 + log.age | pop.cohort), data = data, REML = T)

# Model comparison
models <- list(m1a, m2a, m2b, m2e, m2k)
Modnames <- c('m1a', 'm2a', 'm2b', 'm2e', 'm2k')
m_compare <- as.data.frame(aictab(cand.set = models, modnames = Modnames, sort = TRUE))

#### table
df <- tibble(Modnames = c("m1a", "m2a", "m2b", "m2e", "m2k"),
             FishID = c(rep("+", 5)),
             Population = c("", rep("+", 4)),
             `Population:Year` = c("", "+", "+", "", "+"),
             `Age|Population:Year` = c("", "", "+", "", ""),
             `Population:Cohort` = c("", "", "", "+", "+"),
             `Age|Population:Cohort` = c(rep("", 4), "+")) %>%
  left_join(m_compare) %>%
  select(-K, -(ModelLik:Cum.Wt)) %>%
  select(- Modnames) %>%
  arrange(AICc) %>%
  mutate(AICc = round(AICc, 2),
         Delta_AICc = round(Delta_AICc, 2)) 

#### save table
col.header <- c("FishID", 
                "Population", 
                "Population:Year", 
                "Age|Population:Year", 
                "Population:Cohort", 
                "Age|Population:Cohort", 
                "AICc",
                "\u0394AIC")

tab_df(df,
       CSS = css_list,
       col.header = col.header,
       file = file.path(dir_report, "tableS4_random structure.html"))

## table S5 - model - intrinsic structure ----
#### setup
m_compare <- read_rds(file.path(dir_output, "intrinsic.model_comparison.rds"))

df_aic <- as.data.frame(m_compare) %>% 
  select(AICc, delta) %>%
  round(2)

df_var <- as.data.frame(m_compare) %>% 
  select(1:7) %>% 
  mutate_all(fun_plus) 

#### table
df <- bind_cols(df_var, df_aic) %>%
  slice(1:5) %>%
  select(`(Intercept)`, log.age, log.aac, method, datasource, `log.age:method`, `datasource:log.age`, AICc, delta)

#### save table
names(df)
col.header <- c("Intercept", 
                 "Age", 
                 "Age at capture", 
                 "Preparation method", 
                 "Reading institute", 
                 "Age * Preparation method", 
                 "Age * Reading institute", 
                 "AIC", 
                 "\u0394AIC")
tab_df(df,
       CSS = css_list,
       col.header = col.header,
       file = file.path(dir_report, "tableS5_intrinsic structure.html")) 

## table S6 - model - extrinsic structure ----
#### setup
list_model <- tibble(source = c("isimip", "oras5", "nemo-medusa"),
                     id_best = c(4, 13, 9)) # best model id see 2_analyze data.R
df <- tibble()
for(i in 1:3) {
  m_compare <- read_rds(file.path(dir_output, paste0("extrinsic.model_comparison_", list_model$source[i], ".rds")))
  
  df_aic <- as.data.frame(m_compare) %>% 
    select(AICc, delta) %>%
    round(2)
  
  df_info <- as.data.frame(m_compare) %>% 
    select(source, year_range, var_added, id)
  
  df_var <- as.data.frame(m_compare) %>% 
    select(-(df:id)) %>% 
    mutate_all(fun_plus) 
  
  if(list_model$id_best[i] <= 5) {
    df_temp <- bind_cols(df_var, df_aic, df_info) %>%
      slice(1:5)
  } else {
    df_temp <- bind_cols(df_var, df_aic, df_info) %>%
      slice(1:4, list_model$id_best[i])
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
         AICc, delta, id, var_added, source_name, year_range)

#### save table
names(df)
col.header <- c("Spawning stock biomass", 
                "Spawning stock biomass * Age", 
                "Recruitment", 
                "Recruitment * Age", 
                "Fishing mortality", 
                "Temperature_population-average",
                "Temperature_population-anomaly",
                "Temperature_population-anomaly * Age",
                "Temperature_population-anomaly * Spawning stock biomass",
                "Temperature_population-anomaly * Recruitment",
                "Temperature_population-anomaly * Fishing mortality",
                "Temperature_population-anomaly * Temperature_population-average",
                "AICc",
                "\u0394AIC",
                "ModelID",
                "No. variable added",
                "Temperature dataset",
                "Year range")

tab_df(df,
       CSS = css_list,
       col.header = col.header,
       file = file.path(dir_report, "tableS6_extrinsic structure.html")) 

## table S7 - model - extrinsic extended - nutrient ----
list_model <- tibble(source = c("isimip", "oras5", "nemo-medusa"))
df <- tibble()
for(i in 1:3) {
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
         c.temp, `c.temp:log.age`, 
         TN, TP,
         AICc, delta, source_name) %>%
  replace_na(list(c.temp = "", `c.temp:log.age` = ""))

#### save table
names(df)
col.header <- c("Spawning stock biomass", 
                "Recruitment", 
                "Temperature_population-anomaly",
                "Temperature_population-anomaly * Age",
                "Total nitrogen",
                "Total phosphorus",
                "AICc",
                "\u0394AIC",
                "Temperature dataset")

tab_df(df,
       CSS = css_list,
       col.header = col.header,
       file = file.path(dir_report, "tableS7_extrinsic structure_nutrient.html")) 

## table S8 - model - extrinsic extended structure ----
#### setup
list_model <- tibble(source = c("isimip", "oras5", "nemo-medusa"),
                     id_best = c(5, 5, 1)) # best model id see 2_analyze data.R
df <- tibble()
for(i in 1:3) {
  m_compare <- read_rds(file.path(dir_output, paste0("extrinsic.model.ext_comparison_", list_model$source[i], ".rds")))
  
  df_aic <- as.data.frame(m_compare) %>% 
    select(AICc, delta) %>%
    round(2)
  
  df_info <- as.data.frame(m_compare) %>% 
    select(source, year_range, var_added, id)
  
  df_var <- as.data.frame(m_compare) %>% 
    select(-(df:id)) %>% 
    mutate_all(fun_plus) 
  
  if(list_model$id_best[i] <= 5) {
    df_temp <- bind_cols(df_var, df_aic, df_info) %>%
      slice(1:5)
  } else {
    df_temp <- bind_cols(df_var, df_aic, df_info) %>%
      slice(1:4, list_model$id_best[i])
  }
  
  df <- bind_rows(df, df_temp) 
}

df <- df %>%
  left_join(df_temp_name) %>%
  select(ssb.i, 
         recruitment.i, 
         f, 
         c.temp_within, `c.temp_within:log.age`, `c.temp_within:ssb.i`, `c.temp_within:recruitment.i`, `c.temp_within:f`, `c.temp_between:c.temp_within`,
         c.temp_between, `c.temp_between:log.age`,
         AICc, delta, id, var_added, source_name, year_range)

#### save table
names(df)
col.header <- c("Spawning stock biomass", 
                "Recruitment", 
                "Fishing mortality",
                "Temperature_individual-average",
                "Temperature_individual-average * Age",
                "Temperature_individual-anomaly",
                "Temperature_individual-anomaly * Age",
                "Temperature_individual-anomaly * Spawning stock biomass",
                "Temperature_individual-anomaly * Recruitment",
                "Temperature_individual-anomaly * Fishing mortality",
                "Temperature_individual-anomaly * Temperature_individual-average",
                "AICc",
                "\u0394AIC",
                "ModelID",
                "No. variable added",
                "Temperature dataset",
                "Year range")

tab_df(df,
       CSS = css_list,
       col.header = col.header,
       file = file.path(dir_report, "tableS8_extrinsic extended structure.html")) 

## table S9 - model - extrinsic extended structure - nutrient ----
list_model <- tibble(source = c("isimip", "oras5", "nemo-medusa"))
df <- tibble()
for(i in 1:3) {
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
         c.temp_within, `c.temp_within:log.age`, `c.temp_within:ssb.i`, `c.temp_within:recruitment.i`, `c.temp_within:f`, 
         c.temp_between, `c.temp_between:log.age`,
         TN, TP,
         AICc, delta, source_name) %>%
  replace_na(list(f = "",
                  `c.temp_between:log.age` = "", 
                  `c.temp_within:ssb.i` = "",
                  `c.temp_within:recruitment.i` = "",
                  `c.temp_within:f` = ""
                  ))

#### save table
names(df)
col.header <- c("Spawning stock biomass", 
                "Recruitment", 
                "Fishing mortality",
                "Temperature_individual-average",
                "Temperature_individual-average * Age",
                "Temperature_individual-anomaly",
                "Temperature_individual-anomaly * Age",
                "Temperature_individual-anomaly * Spawning stock biomass",
                "Temperature_individual-anomaly * Recruitment",
                "Temperature_individual-anomaly * Fishing mortality",
                "Total nitrogen",
                "Total phosphorus",
                "AICc",
                "\u0394AIC",
                "Temperature dataset")

tab_df(df,
       CSS = css_list,
       col.header = col.header,
       file = file.path(dir_report, "tableS9_extrinsic structure_nutrient.html")) 

## table S10 - model summary  ----
pred_labels <- c(
  "Intercept",
  "Age",
  "Age at capture",
  "Reading Institute (WUR)",
  "Reading Institute (WUR) * Age",
  "Spawning Stock Biomass",
  "Recruitment",
  "Temperature_population-anomaly",
  "Temperature_population-anomaly * Age"
)

# summary table
tab_model(m3, m4_isimip, m4_oras5, m4_nm, show.se = TRUE, 
                 show.ci = NULL, show.p = FALSE,
                 dv.labels = c("Intrinsic model", 
                               "Population-level extrinsic model (ISIMIP)", 
                               "Population-level extrinsic model (ORAS5)", 
                               "Population-level extrinsic model (NEMO-MEDUSA)"),
                 pred.labels = pred_labels,
                 string.pred = "Fixed Effects",
                 string.est = "Estimate (SE)",
                 collapse.se = TRUE,
                 show.icc = FALSE,
                 digits.rsq = 2,
                 CSS = css_list,
                 file = file.path(dir_report, "tableS10_model_summary.html")
)


## table S11 - model summary - nutrient ----
# labels
pred_labels <- c(
  "Spawning Stock Biomass",
  "Recruitment",
  "Total Phosphorus",
  "Temperature",
  "Temperature * Age"
)

# summary table
tab_model(m4_nu_isimip, m4_nu_oras5, m4_nu_nm, show.se = TRUE, 
          terms = c("ssb.i", "recruitment.i", "c.temp", "log.age:c.temp", "TP"),
          order.terms = c(1,2,3,4,5),
          show.ci = NULL, show.p = FALSE,
          show.re.var = FALSE,
          dv.labels = c("Population-level extrinsic model (ISIMIP)", 
                        "Population-level extrinsic model (ORAS5)", 
                        "Population-level extrinsic model (NEMO-MEDUSA)"),
          pred.labels = pred_labels,
          string.pred = "Fixed Effects",
          string.est = "Estimate (SE)",
          collapse.se = TRUE,
          show.icc = FALSE,
          digits = 3,
          digits.rsq = 2,
          CSS = css_list,
          file = file.path(dir_report, "tableS11_model_summary_nu.html"))

## table S12 - model summary - extrinsic extended - nutrient  ----
# labels
pred_labels <- c(
  "Spawning Stock Biomass",
  "Recruitment",
  "Temperature_individual-anomaly",
  "Temperature_individual-average",
  "Total Phosphorus",
  "Temperature_individual-anomaly * Age",
  "Temperature_individual-average * Age",
  "Temperature_individual-anomaly * Spawning Stock Biomass",
  "Fishing mortality",
  "Temperature_individual-anomaly * Fishing mortality",
  "Temperature_individual-anomaly * Recruitment"
)

# summary table
tab_model(m5_nu_isimip, m5_nu_oras5, m5_nu_nm, show.se = TRUE, 
                 terms = c("ssb.i", "recruitment.i", "f",
                           "c.temp_within", "log.age:c.temp_within",
                           "c.temp_between", "log.age:c.temp_between",
                           "ssb.i:c.temp_within", "recruitment.i:c.temp_within", "c.temp_within:f", 
                           "TP"),
                 order.terms = c(1,2,9,5,3,6,8,11,10,4,7),
                 pred.labels = pred_labels,
                 show.ci = NULL, show.p = FALSE,
                 show.re.var = FALSE,
                 dv.labels = c("Individual-level extrinsic model (ISIMIP)", 
                               "Individual-level extrinsic model (ORAS5)", 
                               "Individual-level extrinsic model (NEMO-MEDUSA)"),
                 string.pred = "Fixed Effects",
                 string.est = "Estimate (SE)",
                 collapse.se = TRUE,
                 show.icc = FALSE,
                 digits = 3,
                 digits.rsq = 2,
                 CSS = css_list,
                 file = file.path(dir_report, "tableS12_model_ext_summary_nu.html")
)

## table S13 - c.temp variance ratio ----
list_pair <- tibble(pop_pair = c("4bc/8ab", "7a/8ab", "4bc/7a"),
                    pop_pair_name = factor(c("North Sea/Bay of Biscay", 
                                             "Irish Sea/Bay of Biscay", 
                                             "North Sea/Irish Sea"),
                                           levels = c("North Sea/Bay of Biscay", 
                                                      "Irish Sea/Bay of Biscay", 
                                                      "North Sea/Irish Sea")))
pop_pair <- tibble(pop1 = c("4bc", "4bc", "7a"),
                   pop2 = c("7a", "8ab", "8ab"))
source_list = c("isimip", "oras5", "nemo-medusa")

df <- tibble()
for(i in 1:3) {
  for(p in 1:3) {
    print(paste0("processing ", pop_pair$pop1[p], "/", pop_pair$pop2[p]))
    
    pop1 <- pop_pair$pop1[p]
    pop2 <- pop_pair$pop2[p]
    
    df_pop1 <- data_temp %>% filter(pop == pop1, source == source_list[i])
    df_pop2 <- data_temp %>% filter(pop == pop2, source == source_list[i])
    
    var_test <- var.test(df_pop1$c.temp, df_pop2$c.temp)
    
    df_temp <- tibble(source = source_list[i],
                      pop_pair = paste0(pop1, "/", pop2),
                      var_ratio = sprintf("%.*f", 2, var_test$statistic),
                      ci95 = paste0(sprintf("%.*f", 2, var_test$conf.int[1]),
                                    " - ",
                                    sprintf("%.*f", 2, var_test$conf.int[2])),
                      p_value = sprintf("%.*f", 2, var_test$p.value))
    
    df <- bind_rows(df, df_temp)
  }
}

df <- df %>% 
  left_join(df_temp_name) %>%
  left_join(list_pair) %>%
  select(source_name, pop_pair_name, var_ratio, ci95, p_value) %>% 
  arrange(source_name, pop_pair_name)

#### save table
col.header <- c("Temperature dataset", 
                "Population pair", 
                "Variance ratio",
                "95% confidence interval",
                "p-value")

tab_df(df,
       CSS = css_list,
       col.header = col.header,
       file = file.path(dir_report, "tableS13_c.temp variance.html"))
