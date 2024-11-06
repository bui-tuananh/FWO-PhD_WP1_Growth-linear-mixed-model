# BEFORE: 
#       2_analyze data
# AFTER:
#       WP1/report

# 1. SETUP --------------
library(tidyverse)
library(sf)
library(stars)
library(RColorBrewer)
library(patchwork)

# FUNCTION
#centring function
c. <- function (x) {(x - mean(x))} 

# WORKING DIRECTORY
dir_report <- "./report" # not indicate dir_report 

# POPULATION NAME
# create df_pop
df_pop <- tibble(IcesAreaGroup = c("4bc", "7a", "8ab"),
                 pop = c("4bc", "7a", "8ab"),
                 pop.name = factor(c("North Sea", "Irish Sea", "Bay of Biscay"),
                                   levels = c("North Sea", "Irish Sea", "Bay of Biscay")))

# THEME
theme_set(theme_classic()) 
theme_update(panel.border = element_rect(colour = "black", fill=NA))

# 2. LOAD DATA ------------
## 2.1. Otolith data ---------------------------------------------------

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

## 2.2. Temperature data ----
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

## 2.3. Fishing data ----
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

## 2.4. Nutrient data ----
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

## 2.2. GIS data -----------

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

# distribution area - datras
dir_ices <- "./data/ices"
ices_datras <- read_sf(file.path(dir_ices, "hl_loc_4abc7a8ab.gpkg"))

# sbt - oras5 (1958-2019) 
dir_temp <- "./data/temp"
oras5 <- read_stars(file.path(dir_temp, 
                              "oras5.tif")) 
oras5_df <- oras5 %>% 
  st_crop(ices_datras) %>%
  slice(band, 1:744) %>%
  st_apply(1:2, mean) %>% 
  as.data.frame() %>%
  filter(is.na(mean) == F)

#note: 2019 only - 733-744

# 3. PLOT DATA -------------

## main study area --------

# setup
# use quantile scale to aid visualization
## setup breaks and labels 
legend <- tibble(breaks = seq(0, 1, 0.1),
                 labels = quantile(oras5_df$mean, probs = seq(0, 1, 0.1))) %>%
  mutate(labels = round(labels, 1))
## create color scale
fill_colours = rev(brewer.pal(nrow(legend), "RdYlBu"))
## add quantile value to dataframe to plot
oras5_df$quant <- ecdf(oras5_df$mean)(oras5_df$mean)  

# pop name
sf_pop <- data.frame(name = c("North Sea", "Irish Sea", "Bay of Biscay"),
                     x = c(3, -5, -6.6),
                     y = c(57.75, 55.25, 48.25)) %>% 
  st_as_sf(coords = c("x", "y"), crs = 4326) 

# plot
st_bbox(ices_area)
ggplot() +
  geom_tile(data = oras5_df, aes(x = x, y = y, fill = quant)) +
  geom_sf(data = ices_datras, fill = NA, color = "black", linewidth = 0.5) +
  geom_sf(data = countries, fill = "white", color = "grey",linewidth = 0.25) + 
  geom_sf(data = continents, fill = NA, color = "grey", linewidth = 0.5) + 
  geom_sf(data = ices_area, fill = NA, color = "black", linewidth = 0.5, linetype = "dashed") +
  geom_sf_text(data = sf_pop, aes(label = name), size = 6/.pt) +
  coord_sf(xlim = c(-9, 9.5), ylim = c(43, 62.5), expand = FALSE) +
  theme_bw() +
  scale_fill_gradientn(
    colours = fill_colours,
    breaks = slice(legend, c(1,3,5,7,9,11))$breaks,
    labels = sprintf("%.1f", slice(legend, c(1,3,5,7,9,11))$labels),
    limits = c(0,1),
    guide = guide_colorbar(title = "Temperature (°C)\n(1958-2019)", 
                           title.position = "top", #left
                           #title.vjust = 0.7,
                           #barwidth = 6,
                           title.hjust = 0.5)) +
  labs(x = "Longitude",
       y = "Latitude") +
  theme(legend.position = c(0.99, 0.001),
        legend.justification = c(1, 0),
        legend.direction="horizontal",
        legend.key.height = unit(10, "pt"),
        legend.key.width = unit(15, "pt"),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 5)
  ) 

p_main <- last_plot()

## temp ----
data_temp_oras5 <- data_temp %>% 
  filter(source == "oras5") %>%
  left_join(df_pop)

## temp
ggplot(data = data_temp_oras5, aes(x = year, y = temp, color = pop.name)) + 
  geom_line() +
  scale_color_brewer(palette = "Dark2", direction = -1) +
  labs(y = "Temperature\n(°C)") +
  theme(legend.position = c(0, 0.5),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill='transparent'),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.key.height = unit(7, "pt"),
        legend.text = element_text(size = 6)) +
  theme(legend.position = "none") +
  xlim(1958, 2019)
p_temp <- last_plot()

## ssb, recruitment, f ----
data_sol <- data_sol %>% left_join(df_pop)

## ssb.i
ggplot(data = data_sol, aes(x = year, y = ssb.i, color = pop.name)) + 
  geom_line() +
  scale_color_brewer(palette = "Dark2", direction = -1) +
  labs(y = "Spawning stock biomass\n(tonne/km²)") +
  theme(legend.position = c(0, 0.5),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill='transparent'),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.key.height = unit(7, "pt"),
        legend.text = element_text(size = 6)) +
  theme(legend.position = "none") +
  xlim(1958, 2019)
p_ssb <- last_plot()

## recruitment.i
ggplot(data = data_sol, aes(x = year, y = recruitment.i, color = pop.name)) + 
  geom_line() +
  scale_color_brewer(palette = "Dark2", direction = -1) +
  theme(legend.position = "none",
        axis.title.x = element_blank()) +
  labs(y = "Recruitment\n(thoudsand/km²)") +
  xlim(1958, 2019)
p_rec <- last_plot()

## f
ggplot(data = data_sol, aes(x = year, y = f, color = pop.name)) + 
  geom_line() +
  scale_color_brewer(palette = "Dark2", direction = -1) +
  theme(legend.position = "none",
        axis.title.x = element_blank()) +
  labs(y = "Fishing mortality") +
  xlim(1958, 2019)
p_f <- last_plot()

## nutrient ----
data_nu2 <- data_nu %>% 
  pivot_longer(cols = TN:TP) %>%
  group_by(name) %>%
  mutate(s_value = scale(value)) %>%
  left_join(df_pop)

# add Irish Sea, Bay of Biscay for legend
data_nu_add <- data_nu2 %>%
  slice(2) %>%
  select(-pop.name) %>%
  bind_cols(tibble(pop.name = c("Irish Sea", "Bay of Biscay")))

data_nu2 <- bind_rows(data_nu2, data_nu_add) %>%
   mutate(pop.name = factor(pop.name, levels = c("North Sea", "Irish Sea", "Bay of Biscay")))

ggplot(data = data_nu2, aes(x = year, y = s_value, 
                            linetype = name,
                            color = pop.name)) + 
  geom_line() +
  theme(legend.position = c(0, 1.05),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill='transparent'),
        #axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.key.height = unit(7, "pt"),
        legend.text = element_text(size = 6)) +
  labs(x = "Year",
       y = "Scaled value",
       linetype = "Nutrient",
       color = "Nutrient")  +
  scale_color_brewer(palette = "Dark2", direction = -1) +
  scale_linetype_manual(values = c("TN" = "dashed", "TP" = "solid"),
                        labels = c("Total nitrogen", 
                                   "Total phosphorus")) +
  #scale_color_manual(values=c("#7570b3", "#7570b3")) +
  #guides(color = "none") +
  xlim(1958, 2019)

p_nu <- last_plot()

## arrange plot ----
p_sampling <- p_main + 
  (p_temp / p_ssb / p_rec / p_f/ p_nu) + 
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 9),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 7),
        plot.margin = unit(c(1,1,1,1), "pt"))

#p_sampling

## save plot -----
ggsave(p_sampling, file = file.path(dir_report, "fig1_sampling site.pdf"),
       device = cairo_pdf,
       width = 19, height = 19,
       units = "cm")
ggsave(p_sampling, file = file.path(dir_report, "fig1_sampling site.png"),
       width = 19, height = 15, 
       units = "cm",  
       dpi = 1000)

