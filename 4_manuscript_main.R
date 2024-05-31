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

# DEFAULT THEME
theme_set(theme_classic()) 
theme_update(panel.border = element_rect(colour = "black", fill=NA))

# FUNCTION
#centring function
c. <- function (x) {(x - mean(x))} 
# scale function
s. <- function (x) {(x - mean(x))/sd(x)} 

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

## figure 1 - sampling site ----
source("./3_report_sampling site.R")

## figure 2 - otolith measure (ppt) ----
## figure 3 - conceptual diagram temperature components ----
source("./3_report_conceptual temperature effect.R")

## figure 4 - intrinsic effects: age, age at capture ----

list_plot <- list()
data <- m3@frame

### age ----
# setup
pred <- as.data.frame(Effect(c("log.age"), 
                             m3, 
                             xlevels = list(log.age = unique(data$log.age)))) %>%
  mutate(age = round(exp(log.age),0)) %>%
  arrange(age)

# plot
ggplot(data = pred) + 
  geom_line(aes(x = age, 
                y = exp(fit))) +
  geom_ribbon(aes(x = age,
                  ymin = exp(lower),
                  ymax = exp(upper)), 
              alpha = 0.5) +
  labs(x = "Age (years)",
       y = "Predicted increment growth (μm)")

list_plot[[1]] <- last_plot()

### aac*pop ----
pred <- as.data.frame(Effect(c("log.aac"), 
                             m3, 
                             xlevels = list(log.aac = unique(data$log.aac)))) %>%
  mutate(aac = round(exp(log.aac),0)) %>%
  arrange(aac)

#### plot
ggplot(data = pred) + 
  geom_line(aes(x = aac, 
                y = exp(fit)),
            linewidth = 1) +
  geom_ribbon(aes(x = aac,
                  ymin = exp(lower),
                  ymax = exp(upper)), 
              alpha = 0.5) +
  labs(x = "Age at capture (years)",
       y = "Predicted increment growth (μm)",
       color = "Population",
       fill = "Population")

list_plot[[2]] <- last_plot()

### merge and save plot ----
# merge plot
(list_plot[[1]] | list_plot[[2]]) + 
  plot_annotation(tag_levels = 'A') + 
  plot_layout(guides = "collect", axis_titles = "collect") &
  theme(plot.tag.position  = c(0.95, 0.9),
        axis.title = element_text(size = 9))

# save file
ggsave(last_plot(), file = file.path(dir_report, "fig4_intrinsic effects.pdf"),
       device = cairo_pdf,
       width =  19, height = 6.3,
       units = "cm"
       ) 
ggsave(last_plot(), file = file.path(dir_report, "fig4_intrinsic effects.tiff"),
       width =  19, height = 6.3,
       units = "cm",
       dpi = 1000)

## figure 5 - growth temporal trend (year random effect) ----
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

age1 <- pred_year %>% mutate(ranef = intercept + log.age*df_age$log.age[1],
                             ranef_upper = intercept + intercept_se*1.96 + (log.age + log.age_se*1.96)*df_age$log.age[1],
                             ranef_lower = intercept - intercept_se*1.96 + (log.age - log.age_se*1.96)*df_age$log.age[1],
                             age = "Age 1")
age5 <- pred_year %>% mutate(ranef = intercept + log.age*df_age$log.age[5],
                             ranef_upper = intercept + intercept_se*1.96 + (log.age + log.age_se*1.96)*df_age$log.age[5],
                             ranef_lower = intercept - intercept_se*1.96 + (log.age - log.age_se*1.96)*df_age$log.age[5],
                             age = "Age 5")
age_all <- bind_rows(age1, age5)

#### plot
ggplot() + 
  geom_line(data = age_all, 
            aes(x = year, y = ranef)) +
  geom_ribbon(data = age_all,
              aes(x = year, 
                  ymin = ranef_lower,
                  ymax = ranef_upper),
              alpha = 0.5) +
  geom_hline(yintercept = 0, 
             linetype="dashed") +
  facet_grid(pop.name ~ age)  +
  labs(x = "Year",
       y = "Year random effect") +
  theme(axis.title = element_text(size = 9))

#### save plot
ggsave(last_plot(), file = file.path(dir_report, "fig5_growth temporal trend.pdf"),
       device = cairo_pdf,
       width =  19, height = 12.6,
       units = "cm") 
ggsave(last_plot(), file = file.path(dir_report, "fig5_growth temporal trend.tiff"),
       width =  19, height = 12.6,
       units = "cm",
       dpi = 1000)

## figure 6 - c.temp, c.temp_within, c.temp_between ----

list_plot <- list()

### c.temp ----
#### setup
# note: fit within range [-1,1] as temperatures has different range 

list_model <- list(#"isimip" = m4_isimip,
  "oras5" = m4_oras5,
  "nemo-medusa" = m4_nm
)

pred <- tibble()
for (i in 1:2) {
  model <- list_model[[i]]
  data <- model@frame
  temp_source <- names(list_model[i])
  
  pred_temp <- as.data.frame(Effect(c("log.age", "c.temp"), 
                                    model, 
                                    xlevels = list(log.age = unique(data$log.age),
                                                   c.temp = seq(-1, 1, 0.1)))) %>%
    mutate(age = round(exp(log.age),1)) %>%
    mutate(source = temp_source)
  
  pred <- bind_rows(pred, pred_temp)
}

df_age <- tibble(age = c(1,5),
                 age_name = c("Age 1", "Age 5"))

pred <- pred %>% 
  left_join(df_age) %>%
  left_join(df_temp_name)

#### plot
# add a point of ISIMIP to pred to merge legend 
pred_add <- pred %>% 
  filter(age %in% c(1,5), 
         c.temp == 0,
         source_name == "ORAS5") %>%
  mutate(source_name = "ISIMIP")

pred <- bind_rows(pred, pred_add) %>%
  mutate(source_name = factor(source_name, levels = c("ISIMIP", "ORAS5", "NEMO-MEDUSA")))

col_scale <- c("ISIMIP" = "#F8766D", "ORAS5" = "#00BFC4", "NEMO-MEDUSA" = "#7CAE00")

p1 <- ggplot(data = pred %>% filter(age %in% c(1))) + 
  geom_line(aes(x = c.temp, 
                y = exp(fit),
                color = source_name),
            linewidth = 1) +
  geom_ribbon(aes(x = c.temp,
                  ymin = exp(lower),
                  ymax = exp(upper),
                  fill = source_name), 
              alpha = 0.1) +
  facet_grid(. ~ age_name)  +
  labs(x = expression('T'['population-anomaly'] * ' (°C)'),
       y = "Predicted increment growth (μm)",
       color = "Dataset",
       fill = "Dataset") +
  scale_color_manual(values = col_scale) +
  scale_fill_manual(values = col_scale)

p2 <- ggplot(data = pred %>% filter(age %in% c(5))) + 
  geom_line(aes(x = c.temp, 
                y = exp(fit),
                color = source_name),
            linewidth = 1) +
  geom_ribbon(aes(x = c.temp,
                  ymin = exp(lower),
                  ymax = exp(upper),
                  fill = source_name), 
              alpha = 0.1) +
  facet_grid(. ~ age_name) +
  labs(x = expression('T'['population-anomaly'] * ' (°C)'),
       y = "Predicted increment growth (μm)",
       color = "Dataset",
       fill = "Dataset") +
  scale_color_manual(values = col_scale) +
  scale_fill_manual(values = col_scale)

list_plot[[1]] <- p1 
list_plot[[2]] <- p2 


### c.temp_between ----
#### setup
# note: fit within range [-1,1] as temperatures has different range 
# note: oras5 has smaller range [-0.8,0.9]

list_model <- list("isimip" = m5_isimip,
                   "oras5" = m5_oras5,
                   "nemo-medusa" = m5_nm
)

pred <- tibble()
for (i in 1:3) {
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

df_age <- tibble(age = c(1,5),
                 age_name = c("Age 1", "Age 5"))

pred <- pred %>% 
  left_join(df_age) %>%
  left_join(df_temp_name)

#### plot
col_scale <- c("ISIMIP" = "#F8766D", "ORAS5" = "#00BFC4", "NEMO-MEDUSA" = "#7CAE00")

p1 <- ggplot(data = pred %>% filter(age %in% c(1))) + 
  geom_line(aes(x = c.temp_between, 
                y = exp(fit),
                color = source_name),
            linewidth = 1) +
  geom_ribbon(aes(x = c.temp_between,
                  ymin = exp(lower),
                  ymax = exp(upper),
                  fill = source_name), 
              alpha = 0.1) +
  facet_grid(. ~ age_name)  +
  labs(x = expression('T'['individual-average'] * ' (°C)'),
       y = "Predicted increment growth (μm)",
       color = "Dataset",
       fill = "Dataset") +
  scale_color_manual(values = col_scale) +
  scale_fill_manual(values = col_scale)

p2 <- ggplot(data = pred %>% filter(age %in% c(5))) + 
  geom_line(aes(x = c.temp_between, 
                y = exp(fit),
                color = source_name),
            linewidth = 1) +
  geom_ribbon(aes(x = c.temp_between,
                  ymin = exp(lower),
                  ymax = exp(upper),
                  fill = source_name), 
              alpha = 0.1) +
  facet_grid(. ~ age_name) +
  labs(x = expression('T'['individual-average'] * ' (°C)'),
       y = "Predicted increment growth (μm)",
       color = "Dataset",
       fill = "Dataset") +
  scale_color_manual(values = col_scale) +
  scale_fill_manual(values = col_scale)

list_plot[[3]] <- p1 
list_plot[[4]] <- p2 

### c.temp_within ----
#### setup
# note: fit within range [-1,1] as temperatures has different range 

list_model <- list("isimip" = m5_isimip,
                   "oras5" = m5_oras5,
                   "nemo-medusa" = m5_nm
)

pred <- tibble()
for (i in 1:3) {
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

df_age <- tibble(age = c(1,5),
                 age_name = c("Age 1", "Age 5"))

pred <- pred %>% 
  left_join(df_age) %>%
  left_join(df_temp_name)

#### plot
col_scale <- c("ISIMIP" = "#F8766D", "ORAS5" = "#00BFC4", "NEMO-MEDUSA" = "#7CAE00")

p1 <- ggplot(data = pred %>% filter(age %in% c(1))) + 
  geom_line(aes(x = c.temp_within, 
                y = exp(fit),
                color = source_name),
            linewidth = 1) +
  geom_ribbon(aes(x = c.temp_within,
                  ymin = exp(lower),
                  ymax = exp(upper),
                  fill = source_name), 
              alpha = 0.1) +
  facet_grid(. ~ age_name)  +
  labs(x = expression('T'['individual-anomaly'] * ' (°C)'),
       y = "Predicted increment growth (μm)",
       color = "Dataset",
       fill = "Dataset") +
  scale_color_manual(values = col_scale) +
  scale_fill_manual(values = col_scale) 

p2 <- ggplot(data = pred %>% filter(age %in% c(5))) + 
  geom_line(aes(x = c.temp_within, 
                y = exp(fit),
                color = source_name),
            linewidth = 1) +
  geom_ribbon(aes(x = c.temp_within,
                  ymin = exp(lower),
                  ymax = exp(upper),
                  fill = source_name), 
              alpha = 0.1) +
  facet_grid(. ~ age_name) +
  labs(x = expression('T'['individual-anomaly'] * ' (°C)'),
       y = "Predicted increment growth (μm)",
       color = "Dataset",
       fill = "Dataset") +
  scale_color_manual(values = col_scale) +
  scale_fill_manual(values = col_scale) 

list_plot[[5]] <- p1 
list_plot[[6]] <- p2 

### merge and save plot ----
# merge plot
(list_plot[[1]] + list_plot[[2]] + 
   list_plot[[3]] + list_plot[[4]] +
    list_plot[[5]] + list_plot[[6]]) +
  plot_layout(nrow = 3, axis_titles = "collect", guides = "collect") +
  plot_annotation(tag_levels = 'A') &  
  theme(legend.position = "bottom", 
        plot.tag.position  = c(0.92, 0.82), # c(0.95, 0.85)
        axis.title = element_text(size = 9),
        legend.key.size = unit(10, "pt"),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9))

# save file
ggsave(last_plot(), file = file.path(dir_report, "fig6_temp effects.pdf"),
       device = cairo_pdf,
       width = 9, height = 18,
       units = "cm")
ggsave(last_plot(), file = file.path(dir_report, "fig6_temp effects.tiff"),
       width = 9, height = 18,
       units = "cm",
       dpi = 1000)

## figure 7 - extrinsic effects ----
list_plot <- list()

### ssb.i ----
list_model <- list("isimip" = m5_isimip,
                   "oras5" = m5_oras5,
                   "nemo-medusa" = m5_nm
)
# 0.02-0.46 - nemo-medusa different range ssb.i
pred <- tibble()
for (i in 1:3) {
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

list_plot[[1]] <- last_plot()

### recruitment.i ----
list_model <- list("isimip" = m5_isimip,
                   "oras5" = m5_oras5,
                   "nemo-medusa" = m5_nm
)
pred <- tibble()
for (i in 1:3) {
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

list_plot[[2]] <- last_plot()

### f ----
list_model <- list(
  "oras5" = m5_oras5
)

pred <- tibble()

for (i in 1:1) {
  model <- list_model[[i]]
  data <- model@frame
  temp_source <- names(list_model[i])
  
  pred_temp <- as.data.frame(Effect(c("f"), 
                                    model, 
                                    xlevels = list(f = unique(data$f)))) %>%
    mutate(source = temp_source)
  
  pred <- bind_rows(pred, pred_temp)
}

pred <- pred %>% 
  left_join(df_temp_name)

#### plot
# add a point of ISIMIP, NEMO-MEDUSA to pred to merge legend 
pred_add <- pred %>% 
  slice(1:2) %>%
  mutate(source_name = c("ISIMIP", "NEMO-MEDUSA"))

pred <- bind_rows(pred, pred_add) %>%
  mutate(source_name = factor(source_name, levels = c("ISIMIP", "ORAS5", "NEMO-MEDUSA")))

col_scale <- c("ISIMIP" = "#F8766D", "ORAS5" = "#00BFC4", "NEMO-MEDUSA" = "#7CAE00")

ggplot(data = pred) + 
  geom_line(aes(x = f, 
                y = exp(fit),
                color = source_name),
            linewidth = 1) +
  geom_ribbon(aes(x = f,
                  ymin = exp(lower),
                  ymax = exp(upper),
                  fill = source_name), 
              alpha = 0.1) +
  labs(x = "Fishing mortality",
       y = "Predicted increment growth (μm)",
       color = "Dataset",
       fill = "Dataset") +
  scale_color_manual(values = col_scale) +
  scale_fill_manual(values = col_scale)

list_plot[[3]] <- last_plot()

### nutrient ----

list_model <- list("isimip" = m5_nu_isimip,
                   "oras5" = m5_nu_oras5
)
# 
pred <- tibble()
for (i in 1:2) {
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
# add a point of NEMO-MEDUSA to pred to merge legend 
pred_add <- pred %>% 
  slice(1) %>%
  mutate(source_name = c("NEMO-MEDUSA"))

pred <- bind_rows(pred, pred_add) %>%
  mutate(source_name = factor(source_name, levels = c("ISIMIP", "ORAS5", "NEMO-MEDUSA")))

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

list_plot[[4]] <- last_plot()

### merge and save plot ----
# merge plot
(list_plot[[1]] + list_plot[[2]] + 
   list_plot[[3]] + list_plot[[4]]) +
  plot_layout(nrow = 2, axis_titles = "collect", guides = "collect") +
  plot_annotation(tag_levels = 'A') &  
  theme(legend.position = "bottom", 
        plot.tag.position  = c(0.95, 0.9),
        axis.title = element_text(size = 9),
        legend.key.size = unit(10, "pt"),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9)
        )

# save file
ggsave(last_plot(), file = file.path(dir_report, "fig7_extrinsic effects.pdf"),
       device = cairo_pdf,
       width =  19, height = 12.6,
       units = "cm")
ggsave(last_plot(), file = file.path(dir_report, "fig7_extrinsic effects.tiff"),
       width =  19, height = 12.6,
       units = "cm",
       dpi = 1000)

## figure 8 - interaction c.temp_within vs env ----

list_plot <- list()

### ssb.i ----
# range: p25, 50, 75
df_age <- tibble(age = c(1,5),
                 age_name = c("Age 1", "Age 5"))

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
  filter(age %in% c(1,5), 
         c.temp_within == 0) %>%
  mutate(source_name = c("ORAS5", "ORAS5", "NEMO-MEDUSA", "NEMO-MEDUSA"))

pred <- bind_rows(pred, pred_add) %>%
  mutate(source_name = factor(source_name, levels = c("ISIMIP", "ORAS5", "NEMO-MEDUSA"))) %>%
  mutate(var = "Spawning stock biomass")

col_scale <- c("ISIMIP" = "#F8766D", "ORAS5" = "#00BFC4", "NEMO-MEDUSA" = "#7CAE00")

p1 <- ggplot(data = pred %>% filter(age %in% c(1))) + 
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

p2 <- ggplot(data = pred %>% filter(age %in% c(5))) + 
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

list_plot[[1]] <- p1
list_plot[[2]] <- p2

### recuitment.i ----
# recruitment - nemo-medusa
# range: p25, 50, 75
df_age <- tibble(age = c(1,5),
                 age_name = c("Age 1", "Age 5"))

model <- m5_nm
data <- model@frame

df_recruitment <- tibble(p_name = c("P25", "P75"),
                         recruitment.i = round(quantile(data$recruitment.i, probs = c(0.25, 0.75)),3)) 

pred <- as.data.frame(Effect(c("log.age", "c.temp_within", "recruitment.i"), 
                             model, 
                             xlevels = list(log.age = unique(data$log.age),
                                            c.temp_within = seq(-1, 1, 0.1),
                                            recruitment.i = df_recruitment$recruitment.i))) %>%
  mutate(age = round(exp(log.age),1)) %>%
  mutate(source = "nemo-medusa")

pred <- pred %>% 
  left_join(df_age) %>%
  left_join(df_temp_name) %>%
  left_join(df_recruitment)

#### plot
# add a point of ISIMIP, ORAS5 to pred to merge legend 
pred_add <- pred %>% 
  filter(age %in% c(1,5), 
         c.temp_within == 0) %>%
  mutate(source_name = c("ISIMIP", "ISIMIP", "ORAS5", "ORAS5"))

pred <- bind_rows(pred, pred_add) %>%
  mutate(source_name = factor(source_name, levels = c("ISIMIP", "ORAS5", "NEMO-MEDUSA"))) %>%
  mutate(var = "Recruitment")

col_scale <- c("ISIMIP" = "#F8766D", "ORAS5" = "#00BFC4", "NEMO-MEDUSA" = "#7CAE00")

p1 <- ggplot(data = pred %>% filter(age %in% c(1))) + 
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
  scale_fill_manual(values = col_scale)

p2 <- ggplot(data = pred %>% filter(age %in% c(5))) + 
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
  facet_grid(var ~ age_name) +
  labs(x = expression('T'['individual-anomaly'] * ' (°C)'),
       y = "Predicted increment growth (μm)",
       color = "Dataset",
       fill = "Dataset",
       linetype = "Variable value") +
  scale_color_manual(values = col_scale) +
  scale_fill_manual(values = col_scale)

list_plot[[3]] <- p1
list_plot[[4]] <- p2

### f ----
# f - oras5
# range: p25, 50, 75
df_age <- tibble(age = c(1,5),
                 age_name = c("Age 1", "Age 5"))

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
  filter(age %in% c(1,5), 
         c.temp_within == 0) %>%
  mutate(source_name = c("ISIMIP", "ISIMIP", "NEMO-MEDUSA", "NEMO-MEDUSA"))

pred <- bind_rows(pred, pred_add) %>%
  mutate(source_name = factor(source_name, levels = c("ISIMIP", "ORAS5", "NEMO-MEDUSA"))) %>%
  mutate(var = "Fishing mortality")

col_scale <- c("ISIMIP" = "#F8766D", "ORAS5" = "#00BFC4", "NEMO-MEDUSA" = "#7CAE00")

p1 <- ggplot(data = pred %>% filter(age %in% c(1))) + 
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
  facet_grid(. ~ age_name) +
  labs(x = expression('T'['individual-anomaly'] * ' (°C)'),
       y = "Predicted increment growth (μm)",
       color = "Dataset",
       fill = "Dataset",
       linetype = "Variable value") +
  scale_color_manual(values = col_scale) +
  scale_fill_manual(values = col_scale)

p2 <- ggplot(data = pred %>% filter(age %in% c(5))) + 
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
  facet_grid(var ~ age_name) +
  labs(x = expression('T'['individual-anomaly'] * ' (°C)'),
       y = "Predicted increment growth (μm)",
       color = "Dataset",
       fill = "Dataset",
       linetype = "Variable value") +
  scale_color_manual(values = col_scale) +
  scale_fill_manual(values = col_scale)

list_plot[[5]] <- p1
list_plot[[6]] <- p2

### merge and save plot ----
# merge plot
(list_plot[[1]] + list_plot[[2]] + 
   list_plot[[3]] + list_plot[[4]] +
   list_plot[[5]] + list_plot[[6]]) +
  plot_layout(nrow = 3, axis_titles = "collect", guides = "collect") +
  plot_annotation(tag_levels = 'A') &  
  theme(legend.position = "bottom", 
        legend.box = "vertical",
        axis.title = element_text(size = 9),
        plot.tag.position  = c(0.81, 0.82),
        #plot.tag = element_text(size = 9),
        legend.key.size = unit(10, "pt"),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        strip.text.y = element_text(size = 7))

# save file
ggsave(last_plot(), file = file.path(dir_report, "fig8_temp vs env.pdf"),
       device = cairo_pdf,
       width = 9, height = 18,
       units = "cm")
ggsave(last_plot(), file = file.path(dir_report, "fig8_temp vs env.tiff"),
       width = 9, height = 18,
       units = "cm",
       dpi = 1000)

## figure 9 - ratio variance ----
#### data 
list_source <- c("isimip", "oras5", "nemo-medusa")

list_pair <- tibble(pop_pair = c("4bc/8ab", "7a/8ab", "4bc/7a"),
                    pop_pair_name = factor(c("North Sea/Bay of Biscay", 
                                             "Irish Sea/Bay of Biscay", 
                                             "North Sea/Irish Sea"),
                                           levels = c("North Sea/Bay of Biscay", 
                                                      "Irish Sea/Bay of Biscay", 
                                                      "North Sea/Irish Sea")))

var_ratio <- tibble()
for (i in 1:3) {
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
             position = position_dodge(width=0.5),
             size = rel(0.8)) +
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
        strip.text = element_text(size = 5.5),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 9))

#### save plot
ggsave(last_plot(), file = file.path(dir_report, "fig9_variance ratio.pdf"),
       device = cairo_pdf,
       width = 9, height = 9,
       units = "cm")
ggsave(last_plot(), file = file.path(dir_report, "fig9_variance ratio.tiff"),
       width = 9, height = 9,
       units = "cm",
       dpi = 1000)

# 4. TABLE ----
## table 1 - otolith sampling summary (doc) ----
## table 2 - list of predictors (doc) ----
## table 3 - intrinsic, extrinsic model (without nutrient) ----

#### setup
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

#### ref:
#https://cran.r-project.org/web/packages/sjPlot/vignettes/tab_mixed.html
#ref: https://cran.r-project.org/web/packages/sjPlot/vignettes/table_css.html

# labels
pred_labels <- c(
  "Intercept",
  "Age",
  "Age at capture",
  "Reading Institute (WUR)",
  "Reading Institute (WUR) * Age",
  "Spawning Stock Biomass",
  "Recruitment",
  "Temperature_individual-anomaly",
  "Temperature_individual-average",
  "Temperature_individual-anomaly * Age",
  "Temperature_individual-average * Age",
  "Temperature_individual-anomaly * Spawning Stock Biomass",
  "Fishing mortality",
  "Temperature_individual-anomaly * Fishing mortality",
  "Temperature_individual-anomaly * Recruitment"
)

# summary table
tab_model(m3, m5_isimip, m5_oras5, m5_nm, show.se = TRUE, 
          show.ci = NULL, show.p = FALSE,
          dv.labels = c("Intrinsic model", 
                        "Individual-level extrinsic model (ISIMIP)", 
                        "Individual-level extrinsic model (ORAS5)", 
                        "Individual-level extrinsic model (NEMO-MEDUSA)"),
          order.terms = c(1,2,3,4,5,6,7,13,8,10,12,15,14,9,11),
          pred.labels = pred_labels,
          string.pred = "Fixed Effects",
          string.est = "Estimate (SE)",
          collapse.se = TRUE,
          show.icc = FALSE,
          digits.rsq = 2,
          CSS = css_list,
          file = file.path(dir_report, "table3_model_ext_summary.html")
)

#### note: 
# tab_model cannot customised the random effect names
# things to be changed manually after tab_model
# 1. names random effects
# 2. number of years (tab_model returns number of pop.year)
# get no. year instead of no. pop.year
year_all = as.numeric(sub(".*:", "", m3@frame$pop.year))
year_nm = as.numeric(sub(".*:", "", m4_nm@frame$pop.year))
tibble(data = c("full", "nemo-medusa"),
       n_year = c(length(unique(year_all)), length(unique(year_nm))))









