# BEFORE: 
#        1_explore data.R

# AFTER: 
#        output

# 1. SETUP ----------------------------------------------------------------------
# PACKAGES
library(tidyverse)  # process data frame
library(lme4)       # fit model
library(AICcmodavg) # calculate AIC
library(effects)    # display effect
library(MuMIn)      # compare models
library(sf)         # process geospatial data

# FUNCTION
# centring function
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

# dir to save output files
dir_output <- "./output_revised"

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

# 3. ANALYZE DATA ------------------------------------------------------------

## 3.1. Intrinsic structure --------------------------------------------

data <- data_otl %>% 
  mutate(method2 = if_else(method == "broken/burned", "broken", "sectioned"))

### 3.1.1. Random structure ---------

# Workflow
# 1. fit maximal intrinsic structure with increasing complex random effect structure (REML = T) (Morrongiello et al., 2019)
# 2. compare models using AICc

#### fishid
m1a <- lmer(log.increment ~ 1 + log.age + log.aac + log.age*pop + 
              (1 | fishid), data = data, REML = T)
m1b <- lmer(log.increment ~ 1 + log.age + log.aac + log.age*pop +
              (1 + log.age | fishid), data = data, REML = T) #singular
# m1b is singular -> continue with m1a

#### pop + pop.year (year nested in population) 
m2a <- lmer(log.increment ~ 1 + log.age + log.aac + log.age*pop +  
              (1 | fishid) + (1 | pop.year), data = data, REML = T)
m2b <- lmer(log.increment ~ 1 + log.age + log.aac + log.age*pop + 
              (1 | fishid) + (1 + log.age | pop.year), data = data, REML = T)

#### pop + pop.cohort (year nested in population) 
m2c <- lmer(log.increment ~ 1 + log.age + log.aac + log.age*pop + 
              (1 | fishid) + (1 | pop.cohort), data = data, REML = T)
m2d <- lmer(log.increment ~ 1 + log.age + log.aac + log.age*pop + 
              (1 | fishid) + (1 + log.age | pop.cohort), data = data, REML = T)

# m2d singular -> not compared
models <- list(m1a, m2a, m2b, m2c)
Modnames <- c('m1a', 'm2a', 'm2b', 'm2c')
aictab(cand.set = models, modnames = Modnames, sort = TRUE)

# m2b is the best model

#### pop + pop.year + pop.cohort
m2e <- lmer(log.increment ~ 1 + log.age + log.aac + log.age*pop +  
              (1 | fishid) + (1 | pop.year) + (1 | pop.cohort), data = data, REML = T)
m2f <- lmer(log.increment ~ 1 + log.age + log.aac + log.age*pop +  
              (1 | fishid) + (1 + log.age | pop.year) + (1 | pop.cohort), data = data, REML = T)
m2g <- lmer(log.increment ~ 1 + log.age + log.aac + log.age*pop +  
              (1 | fishid) + (1 | pop.year) + (1 + log.age | pop.cohort), data = data, REML = T)
m2h <- lmer(log.increment ~ 1 + log.age + log.aac + log.age*pop +  
              (1 | fishid) + (1 + log.age | pop.year) + (1 + log.age | pop.cohort), data = data, REML = T)

# m2h singular -> not compared
models <- list(m2b, m2e, m2f, m2g)
Modnames <- c('m2b', 'm2e', 'm2f', 'm2g')
aictab(cand.set = models, modnames = Modnames, sort = TRUE)

# m2f is the best model

# BEST RANDOM STRUCTURE: m2b - (1 | fishid) + (1 | pop) + (1 + log.age | pop.year)

### 3.1.2. Intrinsic structure ----------

#### Workflow
# 1. fit maximal intrinsic structure (with optimal random structure) (REML = F - using ML) (Morrongiello et al., 2019)
# 2. compare models using MuMIn dredge() (setting: na.action = "na.fail", REML = F - to compare fixed effects)
# 3. refit using REML = T for optimal structure (least variables and lowest AICc, consider more variables if delta_AICc >= 2/additional variable)  

#### fit maximal model
m3 <- lmer(log.increment ~ 1 + log.age + log.aac + log.age*pop +  
             (1 | fishid) + (1 + log.age | pop.year) + (1 | pop.cohort), data = data, 
           REML = F, na.action = "na.fail")

#### compare model
model.compare.intrinsic <- dredge(m3,  fixed = c("log.age"))

model.compare.intrinsic <- as.data.frame(model.compare.intrinsic) %>%
  mutate(var_added = 4 - rowSums(is.na(.)),
         id = seq(1,nrow(model.compare.intrinsic)))

write_rds(model.compare.intrinsic, file.path(dir_output, "intrinsic.model_comparison.rds"))

#### OPTIMAL INTRINSIC STRUCTURE: log.age + log.aac + 
m3.reml <- lmer(log.increment ~ 1 + log.age + log.aac +
                 (1 | fishid) + (1 + log.age | pop.year) + (1 | pop.cohort), data = data, REML = T)
write_rds(m3.reml, file.path(dir_output, "intrinsic.model_best.rds"))

### OPTIMAL INTRINSIC STRUCTURE - scaled variables
data <- data %>%
  mutate(s.log.age = s.(log.age),
         s.log.aac = s.(log.aac))

m3.reml_s <- lmer(log.increment ~ 1 + s.log.age + s.log.aac +
                    (1 | fishid) + (1 + s.log.age | pop.year) + (1 | pop.cohort), data = data, REML = T)
write_rds(m3.reml_s, file.path(dir_output, "intrinsic.model_best_scaled.rds"))

## 3.2. Extrinsic structure --------------------------------------------

#### Workflow
# for each temperature, fith data using the following workflow
# 1. fit maximal extirinsic structure (with optimal intrinsic structure) (REML = F - using ML) (Morrongiello et al., 2019)
# 2. compare models using MuMIn dredge() (setting: na.action = "na.fail", REML = F - to compare fixed effects)
# 3. refit using REML = T for optimal structure (least variables and lowest AICc, consider more variables if delta_AICc >= 2/additional variable)  

### 3.2.1. without nutrient data ----
#### isimip ----
#### data
data <- data_otl %>% 
  left_join(filter(data_temp, source == "isimip"),
            by = join_by(pop, year)) %>%
  filter(is.na(c.temp) == F) %>%
  left_join(data_sol, by = join_by(pop, year)) %>%
  mutate(method2 = if_else(method == "broken/burned", "broken", "sectioned"))

#### fit maximal model
m4 <- lmer(log.increment ~ 1 + log.age + log.aac +
             c.temp*log.age + 
             c.temp*ave.temp + 
             c.temp*f + c.temp*ssb.i + c.temp*recruitment.i +
             ssb.i*log.age + recruitment.i*log.age +
             (1 | fishid) + (1 + log.age | pop.year) + (1 | pop.cohort), 
           data = data, REML = F, na.action = "na.fail")

#### compare model
## note: keep intrinsic structure using subset argument
## the variables in `` should be in alphabet order, i.e., c.log.age:datasource instead of datasource:c.log.age
model.compare.extrinsic <- dredge(m4, subset = `log.age` & `log.aac`)

# add source, and notation of null model, number of var_added, and model id
model.compare.extrinsic <- as.data.frame(model.compare.extrinsic) %>% 
  mutate(source = unique(data$source)) %>%
  mutate(year_range = paste(min(data$year), "-", max(data$year)))
model.compare.extrinsic <- model.compare.extrinsic %>%
  mutate(var_added = 12 - rowSums(is.na(.)),
         id = seq(1,nrow(model.compare.extrinsic)))

# save model comparison
write_rds(model.compare.extrinsic, file.path(dir_output, paste0("extrinsic.model_comparison_", unique(data$source), ".rds")))

#### OPTIMAL EXTRINSIC STRUCTURE: c.temp + ave.temp + c.temp*ssb.i + recruitment.i (model id 1)
m4.reml <- lmer(log.increment ~ 1 + log.age + log.aac + 
                  c.temp +
                  ave.temp + 
                  c.temp*ssb.i + recruitment.i +
                  (1 | fishid) + (1 + log.age | pop.year) + (1 | pop.cohort), 
                data = data, REML = T)
write_rds(m4.reml, file.path(dir_output, paste0("extrinsic.model_best_", unique(data$source), ".rds")))

### OPTIMAL EXTRINSIC STRUCTURE - scaled variables
data <- data %>%
  mutate(s.log.age = s.(log.age),
         s.log.aac = s.(log.aac),
         s.c.temp = s.(c.temp),
         s.ave.temp = s.(ave.temp),
         s.ssb.i = s.(ssb.i),
         s.recruitment.i = s.(recruitment.i))

m4.reml_s <- lmer(log.increment ~ 1 + s.log.age + s.log.aac + 
                    s.c.temp +
                    s.ave.temp + 
                    s.c.temp*s.ssb.i + s.recruitment.i +
                    (1 | fishid) + (1 + log.age | pop.year) + (1 | pop.cohort), 
                  data = data, REML = T,
                  control = lmerControl(optimizer = "bobyqa")) 
# note: model failed to converge, but converge with optimizer "bobyqa" with very similar coef
write_rds(m4.reml_s, file.path(dir_output, paste0("extrinsic.model_best_scaled_", unique(data$source), ".rds")))

#### oras5 ----
#### data
data <- data_otl %>% 
  left_join(filter(data_temp, source == "oras5"),
            by = join_by(pop, year)) %>%
  filter(is.na(c.temp) == F) %>%
  left_join(data_sol, by = join_by(pop, year)) %>%
  mutate(method2 = if_else(method == "broken/burned", "broken", "sectioned"))

#### fit maximal model
m4 <- lmer(log.increment ~ 1 + log.age + log.aac + 
             c.temp*log.age + 
             c.temp*ave.temp + 
             c.temp*f + c.temp*ssb.i + c.temp*recruitment.i +
             ssb.i*log.age + recruitment.i*log.age +
             (1 | fishid) + (1 + log.age | pop.year) + (1 | pop.cohort), 
           data = data, REML = F, na.action = "na.fail")

#### compare model
## note: keep intrinsic structure using subset argument
## the variables in `` should be in alphabet order, i.e., c.log.age:datasource instead of datasource:c.log.age
model.compare.extrinsic <- dredge(m4, subset = `log.age` & `log.aac`)

# add source, and notation of null model, number of var_added, and model id
model.compare.extrinsic <- as.data.frame(model.compare.extrinsic) %>% 
  mutate(source = unique(data$source)) %>%
  mutate(year_range = paste(min(data$year), "-", max(data$year)))
model.compare.extrinsic <- model.compare.extrinsic %>%
  mutate(var_added = 12 - rowSums(is.na(.)),
         id = seq(1,nrow(model.compare.extrinsic)))

# save model comparison
write_rds(model.compare.extrinsic, file.path(dir_output, paste0("extrinsic.model_comparison_", unique(data$source), ".rds")))

#### OPTIMAL EXTRINSIC STRUCTURE:  ssb.i + recruitment.i (model id 65) - (candidates: id 65, 44, 12, 2, 1)
m4.reml <- lmer(log.increment ~ 1 + log.age + log.aac +
                  ssb.i + recruitment.i +
                  (1 | fishid) + (1 + log.age | pop.year) + (1 | pop.cohort), 
                data = data, REML = T)

# note: optimal sturcture - without random c.temp | pop
write_rds(m4.reml, file.path(dir_output, paste0("extrinsic.model_best_", unique(data$source), ".rds")))

# save model with ave.temp + c.temp for visualisation
m4.reml_temp <- lmer(log.increment ~ 1 + log.age + log.aac + 
                  ssb.i + recruitment.i +
                    ave.temp + c.temp + 
                  (1 | fishid) + (1 + log.age | pop.year) + (1 | pop.cohort), 
                data = data, REML = T)
write_rds(m4.reml_temp, file.path(dir_output, paste0("extrinsic.model_best_", unique(data$source), "_temp.rds")))

### OPTIMAL EXTRINSIC STRUCTURE - scaled variables
data <- data %>%
  mutate(s.log.age = s.(log.age),
         s.log.aac = s.(log.aac),
         s.c.temp = s.(c.temp),
         s.ssb.i = s.(ssb.i),
         s.recruitment.i = s.(recruitment.i))

m4.reml_s <- lmer(log.increment ~ 1 + s.log.age + s.log.aac + 
                    s.ssb.i + s.recruitment.i +
                    (1 | fishid) + (1 + s.log.age | pop.year) + (1 | pop.cohort), 
                  data = data, REML = T)
write_rds(m4.reml_s, file.path(dir_output, paste0("extrinsic.model_best_scaled_", unique(data$source), ".rds")))

### 3.2.2 with nutrient data ----
#### isimip ----
#### data
data <- data_otl %>% 
  left_join(filter(data_temp, source == "isimip"),
            by = join_by(pop, year)) %>%
  filter(is.na(c.temp) == F) %>%
  left_join(data_sol, by = join_by(pop, year)) %>%
  left_join(data_nu, by = join_by(pop, year)) %>%
  mutate(method2 = if_else(method == "broken/burned", "broken", "sectioned"))
data <- data %>% filter(!is.na(TN) & !is.na(TP))

#### fit maximal model - m4.reml + TN/TP
m4_nu <- lmer(log.increment ~ 1 + log.age + log.aac + 
                c.temp +
                #ave.temp + # only Nort Sea population -> no ave.temp
                c.temp*ssb.i + recruitment.i +
                TN + TP +
                (1 | fishid) + (1 + log.age | pop.year) + (1 | pop.cohort), 
              data = data, REML = F, na.action = "na.fail")

# unrealistic estimate of log.age vs intrinsic model -> fit without pop.cohort
# similar estimate of all variables to extrinsic model with full data
m4_nu <- lmer(log.increment ~ 1 + log.age + log.aac + 
                c.temp +
                #ave.temp + # only Nort Sea population -> no ave.temp
                c.temp*ssb.i + recruitment.i +
                TN + TP +
                (1 | fishid) + (1 + log.age | pop.year), 
              data = data, REML = F, na.action = "na.fail")

m4 <- read_rds(file.path(dir_output, "extrinsic.model_best_isimip.rds"))
fixef(m4)
fixef(m4_nu)

#### compare model
## note: keep intrinsic structure using subset argument
## the variables in `` should be in alphabet order, i.e., c.log.age:datasource instead of datasource:c.log.age
model.compare.extrinsic.nu <- dredge(m4_nu, subset = `log.age` & `log.aac` &
                                       c.temp & `c.temp:ssb.i` & recruitment.i &
                                        !(TN && TP))

# save model comparison
write_rds(model.compare.extrinsic.nu, file.path(dir_output, paste0("extrinsic.model.nu_comparison_", unique(data$source), ".rds")))

#### OPTIMAL EXTRINSIC STRUCTURE:  TP
m4_nu.reml <- lmer(log.increment ~ 1 + log.age + log.aac + 
                     c.temp +
                     #ave.temp + 
                     c.temp*ssb.i + recruitment.i +
                     TP +
                     (1 | fishid) + (1 + log.age | pop.year),
                   data = data, REML = T)

write_rds(m4_nu.reml, file.path(dir_output, paste0("extrinsic.model.nu_best_", unique(data$source), ".rds")))

### OPTIMAL EXTRINSIC STRUCTURE - scaled variables
data <- data %>%
  mutate(s.log.age = s.(log.age),
         s.log.aac = s.(log.aac),
         s.c.temp = s.(c.temp),
         #s.ave.temp = s.(ave.temp),
         s.ssb.i = s.(ssb.i),
         s.recruitment.i = s.(recruitment.i),
         s.TP = s.(TP))

m4_nu.reml_s <- lmer(log.increment ~ 1 + s.log.age + s.log.aac + 
                       s.c.temp +
                       #s.ave.temp + 
                       s.c.temp*ssb.i + s.recruitment.i +
                       s.TP +
                       (1 | fishid) + (1 + s.log.age | pop.year),
                     data = data, REML = T)
write_rds(m4_nu.reml_s, file.path(dir_output, paste0("extrinsic.model.nu_best_scaled_", unique(data$source), ".rds")))

#### oras5 ----
#### data
data <- data_otl %>% 
  left_join(filter(data_temp, source == "oras5"),
            by = join_by(pop, year)) %>%
  filter(is.na(c.temp) == F) %>%
  left_join(data_sol, by = join_by(pop, year)) %>%
  left_join(data_nu, by = join_by(pop, year)) %>%
  mutate(method2 = if_else(method == "broken/burned", "broken", "sectioned"))
data <- data %>% filter(!is.na(TN) & !is.na(TP))

#### fit maximal model - m4.reml + TN/TP
m4_nu <- lmer(log.increment ~ 1 + log.age + log.aac +
                ssb.i + recruitment.i +
                TN + TP +
                (1 | fishid) + (1 + log.age | pop.year) + (1 | pop.cohort), 
              data = data, REML = F, na.action = "na.fail")

# unrealistic estimate of log.age vs intrinsic model -> fit without pop.cohort
# similar estimate of all variables to extrinsic model with full data
m4_nu <- lmer(log.increment ~ 1 + log.age + log.aac + 
                ssb.i + recruitment.i +
                TN + TP +
                (1 | fishid) + (1 + log.age | pop.year), 
              data = data, REML = F, na.action = "na.fail")

m4 <- read_rds(file.path(dir_output, "extrinsic.model_best_oras5.rds"))
fixef(m4)
fixef(m4_nu)

#### compare model
## note: keep intrinsic structure using subset argument
## the variables in `` should be in alphabet order, i.e., c.log.age:datasource instead of datasource:c.log.age
model.compare.extrinsic.nu <- dredge(m4_nu, subset = `log.age` & `log.aac` &
                                       ssb.i & recruitment.i &
                                       !(TN && TP))

# save model comparison
write_rds(model.compare.extrinsic.nu, file.path(dir_output, paste0("extrinsic.model.nu_comparison_", unique(data$source), ".rds")))

#### OPTIMAL EXTRINSIC STRUCTURE:  TP
m4_nu.reml <- lmer(log.increment ~ 1 + log.age + log.aac + 
                     ssb.i + recruitment.i +
                     TP +
                     (1 | fishid) + (1 + log.age | pop.year), 
                   data = data, REML = T)

write_rds(m4_nu.reml, file.path(dir_output, paste0("extrinsic.model.nu_best_", unique(data$source), ".rds")))

### OPTIMAL EXTRINSIC STRUCTURE - scaled variables
data <- data %>%
  mutate(s.log.age = s.(log.age),
         s.log.aac = s.(log.aac),
         s.c.temp = s.(c.temp),
         s.ssb.i = s.(ssb.i),
         s.recruitment.i = s.(recruitment.i),
         s.TP = s.(TP))

m4_nu.reml_s <- lmer(log.increment ~ 1 + s.log.age + s.log.aac + 
                       s.ssb.i + s.recruitment.i +
                       s.TP +
                       (1 | fishid) + (1 + s.log.age | pop.year), 
                     data = data, REML = T)
write_rds(m4_nu.reml_s, file.path(dir_output, paste0("extrinsic.model.nu_best_scaled_", unique(data$source), ".rds")))

## 3.3. Extrinsic structure extended ----
### 3.3.1. without nutrient data ----
#### isimip ----
#### data
data <- data_otl %>% 
  left_join(filter(data_temp, source == "isimip"),
            by = join_by(pop, year)) %>%
  filter(is.na(c.temp) == F) %>%
  left_join(data_sol, by = join_by(pop, year)) %>%
  mutate(method2 = if_else(method == "broken/burned", "broken", "sectioned"))

data <- data %>% 
  group_by(fishid) %>%
  mutate(c.temp_between = mean(c.temp),
         c.temp_within = c.temp - c.temp_between) %>%
  ungroup()

#### difference c.temp_within and c.temp_between (see van de Pol and Wright 2009)
# add c.temp and c.temp_between to the best m4 model
m <- lmer(log.increment ~ 1 + log.age + log.aac + 
            c.temp +
            ave.temp + 
            c.temp*ssb.i + recruitment.i +
            c.temp_between + 
            (1 | fishid) + (1 + log.age | pop.year) + (1 | pop.cohort), 
          data = data, REML = T)
summary(m) #c.temp_between 0.22 +- 0.03 -> c.temp_within differs c.temp_between

#### fit maximal model
m5 <- lmer(log.increment ~ 1 + log.age + log.aac + 
             c.temp_within*log.age + c.temp_between*log.age + 
             c.temp_within*c.temp_between +
             c.temp_within*ave.temp +
             c.temp_within*f + c.temp_within*ssb.i + c.temp_within*recruitment.i +
             (1 | fishid) + (1 + log.age | pop.year) + (1 | pop.cohort), 
           data = data, REML = F, na.action = "na.fail")

#### compare model
## note: keep intrinsic structure using subset argument
## the variables in `` should be in alphabet order, i.e., c.log.age:datasource instead of datasource:c.log.age
model.compare.extrinsic <- dredge(m5, subset = `log.age` & `log.aac` & 
                                    ave.temp &
                                    c.temp_within & c.temp_between & ssb.i & recruitment.i)

# add source, and notation of null model, number of var_added, and model id
model.compare.extrinsic <- as.data.frame(model.compare.extrinsic) %>% 
  mutate(source = unique(data$source)) %>%
  mutate(year_range = paste(min(data$year), "-", max(data$year)))

model.compare.extrinsic <- model.compare.extrinsic %>%
  mutate(var_added = 13 - rowSums(is.na(.)), #12 with c.temp_between; 13 with c.temp_between*log.age
         id = seq(1,nrow(model.compare.extrinsic)))

# save model comparison
write_rds(model.compare.extrinsic, file.path(dir_output, paste0("extrinsic.model.ext_comparison_", unique(data$source), ".rds")))

#### OPTIMAL EXTRINSIC EXTENDED STRUCTURE - FIXED EFFECT: 
# c.temp_within*log.age + c.temp_between*log.age + ave.temp + c.temp_within*ssb.i + recruitment.i (model id 1)

# compare extrinsic extended structure with extended structure
m4 <- lmer(log.increment ~ 1 + log.age + log.aac + 
             c.temp +
             ave.temp + 
             c.temp*ssb.i + recruitment.i +
             (1 | fishid) + (1 + log.age | pop.year) + (1 | pop.cohort), 
           data = data, REML = F)
m5 <- lmer(log.increment ~ 1 + log.age + log.aac +
             c.temp_within*log.age + c.temp_between*log.age + 
             ave.temp +
             c.temp_within*ssb.i + recruitment.i +
             (1 | fishid) + (1 + log.age | pop.year) + (1 | pop.cohort), 
           data = data, REML = F)
# note: 
# m5 failed to converge with default optimizer
# model comparison
aictab(list("m4" = m4, "m5" = m5))

# note: extrinsic extended structure significantly improve the AICc

#### OPTIMAL EXTRINSIC EXTENDED STRUCTURE - RANDOM EFFECT: 
m5a <- lmer(log.increment ~ 1 + log.age + log.aac + 
              c.temp_within*log.age + c.temp_between*log.age + 
              ave.temp +
              c.temp_within*ssb.i + recruitment.i +
              (1 | fishid) + (1 + log.age | pop.year) + (1 | pop.cohort), 
            data = data, REML = T)

m5b <- lmer(log.increment ~ 1 + log.age + log.aac + 
              c.temp_within*log.age + c.temp_between*log.age + 
              ave.temp +
              c.temp_within*ssb.i + recruitment.i +
              (1 + c.temp_within | fishid) + (1 + log.age | pop.year) + (1 | pop.cohort), 
            data = data, REML = T, 
            control = lmerControl(optimizer = "bobyqa")) 

# note: 
# m5b failed to converge with default optimizer
# model comparison
aictab(list("m5a" = m5a, "m5b" = m5b))

#### OPTIMAL EXTRINSIC EXTENDED STRUCTURE - m5b
m5.reml <- lmer(log.increment ~ 1 + log.age + log.aac + 
                  c.temp_within*log.age + c.temp_between*log.age + 
                  ave.temp +
                  c.temp_within*ssb.i + recruitment.i +
                  (1 + c.temp_within | fishid) + (1 + log.age | pop.year) + (1 | pop.cohort), 
                data = data, REML = T, 
                control = lmerControl(optimizer = "bobyqa"))
write_rds(m5.reml, file.path(dir_output, paste0("extrinsic.model.ext_best_", unique(data$source), ".rds")))

### OPTIMAL EXTRINSIC STRUCTURE - scaled variables
data <- data %>%
  mutate(s.log.age = s.(log.age),
         s.log.aac = s.(log.aac),
         s.c.temp_within = s.(c.temp_within),
         s.c.temp_between = s.(c.temp_between),
         s.ave.temp = s.(ave.temp),
         s.ssb.i = s.(ssb.i),
         s.recruitment.i = s.(recruitment.i))

m5.reml_s <- lmer(log.increment ~ 1 + s.log.age + s.log.aac +
                    s.c.temp_within*s.log.age + s.c.temp_between*s.log.age + 
                    s.ave.temp +
                    s.c.temp_within*s.ssb.i + s.recruitment.i +
                    (1 + s.c.temp_within | fishid) + (1 + s.log.age | pop.year) + (1 | pop.cohort), 
                  data = data, REML = T, 
                  control = lmerControl(optimizer = "bobyqa"))

write_rds(m5.reml_s, file.path(dir_output, paste0("extrinsic.model.ext_best_scaled_", unique(data$source), ".rds")))

#### oras5 ----
#### data
data <- data_otl %>% 
  left_join(filter(data_temp, source == "oras5"),
            by = join_by(pop, year)) %>%
  filter(is.na(c.temp) == F) %>%
  left_join(data_sol, by = join_by(pop, year)) %>%
  mutate(method2 = if_else(method == "broken/burned", "broken", "sectioned"))

data <- data %>% 
  group_by(fishid) %>%
  mutate(c.temp_between = mean(c.temp),
         c.temp_within = c.temp - c.temp_between) %>%
  ungroup()

#### difference c.temp_within and c.temp_between (see van de Pol and Wright 2009)
# add c.temp and c.temp_between to the best m4 model
m <- lmer(log.increment ~ 1 + log.age + log.aac + 
            ssb.i + recruitment.i +
            c.temp*log.age + c.temp_between*log.age +
            (1 | fishid) + (1 + log.age | pop.year) + (1 | pop.cohort), 
          data = data, REML = T)
summary(m) 
# withou age: c.temp_between -0.17 +- 0.04 -> c.temp_within differ c.temp_between

#### fit maximal model
m5 <- lmer(log.increment ~ 1 + log.age + log.aac + 
             c.temp_within*log.age + c.temp_between*log.age + 
             c.temp_within*c.temp_between +
             c.temp_within*ave.temp +
             c.temp_within*f + c.temp_within*ssb.i + c.temp_within*recruitment.i +
             (1 | fishid) + (1 + log.age | pop.year) + (1 | pop.cohort), 
           data = data, REML = F, na.action = "na.fail")

#### compare model
## note: keep intrinsic structure using subset argument
## the variables in `` should be in alphabet order, i.e., c.log.age:datasource instead of datasource:c.log.age
model.compare.extrinsic <- dredge(m5, subset = `log.age` & `log.aac` &
                                    c.temp_within & c.temp_between & ssb.i & recruitment.i)

# add source, and notation of null model, number of var_added, and model id
model.compare.extrinsic <- as.data.frame(model.compare.extrinsic) %>% 
  mutate(source = unique(data$source)) %>%
  mutate(year_range = paste(min(data$year), "-", max(data$year)))

model.compare.extrinsic <- model.compare.extrinsic %>%
  mutate(var_added = 13 - rowSums(is.na(.)), #12 with c.temp_between; 13 with c.temp_between*log.age
         id = seq(1,nrow(model.compare.extrinsic)))

# save model comparison
write_rds(model.compare.extrinsic, file.path(dir_output, paste0("extrinsic.model.ext_comparison_", unique(data$source), ".rds")))

#### OPTIMAL EXTRINSIC EXTENDED STRUCTURE - FIXED EFFECT: (candidates id 14, 1)
# c.temp_within*log.age + c.temp_between + c.temp_within*c.temp_between + c.temp_within*f + ssb.i + c.temp_within*recruitment.i (model id 1)

# compare extrinsic extended structure with extended structure
m4 <- lmer(log.increment ~ 1 + log.age + log.aac + 
             ssb.i + recruitment.i +
             (1 | fishid) + (1 + log.age | pop.year) + (1 | pop.cohort), 
           data = data, REML = F)
m5 <- lmer(log.increment ~ 1 + log.age + log.aac +
             c.temp_within*log.age + c.temp_between + 
             c.temp_within*c.temp_between +
             c.temp_within*f + ssb.i + c.temp_within*recruitment.i +
             (1 | fishid) + (1 + log.age | pop.year) + (1 | pop.cohort), 
           data = data, REML = F)

# model comparison
aictab(list("m4" = m4, "m5" = m5))

# note: extrinsic extended structure significantly improve the AICc

#### OPTIMAL EXTRINSIC EXTENDED STRUCTURE - RANDOM EFFECT: 
m5a <- lmer(log.increment ~ 1 + log.age + log.aac +
              c.temp_within*log.age + c.temp_between +
              c.temp_within*c.temp_between +
              c.temp_within*f + ssb.i + c.temp_within*recruitment.i +
              (1 | fishid) + (1 + log.age | pop.year) + (1 | pop.cohort), 
            data = data, REML = T)

m5b <- lmer(log.increment ~ 1 + log.age + log.aac + 
              c.temp_within*log.age + c.temp_between + 
              c.temp_within*c.temp_between +
              c.temp_within*f + ssb.i + c.temp_within*recruitment.i +
              (1 + c.temp_within | fishid) + (1 + log.age | pop.year) + (1 | pop.cohort), 
            data = data, REML = T)

# note: 
# model comparison
aictab(list("m5a" = m5a, 
            "m5b" = m5b))

#### OPTIMAL EXTRINSIC EXTENDED STRUCTURE - m5b 
m5.reml <- lmer(log.increment ~ 1 + log.age + log.aac +
                  c.temp_within*log.age + c.temp_between + 
                  c.temp_within*c.temp_between +
                  c.temp_within*f + ssb.i + c.temp_within*recruitment.i +
                  (1 + c.temp_within | fishid) + (1 + log.age | pop.year) + (1 | pop.cohort), 
                data = data, REML = T)
write_rds(m5.reml, file.path(dir_output, paste0("extrinsic.model.ext_best_", unique(data$source), ".rds")))

### OPTIMAL EXTRINSIC STRUCTURE - scaled variables
data <- data %>%
  mutate(s.log.age = s.(log.age),
         s.log.aac = s.(log.aac),
         s.c.temp_within = s.(c.temp_within),
         s.c.temp_between = s.(c.temp_between),
         s.f = s.(f),
         s.ssb.i = s.(ssb.i),
         s.recruitment.i = s.(recruitment.i))

m5.reml_s <- lmer(log.increment ~ 1 + s.log.age + s.log.aac + method2 +
                    s.c.temp_within*s.log.age + s.c.temp_between + 
                    s.c.temp_within*s.c.temp_between +
                    s.c.temp_within*s.f + s.ssb.i + s.c.temp_within*s.recruitment.i +
                    (1 + s.c.temp_within | fishid) + (1 + s.log.age | pop.year) + (1 | pop.cohort), 
                  data = data, REML = T)
write_rds(m5.reml_s, file.path(dir_output, paste0("extrinsic.model.ext_best_scaled_", unique(data$source), ".rds")))

### 3.3.2. with nutrient data ----
#### isimip ----
#### data
data <- data_otl %>% 
  left_join(filter(data_temp, source == "isimip"),
            by = join_by(pop, year)) %>%
  filter(is.na(c.temp) == F) %>%
  left_join(data_sol, by = join_by(pop, year)) %>%
  left_join(data_nu, by = join_by(pop, year)) %>%
  mutate(method2 = if_else(method == "broken/burned", "broken", "sectioned"))

data <- data %>% filter(!is.na(TN) & !is.na(TP))

data <- data %>% 
  group_by(fishid) %>%
  mutate(c.temp_between = mean(c.temp),
         c.temp_within = c.temp - c.temp_between) %>%
  ungroup()

#### fit maximal model - m4.reml + TN/TP
m5_nu <- lmer(log.increment ~ 1 + log.age + log.aac + 
                c.temp_within*log.age + c.temp_between*log.age + 
                #ave.temp +
                c.temp_within*ssb.i + recruitment.i +
                TN + TP +
                (1 + c.temp_within| fishid) + (1 + log.age | pop.year) + (1 | pop.cohort), 
              data = data, REML = F, na.action = "na.fail",
              control = lmerControl(optimizer = "bobyqa"))

# unrealistic estimate of log.age vs intrinsic model -> fit without pop.cohort
# similar estimate of all variables to extrinsic model with full data
m5_nu <- lmer(log.increment ~ 1 + log.age + log.aac +
                c.temp_within*log.age + c.temp_between*log.age + 
                #ave.temp +
                c.temp_within*ssb.i + recruitment.i +
                TN + TP +
                (1 + c.temp_within | fishid) + (1 + log.age | pop.year), 
              data = data, REML = F, na.action = "na.fail")

m5 <- read_rds(file.path(dir_output, "extrinsic.model.ext_best_isimip.rds"))
fixef(m5)
fixef(m5_nu)

#### compare model
## note: keep intrinsic structure using subset argument
## the variables in `` should be in alphabet order, i.e., c.log.age:datasource instead of datasource:c.log.age
model.compare.extrinsic.nu <- dredge(m5_nu, subset = `log.age` & `log.aac` &
                                       `c.temp_within:log.age` & `c.temp_between:log.age` &
                                       `c.temp_within:ssb.i` & recruitment.i &
                                       !(TN && TP))

# save model comparison
write_rds(model.compare.extrinsic.nu, file.path(dir_output, paste0("extrinsic.model.ext.nu_comparison_", unique(data$source), ".rds")))

#### OPTIMAL EXTRINSIC STRUCTURE: + TP
m5_nu.reml <- lmer(log.increment ~ 1 + log.age + log.aac + 
                     c.temp_within*log.age + c.temp_between*log.age + 
                     #ave.temp +
                     c.temp_within*ssb.i + recruitment.i +
                     TP +
                     (1 + c.temp_within| fishid) + (1 + log.age | pop.year), 
                   data = data, REML = T)

write_rds(m5_nu.reml, file.path(dir_output, paste0("extrinsic.model.ext.nu_best_", unique(data$source), ".rds")))

### OPTIMAL EXTRINSIC STRUCTURE - scaled variables
data <- data %>%
  mutate(s.log.age = s.(log.age),
         s.log.aac = s.(log.aac),
         s.c.temp_within = s.(c.temp_within),
         s.c.temp_between = s.(c.temp_between),
         s.ssb.i = s.(ssb.i),
         s.recruitment.i = s.(recruitment.i),
         s.TP = s.(TP))

m5_nu.reml_s <- lmer(log.increment ~ 1 + s.log.age + s.log.aac + 
                       s.c.temp_within*s.log.age + s.c.temp_between*s.log.age + 
                       #ave.temp +
                       s.c.temp_within*s.ssb.i + s.recruitment.i +
                       s.TP +
                       (1 + s.c.temp_within| fishid) + (1 + s.log.age | pop.year), 
                     data = data, REML = T)
#control = lmerControl(optimizer = "bobyqa")

write_rds(m5_nu.reml_s, file.path(dir_output, paste0("extrinsic.model.ext.nu_best_scaled_", unique(data$source), ".rds")))

#### oras5 ----
#### data
data <- data_otl %>% 
  left_join(filter(data_temp, source == "oras5"),
            by = join_by(pop, year)) %>%
  filter(is.na(c.temp) == F) %>%
  left_join(data_sol, by = join_by(pop, year)) %>%
  left_join(data_nu, by = join_by(pop, year)) %>%
  mutate(method2 = if_else(method == "broken/burned", "broken", "sectioned"))

data <- data %>% filter(!is.na(TN) & !is.na(TP))

data <- data %>% 
  group_by(fishid) %>%
  mutate(c.temp_between = mean(c.temp),
         c.temp_within = c.temp - c.temp_between) %>%
  ungroup()

#### fit maximal model - m5.reml + TN/TP
m5_nu <- lmer(log.increment ~ 1 + log.age + log.aac + 
                c.temp_within*log.age + c.temp_between + 
                c.temp_within*c.temp_between +
                c.temp_within*f + ssb.i + c.temp_within*recruitment.i +
                TN + TP +
                (1 + c.temp_within | fishid) + (1 + log.age | pop.year) + (1 | pop.cohort), 
              data = data, REML = F, na.action = "na.fail",
              control = lmerControl(optimizer = "bobyqa"))

# unrealistic estimate of log.age vs intrinsic model -> fit without pop.cohort
# similar estimate of all variables to extrinsic model with full data
m5_nu <- lmer(log.increment ~ 1 + log.age + log.aac + 
                c.temp_within*log.age + c.temp_between + 
                c.temp_within*c.temp_between +
                c.temp_within*f + ssb.i + c.temp_within*recruitment.i +
                TN + TP +
                (1 + c.temp_within | fishid) + (1 + log.age | pop.year), 
              data = data, REML = F, na.action = "na.fail")

m5 <- read_rds(file.path(dir_output, "extrinsic.model.ext_best_oras5.rds"))
fixef(m5)
fixef(m5_nu)

#### compare model
## note: keep intrinsic structure using subset argument
## the variables in `` should be in alphabet order, i.e., c.log.age:datasource instead of datasource:c.log.age
model.compare.extrinsic.nu <- dredge(m5_nu, subset = `log.age` & `log.aac` &
                                       `c.temp_within:log.age` & `c.temp_between` & 
                                       `c.temp_between:c.temp_within` &
                                       `c.temp_within:f` & ssb.i & `c.temp_within:recruitment.i` &
                                       !(TN && TP))

# save model comparison
write_rds(model.compare.extrinsic.nu, file.path(dir_output, paste0("extrinsic.model.ext.nu_comparison_", unique(data$source), ".rds")))

#### OPTIMAL EXTRINSIC STRUCTURE:  + TP
m5_nu.reml <- lmer(log.increment ~ 1 + log.age + log.aac + 
                     c.temp_within*log.age + c.temp_between + 
                     c.temp_within*c.temp_between +
                     c.temp_within*f + ssb.i + c.temp_within*recruitment.i +
                     TP +
                     (1 + c.temp_within | fishid) + (1 + log.age | pop.year), 
                   data = data, REML = T,
                   control = lmerControl(optimizer = "bobyqa"))

write_rds(m5_nu.reml, file.path(dir_output, paste0("extrinsic.model.ext.nu_best_", unique(data$source), ".rds")))

### OPTIMAL EXTRINSIC STRUCTURE - scaled variables
data <- data %>%
  mutate(s.log.age = s.(log.age),
         s.log.aac = s.(log.aac),
         s.c.temp_within = s.(c.temp_within),
         s.c.temp_between = s.(c.temp_between),
         s.f = s.(f),
         s.ssb.i = s.(ssb.i),
         s.recruitment.i = s.(recruitment.i),
         s.TP = s.(TP))

m5_nu.reml_s <- lmer(log.increment ~ 1 + s.log.age + s.log.aac + method2 +
                       s.c.temp_within*s.log.age + s.c.temp_between + 
                       s.c.temp_within*s.c.temp_between +
                       s.c.temp_within*s.f + s.ssb.i + s.c.temp_within*s.recruitment.i +
                       s.TP +
                       (1 + c.temp_within | fishid) + (1 + s.log.age | pop.year), 
                     data = data, REML = T,
                     control = lmerControl(optimizer = "bobyqa"))

write_rds(m5_nu.reml_s, file.path(dir_output, paste0("extrinsic.model.ext.nu_best_scaled_", unique(data$source), ".rds")))
