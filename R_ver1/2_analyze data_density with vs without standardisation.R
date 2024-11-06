# BEFORE: 
#        2_analyze data.R

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
library(broom)      # get model coef

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
dir_output <- "./output"

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
         recruitment.i = recruitment/area_km2)

# plaice stock assessment
data_ple <- read_rds(file.path(dir_ices, "ple_stock-assessment_2023.rds")) %>%
  rename(ssb_ple = SSB,
         recruitment_ple = recruitment) %>%
  select(pop, year, ssb_ple, recruitment_ple) %>%
  left_join(datras) %>%
  mutate(ssb_ple.i = ssb_ple/area_km2,
         recruitment_ple.i = recruitment_ple/area_km2)

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

# 3. ANALYZE DATA ----
## isimip ----
#### data
data <- data_otl %>% 
  left_join(filter(data_temp, source == "isimip"),
            by = join_by(pop, year)) %>%
  filter(is.na(c.temp) == F) %>%
  left_join(data_sol, by = join_by(pop, year))

#### best model from 2_analyze data
m4.reml <- lmer(log.increment ~ 1 + log.age*datasource + log.aac +
                  ssb.i + recruitment.i +
                  (1 | fishid) + (1 | pop) + (1 + log.age | pop.year), 
                data = data, REML = T)

df_best <- tidy(m4.reml) %>% 
  filter(term %in% c("ssb.i", "recruitment.i")) %>%
  mutate(pop = "all pop",
         source = unique(data$source))

#### density effect
list_pop <- unique(data$pop)
df <- tibble()

for (i in 1:3) {
  #data
  data_sub <- data %>% filter(pop == list_pop[i])

    if (length(unique(data_sub$datasource)) == 2) {
    
    # fit model
    m4 <- lmer(log.increment ~ 1 + log.age*datasource + log.aac +
                 ssb.i + recruitment.i +
                 (1 | fishid) + (1 + log.age | pop.year), 
               data = data_sub, REML = T)
    
    df_temp <- tidy(m4)
     } else {
    
       m4 <- lmer(log.increment ~ 1 + log.age + log.aac +
                    ssb.i + recruitment.i +
                    (1 | fishid) + (1 + log.age | pop.year), 
                  data = data_sub, REML = T)
       
    df_temp <- tidy(m4)
  }
  
  # add info
  df_temp <- df_temp %>% 
    filter(term %in% c("ssb.i", "recruitment.i")) %>%
    mutate(pop = list_pop[i],
           source = unique(data$source))
  
  df <- bind_rows(df, df_temp)

}

df_isimip <- df %>%
  bind_rows(df_best)

## oras5 ----
#### data
data <- data_otl %>% 
  left_join(filter(data_temp, source == "oras5"),
            by = join_by(pop, year)) %>%
  filter(is.na(c.temp) == F) %>%
  left_join(data_sol, by = join_by(pop, year))

#### best model from 2_analyze data
m4.reml <- lmer(log.increment ~ 1 + log.age*datasource + log.aac +
                  c.temp*log.age + ssb.i + recruitment.i +
                  (1 | fishid) + (1 | pop) + (1 + log.age | pop.year), 
                data = data, REML = T)

df_best <- tidy(m4.reml) %>% 
  filter(term %in% c("ssb.i", "recruitment.i")) %>%
  mutate(pop = "all pop",
         source = unique(data$source))

#### density effect
list_pop <- unique(data$pop)
df <- tibble()

for (i in 1:3) {
  #data
  data_sub <- data %>% filter(pop == list_pop[i])
  
  if (length(unique(data_sub$datasource)) == 2) {
    
    # fit model
    m4 <- lmer(log.increment ~ 1 + log.age*datasource + log.aac +
                 c.temp*log.age + ssb.i + recruitment.i +
                 (1 | fishid) + (1 + log.age | pop.year), 
               data = data_sub, REML = T,
               control = lmerControl(optimizer = "Nelder_Mead"))
    
    df_temp <- tidy(m4)
  } else {
    
    m4 <- lmer(log.increment ~ 1 + log.age + log.aac +
                 c.temp*log.age + ssb.i + recruitment.i +
                 (1 | fishid) + (1 + log.age | pop.year), 
               data = data_sub, REML = T)
    
    df_temp <- tidy(m4)
  }
  
  # add info
  df_temp <- df_temp %>% 
    filter(term %in% c("ssb.i", "recruitment.i")) %>%
    mutate(pop = list_pop[i],
           source = unique(data$source))
  
  df <- bind_rows(df, df_temp)
  
}

# 7a model is singular -> refit without log.age slope
data_sub <- data %>% filter(pop == "7a")
m4 <- lmer(log.increment ~ 1 + log.age + log.aac +
             c.temp*log.age + ssb.i + recruitment.i +
             (1 | fishid) + (1 | pop.year), 
           data = data_sub, REML = T)
df_temp <- tidy(m4) %>% 
  filter(term %in% c("ssb.i", "recruitment.i")) %>%
  mutate(pop = "7a",
         source = unique(data$source))

df <- df %>% 
  filter(pop != "7a") %>%
  bind_rows(df_temp)

df_oras5 <- df %>%
  bind_rows(df_best)

## nemo-medusa ----
#### data
data <- data_otl %>% 
  left_join(filter(data_temp, source == "nemo-medusa"),
            by = join_by(pop, year)) %>%
  filter(is.na(c.temp) == F) %>%
  left_join(data_sol, by = join_by(pop, year))

#### best model from 2_analyze data
m4.reml <- lmer(log.increment ~ 1 + log.age*datasource + log.aac +
                  c.temp*log.age + ssb.i + recruitment.i +
                  (1 | fishid) + (1 | pop) + (1 + log.age | pop.year), 
                data = data, REML = T)

df_best <- tidy(m4.reml) %>% 
  filter(term %in% c("ssb.i", "recruitment.i")) %>%
  mutate(pop = "all pop",
         source = unique(data$source))

#### density effect
list_pop <- unique(data$pop)
df <- tibble()

for (i in 1:3) {
  #data
  data_sub <- data %>% filter(pop == list_pop[i])
  
  if (length(unique(data_sub$datasource)) == 2) {
    
    # fit model
    m4 <- lmer(log.increment ~ 1 + log.age*datasource + log.aac +
                 c.temp*log.age + ssb.i + recruitment.i +
                 (1 | fishid) + (1 + log.age | pop.year), 
               data = data_sub, REML = T,
               control = lmerControl(optimizer = "Nelder_Mead"))
    
    df_temp <- tidy(m4)
  } else {
    
    m4 <- lmer(log.increment ~ 1 + log.age + log.aac +
                 c.temp*log.age + ssb.i + recruitment.i +
                 (1 | fishid) + (1 + log.age | pop.year), 
               data = data_sub, REML = T)
    
    df_temp <- tidy(m4)
  }
  
  # add info
  df_temp <- df_temp %>% 
    filter(term %in% c("ssb.i", "recruitment.i")) %>%
    mutate(pop = list_pop[i],
           source = unique(data$source))
  
  df <- bind_rows(df, df_temp)
  
}

df_nm <- df %>%
  bind_rows(df_best)

# 4. SUMMARY ----
#### setup
df_data <- tibble(pop = c("all pop", "4bc", "7a", "8ab"),
                  data_name = c("Full data", "North Sea subset", "Irish Sea subset", "Bay of Biscay subset")) %>%
  mutate(data_name = factor(data_name, levels = data_name))

df_term <- tibble(term = c("ssb.i", "recruitment.i"),
                  term_name = c("Spawning stock biomass", "Recruitment")) %>%
  mutate(term_name = factor(term_name, levels = term_name))

df_temp_name <- tibble(source = c("isimip", "oras5", "nemo-medusa"),
                       source_name = c("ISIMIP", "ORAS5", "NEMO-MEDUSA")) %>%
  mutate(source_name = factor(source_name, levels = source_name))
                       
df_all <- bind_rows(df_isimip, df_oras5, df_nm) %>%
  left_join(df_data) %>%
  left_join(df_term) %>%
  left_join(df_temp_name)
  
#### save file
write_rds(df_all, file.path(dir_output, "robustness_density without standardisation.rds"))

#### plot
ggplot(data = df_all, 
       aes(x = estimate, 
           y = data_name)) +
  geom_point() +
  geom_linerange(aes(xmin = estimate - std.error,
                     xmax = estimate + std.error)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_grid(term_name ~ source_name)

ggplot(data = df_all, 
       aes(x = estimate, 
           y = data_name)) +
  geom_point() +
  geom_linerange(aes(xmin = estimate - std.error,
                     xmax = estimate + std.error)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_grid(source_name ~ term_name) +
  labs(x = "Parameter estimate",
       y = NULL)

