# BEFORE: 
#        1_explore data.R
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
dir_output <- "./WP1/output"

# 2. LOAD DATA ------------------------------------------------------------

## Otolith data ---------------------------------------------------

dir_otl <- "./WP1/data"
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
dir_gis <- "./WP1/data/admin"
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
dir_temp <- "./WP1/data/temp"

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

dir_ices <- "./WP1/data/ices"

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
dir_nu <- "./WP1/data/nutrient"
data_nu <- read_rds(file.path(dir_nu, "ospar_subset_1978-2017_ices_4abc.rds"))
# summarize all river by year 
data_nu <- data_nu %>% 
  group_by(IcesArea, year) %>% 
  summarize(TN = sum(TN)/1000,
            TP = sum(TP)/1000) %>%
  mutate(IcesArea = "4bc") %>%
  rename(pop = IcesArea)
# unit: 1000tN/year, 1000tP/year

# 3. INDIVIDUAL PLASTICITY VARIANCE ----
## isimip ----
#### data
data <- data_otl %>% 
  left_join(filter(data_temp, source == "isimip"),
            by = join_by(pop, year)) %>%
  filter(is.na(c.temp) == F) %>%
  left_join(data_sol, by = join_by(pop, year))

#### model results
## read model
m <- read_rds(file.path(dir_output, paste0("extrinsic.model.ext_best_", unique(data$source), ".rds")))

## extract thermal reaction norm slope
ind_plastic <- ranef(m)$fishid
ind_plastic <- ind_plastic %>% mutate(fishid = row.names(ind_plastic),
                                      reaction_norm = c.temp_within)
### ratio variance ----
# get pop, fishid, and age
df_fishid <- data %>% select(pop, fishid, age, aac, cohort) %>% 
  group_by(pop, fishid, aac, cohort) %>%
  summarize(age = max(age))

# get ind_plastic for each pop
ind_plastic <- ind_plastic %>% left_join(df_fishid)

## visualisation of reaction norm variation vs age  
ggplot(data = ind_plastic, aes(x = reaction_norm, y = factor(age))) +
  geom_density_ridges(alpha = 0.5) + 
  facet_grid(~ pop) +
  theme_bw()

ggplot(data = ind_plastic, aes(x = age, y = reaction_norm)) +
  geom_jitter(alpha = 0.2) + 
  facet_grid(pop ~ .) +
  theme_bw()

# variation seems to increases with age

# different age group: 1-5, 6-11, 12+ (age >= 5 as suggested in Morrongiello 2019)
# as 8ab data only to age 11
ind_plastic <- ind_plastic %>% mutate(age_group = 
                                        if_else(age <= 5, "Age 2-5",
                                                if_else(age > 5 & age <= 11, "Age 6-11", "Age 12+"))) %>%
  mutate(age_group = factor(age_group, levels = c("Age 2-5", "Age 6-11", "Age 12+")))

ggplot(data = ind_plastic, aes(x = reaction_norm, fill = age_group)) +
  geom_density(alpha = 0.5) + 
  facet_grid(pop ~ .) 

#### variation ratio test
# conduct tests for 4 subsets of data
# age 2 - 11
# age 2 - 22 (all ages)
# age 6 - 11
# age 6 - 22 (all ages)

## setup
# pop_pair - variance ratio test: var_pop1/var_pop2
pop_pair <- tibble(pop1 = c("4bc", "4bc", "7a"),
                   pop2 = c("7a", "8ab", "8ab"))

# age range
range(ind_plastic$age)
age_range <- tibble(min_age = c(2, 2, 6, 6),
                    max_age = c(22, 11, 22, 11))

# note: for each pair, and each age range, 
# sample 100 times the pop with larger sample size to have 2 pops with the same sample size
# for each sample, do variance ratio test with 10000 bootstrap 

n_pop_pair <- nrow(pop_pair)
n_age_range <- nrow(age_range) 
n_sample <- 100
n_boot <- 10000

var_ratio <- tibble()
# set pair
for(p in 1:n_pop_pair) {
  
  print(paste0("processing ", pop_pair$pop1[p], "/", pop_pair$pop2[p]))
  
  pop1 <- pop_pair$pop1[p]
  pop2 <- pop_pair$pop2[p]
  
  df_pop1 <- ind_plastic %>% filter(pop == pop1)
  df_pop2 <- ind_plastic %>% filter(pop == pop2)
  
  # set test 
  for(t in 1:n_age_range) {
    
    print(paste0("processing ", age_range$min_age[t], "-", age_range$max_age[t]))
    
    min_age <- age_range$min_age[t]
    max_age <- age_range$max_age[t]
    
    df_pop1_age <- df_pop1 %>% filter(age >= min_age, age <= max_age)
    df_pop2_age <- df_pop2 %>% filter(age >= min_age, age <= max_age)
    
    # set sample and bootstrap
    results_sample <- tibble()
    for(i in 1:n_sample) {
      set.seed(i)
      
      # subset pop with larger sample size to have 2 pops with the same sample size 
      if (nrow(df_pop1_age) >= nrow(df_pop2_age)) {
        pop1_var <- sample(df_pop1_age$reaction_norm, nrow(df_pop2_age))
        pop2_var <- df_pop2_age$reaction_norm
      } else {
        pop1_var <- df_pop1_age$reaction_norm
        pop2_var <- sample(df_pop2_age$reaction_norm, nrow(df_pop1_age))
      }
      
      df_var <- tibble(pop1_var = pop1_var,
                       pop2_var = pop2_var)
      
      func_var.test <- function(data, i){
        data <- data[i, ]
        results  <- var.test(data$pop1_var, data$pop2_var)
        c(est = results$estimate, CI = results$conf.int, p.value = results$p.value)
      }
      
      boot <- boot(df_var, func_var.test, R = n_boot)
      colnames(boot$t) <- c("var_ratio", "ci_lower", "ci_upper", "p.value")
      var_ratio_temp <- as.data.frame(boot$t) %>% 
        summarize(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>%
        mutate(sample_id = i,
               n_pop1 = nrow(df_pop1_age),
               n_pop2 = nrow(df_pop2_age),
               age_range = paste0("Age ", min_age, " - ", "Age ", max_age),
               pop_pair = paste0(pop1, "/", pop2))
      
      var_ratio <- bind_rows(var_ratio, var_ratio_temp)
    }
  }
}

# save file
write_rds(var_ratio, file.path(dir_output, paste0("ind.plastic_var.ratio_", unique(data$source), ".rds")))

# summary
var_ratio_sum <- var_ratio %>%
  group_by(pop_pair, age_range, n_pop1, n_pop2) %>%
  summarise(across(where(is.numeric), ~ mean(.x))) %>% 
  mutate(across(var_ratio:ci_upper, ~ round(., 2))) %>%
  mutate(p.value = round(p.value, 3))

### cohort specific variance ~ environmental variables ----

# get pop, fishid, age, cohort, env
df_cohort <- data %>% select(pop, fishid, cohort, c.temp, f, ssb.i, recruitment.i) 

# merge ind_plastic and df_cohort
ind_plastic_cohort <- ind_plastic %>%
  full_join(df_cohort, by = join_by(fishid, pop, cohort)) %>%
  group_by(cohort) %>%
  mutate(n_fish = n_distinct(fishid))

# recent cohort - mostly fish with small age -> lower variance
ggplot(data = ind_plastic, aes(x = reaction_norm, y = factor(cohort))) +
  geom_density_ridges(alpha = 0.5) + 
  facet_grid(~ pop) +
  theme_bw()

ggplot(data = ind_plastic, aes(x = cohort, y = reaction_norm)) +
  geom_jitter(alpha = 0.2) + 
  facet_grid(pop ~ .) +
  theme_bw()

# summarize by cohort - only cohort with > 5 fish (Smolinski 2020)
# conduct test for 2 subset of data
# 1. cohorts with all available age
# 2. cohorts with only age > 5 (Morrongiello 2019)

## setup
# list pop
list_pop <- unique(ind_plastic$pop)

# age_range
range(ind_plastic$age)
age_min <- tibble(age_min = c(2, 6),
                  subset = c("all ages", "age 6+"))

df_cor <- tibble()
for(a in 1:2) {
  
  # subset data
  ind_plastic_cohort_sub <- ind_plastic_cohort %>%
    filter(age >= age_min$age_min[a]) %>%
    group_by(cohort) %>%
    mutate(n_fish = n_distinct(fishid)) %>%
    filter(n_fish > 5) 
    
  # calculate var and mean values by cohort
  ind_plastic_cohort_sub <- ind_plastic_cohort_sub %>%
    group_by(pop, cohort) %>%
    mutate(var_cohort = var(reaction_norm)) %>%
    group_by(pop, cohort, var_cohort) %>%
    summarise(across(c.temp:recruitment.i, mean)) %>% 
    ungroup() 
  
  # do correlation test for each population 
  for (i in 1:3) {
    
    df_temp <- ind_plastic_cohort_sub %>% filter(pop == list_pop[i])
    
    df_cor_temp <- cor_test(data = df_temp, 
                            vars = "var_cohort", 
                            vars2 = c("c.temp", "f", "ssb.i", "recruitment.i")) %>%
      mutate(sig_level = if_else(p < 0.05, "significant", "non-significant"),
             pop = list_pop[i],
             subset = age_min$subset[a])
    
    df_cor <- bind_rows(df_cor, df_cor_temp)
  }
}

# save file
write_rds(df_cor, file.path(dir_output, paste0("ind.plastic_cor.test_", unique(data$source), ".rds")))


## oras5 ----
#### data
data <- data_otl %>% 
  left_join(filter(data_temp, source == "oras5"),
            by = join_by(pop, year)) %>%
  filter(is.na(c.temp) == F) %>%
  left_join(data_sol, by = join_by(pop, year))

#### model results
## read model
m <- read_rds(file.path(dir_output, paste0("extrinsic.model.ext_best_", unique(data$source), ".rds")))

## extract thermal reaction norm slope
ind_plastic <- ranef(m)$fishid
ind_plastic <- ind_plastic %>% mutate(fishid = row.names(ind_plastic),
                                      reaction_norm = c.temp_within)
### ratio variance ----
# get pop, fishid, and age
df_fishid <- data %>% select(pop, fishid, age, aac, cohort) %>% 
  group_by(pop, fishid, aac, cohort) %>%
  summarize(age = max(age))

# get ind_plastic for each pop
ind_plastic <- ind_plastic %>% left_join(df_fishid)

## visualisation of reaction norm variation vs age  
ggplot(data = ind_plastic, aes(x = reaction_norm, y = factor(age))) +
  geom_density_ridges(alpha = 0.5) + 
  facet_grid(~ pop) +
  theme_bw()

ggplot(data = ind_plastic, aes(x = age, y = reaction_norm)) +
  geom_jitter(alpha = 0.2) + 
  facet_grid(pop ~ .) +
  theme_bw()

# variation seems to increases with age

# different age group: 1-5, 6-11, 12+ (age >= 5 as suggested in Morrongiello 2019)
# as 8ab data only to age 11
ind_plastic <- ind_plastic %>% mutate(age_group = 
                                        if_else(age <= 5, "Age 2-5",
                                                if_else(age > 5 & age <= 11, "Age 6-11", "Age 12+"))) %>%
  mutate(age_group = factor(age_group, levels = c("Age 2-5", "Age 6-11", "Age 12+")))

ggplot(data = ind_plastic, aes(x = reaction_norm, fill = age_group)) +
  geom_density(alpha = 0.5) + 
  facet_grid(pop ~ .) 

#### variation ratio test
# conduct tests for 4 subsets of data
# age 2 - 11
# age 2 - 22 (all ages)
# age 6 - 11
# age 6 - 22 (all ages)

## setup
# pop_pair - variance ratio test: var_pop1/var_pop2
pop_pair <- tibble(pop1 = c("4bc", "4bc", "7a"),
                   pop2 = c("7a", "8ab", "8ab"))

# age range
range(ind_plastic$age)
age_range <- tibble(min_age = c(2, 2, 6, 6),
                    max_age = c(22, 11, 22, 11))

# note: for each pair, and each age range, 
# sample 100 times the pop with larger sample size to have 2 pops with the same sample size
# for each sample, do variance ratio test with 10000 bootstrap 

n_pop_pair <- nrow(pop_pair)
n_age_range <- nrow(age_range) 
n_sample <- 100
n_boot <- 10000

var_ratio <- tibble()
# set pair
for(p in 1:n_pop_pair) {
  
  print(paste0("processing ", pop_pair$pop1[p], "/", pop_pair$pop2[p]))
  
  pop1 <- pop_pair$pop1[p]
  pop2 <- pop_pair$pop2[p]
  
  df_pop1 <- ind_plastic %>% filter(pop == pop1)
  df_pop2 <- ind_plastic %>% filter(pop == pop2)
  
  # set test 
  for(t in 1:n_age_range) {
    
    print(paste0("processing ", age_range$min_age[t], "-", age_range$max_age[t]))
    
    min_age <- age_range$min_age[t]
    max_age <- age_range$max_age[t]
    
    df_pop1_age <- df_pop1 %>% filter(age >= min_age, age <= max_age)
    df_pop2_age <- df_pop2 %>% filter(age >= min_age, age <= max_age)
    
    # set sample and bootstrap
    results_sample <- tibble()
    for(i in 1:n_sample) {
      set.seed(i)
      
      # subset pop with larger sample size to have 2 pops with the same sample size 
      if (nrow(df_pop1_age) >= nrow(df_pop2_age)) {
        pop1_var <- sample(df_pop1_age$reaction_norm, nrow(df_pop2_age))
        pop2_var <- df_pop2_age$reaction_norm
      } else {
        pop1_var <- df_pop1_age$reaction_norm
        pop2_var <- sample(df_pop2_age$reaction_norm, nrow(df_pop1_age))
      }
      
      df_var <- tibble(pop1_var = pop1_var,
                       pop2_var = pop2_var)
      
      func_var.test <- function(data, i){
        data <- data[i, ]
        results  <- var.test(data$pop1_var, data$pop2_var)
        c(est = results$estimate, CI = results$conf.int, p.value = results$p.value)
      }
      
      boot <- boot(df_var, func_var.test, R = n_boot)
      colnames(boot$t) <- c("var_ratio", "ci_lower", "ci_upper", "p.value")
      var_ratio_temp <- as.data.frame(boot$t) %>% 
        summarize(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>%
        mutate(sample_id = i,
               n_pop1 = nrow(df_pop1_age),
               n_pop2 = nrow(df_pop2_age),
               age_range = paste0("Age ", min_age, " - ", "Age ", max_age),
               pop_pair = paste0(pop1, "/", pop2))
      
      var_ratio <- bind_rows(var_ratio, var_ratio_temp)
    }
  }
}

# save file
write_rds(var_ratio, file.path(dir_output, paste0("ind.plastic_var.ratio_", unique(data$source), ".rds")))

# summary
var_ratio_sum <- var_ratio %>%
  group_by(pop_pair, age_range, n_pop1, n_pop2) %>%
  summarise(across(where(is.numeric), ~ mean(.x))) %>% 
  mutate(across(var_ratio:ci_upper, ~ round(., 2))) %>%
  mutate(p.value = round(p.value, 3))

### cohort specific variance ~ environmental variables ----

# get pop, fishid, age, cohort, env
df_cohort <- data %>% select(pop, fishid, cohort, c.temp, f, ssb.i, recruitment.i) 

# merge ind_plastic and df_cohort
ind_plastic_cohort <- ind_plastic %>%
  full_join(df_cohort, by = join_by(fishid, pop, cohort)) %>%
  group_by(cohort) %>%
  mutate(n_fish = n_distinct(fishid))

# recent cohort - mostly fish with small age -> lower variance
ggplot(data = ind_plastic, aes(x = reaction_norm, y = factor(cohort))) +
  geom_density_ridges(alpha = 0.5) + 
  facet_grid(~ pop) +
  theme_bw()

ggplot(data = ind_plastic, aes(x = cohort, y = reaction_norm)) +
  geom_jitter(alpha = 0.2) + 
  facet_grid(pop ~ .) +
  theme_bw()

# summarize by cohort - only cohort with > 5 fish (Smolinski 2020)
# conduct test for 2 subset of data
# 1. cohorts with all available age
# 2. cohorts with only age > 5 (Morrongiello 2019)

## setup
# list pop
list_pop <- unique(ind_plastic$pop)

# age_min and subset
range(ind_plastic$age)
age_min <- tibble(age_min = c(2, 6),
                  subset = c("all ages", "age 6+"))


df_cor <- tibble()
for(a in 1:2) {
  
  # subset data
  ind_plastic_cohort_sub <- ind_plastic_cohort %>%
    filter(age >= age_min$age_min[a]) %>%
    group_by(cohort) %>%
    mutate(n_fish = n_distinct(fishid)) %>%
    filter(n_fish > 5) 
  
  # calculate var and mean values by cohort
  ind_plastic_cohort_sub <- ind_plastic_cohort_sub %>%
    group_by(pop, cohort) %>%
    mutate(var_cohort = var(reaction_norm)) %>%
    group_by(pop, cohort, var_cohort) %>%
    summarise(across(c.temp:recruitment.i, mean)) %>% 
    ungroup() 
  
  # do correlation test for each population 
  for (i in 1:3) {
    
    df_temp <- ind_plastic_cohort_sub %>% filter(pop == list_pop[i])
    
    df_cor_temp <- cor_test(data = df_temp, 
                            vars = "var_cohort", 
                            vars2 = c("c.temp", "f", "ssb.i", "recruitment.i")) %>%
      mutate(sig_level = if_else(p < 0.05, "significant", "non-significant"),
             pop = list_pop[i],
             subset = age_min$subset[a])
    
    df_cor <- bind_rows(df_cor, df_cor_temp)
  }
}

# save file
write_rds(df_cor, file.path(dir_output, paste0("ind.plastic_cor.test_", unique(data$source), ".rds")))


## nemo-medusa ----
#### data
data <- data_otl %>% 
  left_join(filter(data_temp, source == "nemo-medusa"),
            by = join_by(pop, year)) %>%
  filter(is.na(c.temp) == F) %>%
  left_join(data_sol, by = join_by(pop, year))

#### model results
## read model
m <- read_rds(file.path(dir_output, paste0("extrinsic.model.ext_best_", unique(data$source), ".rds")))

## extract thermal reaction norm slope
ind_plastic <- ranef(m)$fishid
ind_plastic <- ind_plastic %>% mutate(fishid = row.names(ind_plastic),
                                      reaction_norm = c.temp_within)
### ratio variance ----
# get pop, fishid, and age
df_fishid <- data %>% select(pop, fishid, age, aac, cohort) %>% 
  group_by(pop, fishid, aac, cohort) %>%
  summarize(age = max(age))

# get ind_plastic for each pop
ind_plastic <- ind_plastic %>% left_join(df_fishid)

## visualisation of reaction norm variation vs age  
ggplot(data = ind_plastic, aes(x = reaction_norm, y = factor(age))) +
  geom_density_ridges(alpha = 0.5) + 
  facet_grid(~ pop) +
  theme_bw()

ggplot(data = ind_plastic, aes(x = age, y = reaction_norm)) +
  geom_jitter(alpha = 0.2) + 
  facet_grid(pop ~ .) +
  theme_bw()

# variation seems to increases with age

# different age group: 1-5, 6-11, 12+ (age >= 5 as suggested in Morrongiello 2019)
# as 8ab data only to age 11
ind_plastic <- ind_plastic %>% mutate(age_group = 
                                        if_else(age <= 5, "Age 2-5",
                                                if_else(age > 5 & age <= 11, "Age 6-11", "Age 12+"))) %>%
  mutate(age_group = factor(age_group, levels = c("Age 2-5", "Age 6-11", "Age 12+")))

ggplot(data = ind_plastic, aes(x = reaction_norm, fill = age_group)) +
  geom_density(alpha = 0.5) + 
  facet_grid(pop ~ .) 

#### variation ratio test
# conduct tests for 4 subsets of data
# age 2 - 11
# age 2 - 22 (all ages)
# age 6 - 11
# age 6 - 22 (all ages)

## setup
# pop_pair - variance ratio test: var_pop1/var_pop2
pop_pair <- tibble(pop1 = c("4bc", "4bc", "7a"),
                   pop2 = c("7a", "8ab", "8ab"))

# age range
range(ind_plastic$age)
age_range <- tibble(min_age = c(2, 2, 6, 6),
                    max_age = c(22, 11, 22, 11))

# note: for each pair, and each age range, 
# sample 100 times the pop with larger sample size to have 2 pops with the same sample size
# for each sample, do variance ratio test with 10000 bootstrap 

n_pop_pair <- nrow(pop_pair)
n_age_range <- nrow(age_range) 
n_sample <- 100
n_boot <- 10000

var_ratio <- tibble()
# set pair
for(p in 1:n_pop_pair) {
  
  print(paste0("processing ", pop_pair$pop1[p], "/", pop_pair$pop2[p]))
  
  pop1 <- pop_pair$pop1[p]
  pop2 <- pop_pair$pop2[p]
  
  df_pop1 <- ind_plastic %>% filter(pop == pop1)
  df_pop2 <- ind_plastic %>% filter(pop == pop2)
  
  # set test 
  for(t in 1:n_age_range) {
    
    print(paste0("processing ", age_range$min_age[t], "-", age_range$max_age[t]))
    
    min_age <- age_range$min_age[t]
    max_age <- age_range$max_age[t]
    
    df_pop1_age <- df_pop1 %>% filter(age >= min_age, age <= max_age)
    df_pop2_age <- df_pop2 %>% filter(age >= min_age, age <= max_age)
    
    # set sample and bootstrap
    results_sample <- tibble()
    for(i in 1:n_sample) {
      set.seed(i)
      
      # subset pop with larger sample size to have 2 pops with the same sample size 
      if (nrow(df_pop1_age) >= nrow(df_pop2_age)) {
        pop1_var <- sample(df_pop1_age$reaction_norm, nrow(df_pop2_age))
        pop2_var <- df_pop2_age$reaction_norm
      } else {
        pop1_var <- df_pop1_age$reaction_norm
        pop2_var <- sample(df_pop2_age$reaction_norm, nrow(df_pop1_age))
      }
      
      df_var <- tibble(pop1_var = pop1_var,
                       pop2_var = pop2_var)
      
      func_var.test <- function(data, i){
        data <- data[i, ]
        results  <- var.test(data$pop1_var, data$pop2_var)
        c(est = results$estimate, CI = results$conf.int, p.value = results$p.value)
      }
      
      boot <- boot(df_var, func_var.test, R = n_boot)
      colnames(boot$t) <- c("var_ratio", "ci_lower", "ci_upper", "p.value")
      var_ratio_temp <- as.data.frame(boot$t) %>% 
        summarize(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>%
        mutate(sample_id = i,
               n_pop1 = nrow(df_pop1_age),
               n_pop2 = nrow(df_pop2_age),
               age_range = paste0("Age ", min_age, " - ", "Age ", max_age),
               pop_pair = paste0(pop1, "/", pop2))
      
      var_ratio <- bind_rows(var_ratio, var_ratio_temp)
    }
  }
}

# save file
write_rds(var_ratio, file.path(dir_output, paste0("ind.plastic_var.ratio_", unique(data$source), ".rds")))

# summary
var_ratio_sum <- var_ratio %>%
  group_by(pop_pair, age_range, n_pop1, n_pop2) %>%
  summarise(across(where(is.numeric), ~ mean(.x))) %>% 
  mutate(across(var_ratio:ci_upper, ~ round(., 2))) %>%
  mutate(p.value = round(p.value, 3))

### cohort specific variance ~ environmental variables ----

# get pop, fishid, age, cohort, env
df_cohort <- data %>% select(pop, fishid, cohort, c.temp, f, ssb.i, recruitment.i) 

# merge ind_plastic and df_cohort
ind_plastic_cohort <- ind_plastic %>%
  full_join(df_cohort, by = join_by(fishid, pop, cohort)) %>%
  group_by(cohort) %>%
  mutate(n_fish = n_distinct(fishid))

# recent cohort - mostly fish with small age -> lower variance
ggplot(data = ind_plastic, aes(x = reaction_norm, y = factor(cohort))) +
  geom_density_ridges(alpha = 0.5) + 
  facet_grid(~ pop) +
  theme_bw()

ggplot(data = ind_plastic, aes(x = cohort, y = reaction_norm)) +
  geom_jitter(alpha = 0.2) + 
  facet_grid(pop ~ .) +
  theme_bw()

# summarize by cohort - only cohort with > 5 fish (Smolinski 2020)
# conduct test for 2 subset of data
# 1. cohorts with all available age
# 2. cohorts with only age > 5 (Morrongiello 2019)

## setup
# list pop
list_pop <- unique(ind_plastic$pop)

# age_min and subset
range(ind_plastic$age)
age_min <- tibble(age_min = c(2, 6),
                  subset = c("all ages", "age 6+"))


df_cor <- tibble()
for(a in 1:2) {
  
  # subset data
  ind_plastic_cohort_sub <- ind_plastic_cohort %>%
    filter(age >= age_min$age_min[a]) %>%
    group_by(cohort) %>%
    mutate(n_fish = n_distinct(fishid)) %>%
    filter(n_fish > 5) 
  
  # calculate var and mean values by cohort
  ind_plastic_cohort_sub <- ind_plastic_cohort_sub %>%
    group_by(pop, cohort) %>%
    mutate(var_cohort = var(reaction_norm)) %>%
    group_by(pop, cohort, var_cohort) %>%
    summarise(across(c.temp:recruitment.i, mean)) %>% 
    ungroup() 
  
  # do correlation test for each population 
  for (i in 1:3) {
    
    df_temp <- ind_plastic_cohort_sub %>% filter(pop == list_pop[i])
    
    df_cor_temp <- cor_test(data = df_temp, 
                            vars = "var_cohort", 
                            vars2 = c("c.temp", "f", "ssb.i", "recruitment.i")) %>%
      mutate(sig_level = if_else(p < 0.05, "significant", "non-significant"),
             pop = list_pop[i],
             subset = age_min$subset[a])
    
    df_cor <- bind_rows(df_cor, df_cor_temp)
  }
}

# save file
write_rds(df_cor, file.path(dir_output, paste0("ind.plastic_cor.test_", unique(data$source), ".rds")))


