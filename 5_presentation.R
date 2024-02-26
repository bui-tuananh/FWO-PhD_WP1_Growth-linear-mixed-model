
#BEFORE
#   3_report files

#AFTER
#   WP1/presentation


# 1. SETUP ----

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
library(patchwork)  # arrange ggplot 

# FUNCTION
#centring function
c. <- function (x) {(x - mean(x))} 

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
dir_output <- "./WP1/output"
dir_report <- "./WP1/report"
dir_presentation <- "./WP1/presentation"

# COLOR PALETTE
pal <- (brewer.pal(8, "Spectral")) 


# 2. LOAD DATA ----
## 2.1. Otolith data ----

### 2.1.1 Process otl data ---------------------------------------------------

dir_otl <- "D:/OneDrive - UGent/data/WP1/otolith"
otl <- read_rds(file.path(dir_otl, "otl_full.rds"))

# Workflow
# 1. remove 1 fish 8ab cohort < 1985
# 2. add pop.name
otl <- otl %>% 
  filter(FishID != "SOL_52_G1_Z.121_26-03-1990_1602") %>%
  mutate(pop.name = if_else(IcesAreaGroup == "4bc", "North Sea",
                            if_else(IcesAreaGroup == "7a", "Irish Sea", 
                                    "Bay of Biscay")))

# setup
## check 0 increment
incremet_0 <- otl %>% filter(AnnulusDiameterIncrement.um == 0) 
## check number of increment/year
data_sum <- otl %>% group_by(GrowingYear) %>% summarize(n_increment = n())
## create ssb index by dividing ssb by ices area
dir_gis <- "D:/OneDrive - UGent/data/Admin"
ices_4bc <- as.data.frame(st_read(file.path(dir_gis, "ices_areas_sub_group_4326_new.gpkg")))
ices_4bc <- ices_4bc %>% 
  filter(Area_27 %in% c("4bc", "7a", "8ab")) %>% 
  rename(IcesAreaGroup.area.km2 = Area_km2,
         IcesAreaGroup = Area_27) %>%
  select(IcesAreaGroup, IcesAreaGroup.area.km2) 

# Workflow
# 1. remove 1 wur fish with very small age 1 increment (< 200nm)
# 2. 0 increment, add a small value 0.05 to avoid log issue (smallest increment from ILVO is 0.06)
# 3. remove years with less than 30 increments (van der Sleen et al., 2018)
# 4. add ssb index (ssb/ices area)

data <- otl %>% 
  filter(FishID != "sol_fab_0575") %>%
  mutate(AnnulusDiameterIncrement.um = if_else(AnnulusDiameterIncrement.um == 0, 0.05, AnnulusDiameterIncrement.um)) %>%
  filter(GrowingYear > 1957) %>%
  left_join(ices_4bc, by = "IcesAreaGroup")

### 2.1.1 Standardize explanatory variables ----------------------------------------

# standardize and change variables name to formulate models easier
## Age, AgeAtCapture: log-transformed, then centered: x - mean(x)
## sbt, f, ssb.i (ssb/area_sqkm): centered: x - mean(x) 

data <- data %>% mutate(fishid = FishID,
                        increment = AnnulusDiameterIncrement.um,
                        age = Age,
                        c.log.age = c.(log(Age)),
                        aac = AgeAtCapture,
                        c.log.aac = c.(log(AgeAtCapture)), 
                        method = OtolithProcessingMethod,
                        datasource = DataSource,
                        pop = IcesAreaGroup,
                        pop.area.km2 = IcesAreaGroup.area.km2,
                        year = GrowingYear,
                        cohort = Cohort,
                        fyear = factor(GrowingYear),
                        fcohort = factor(Cohort),
                        pop.year = paste0(pop, ":", year),
                        pop.cohort = paste0(pop, ":", cohort),
                        # predictors
                        sbt = SeaBottomTemperature.degC,
                        c.sbt = c.(SeaBottomTemperature.degC),
                        f = FishingMortality,
                        c.f = c.(FishingMortality),
                        ssb.i = SpawningStockBiomass.1000t*1000/IcesAreaGroup.area.km2, #unit: ton/km2
                        c.ssb.i = c.(ssb.i))

data <- data %>% mutate(pop.name = if_else(IcesAreaGroup == "4bc", "North Sea",
                                           if_else(IcesAreaGroup == "7a", "Irish Sea", 
                                                   "Bay of Biscay")))

### 2.1.3 Add climwin outputs -------------------
data <- data %>%
  mutate(sbt.among = ave(SeaBottomTemperature.degC, FishID),
         sbt.within = sbt - sbt.among) 

## 2.4. Get population name
df_pop <- tibble(pop = c("4bc", "7a", "8ab"),
                 pop.name = c("North Sea", "Irish Sea", "Bay of Biscay"))

## 2.2. Model data ----

#BEST INTRINSIC STRUCTURE: c.log.age + c.log.aac*pop + method + datasource*c.log.age
m3 <- read_rds(file.path(dir_output, "intrinsic.model_best.rds"))

# BEST EXTRINSIC STRUCTURE: c.ssb.i*c.log.age
m4 <- read_rds(file.path(dir_output, "extrinsic.model_best.rds"))

# BEST EXTENDED EXTRINSIC STRUCTURE: c.ssb.i*c.log.age + sbt.within*c.log.age + sbt.among + sbt.within*c.ssb.i + sbt.within*sbt.among
m5 <- read_rds(file.path(dir_output, "extrinsic.model.ext_best.rds"))




# 3. FIGURE
## sbt, sbt.within, sbt.among ----
### setup ----
## optimal extrinsic structure + sbt effect
m_sbt <- lmer(log(increment) ~ 1 + c.log.age + c.log.aac*pop + method + datasource*c.log.age + 
                c.ssb.i + sbt + 
                (1 | fishid) + (1 + c.log.age | pop.year), 
              data = data, REML = T)

# sbt
df_sbt <- unique(select(data, sbt)) %>% round(1) %>% unique()
pred_sbt <- as.data.frame(Effect(c("sbt"), 
                                 m_sbt, 
                                 xlevels = list(sbt = df_sbt$sbt))) 

# sbt.among
pred_sbt.among <- as.data.frame(Effect(c("sbt.among"), 
                                       m5, 
                                       xlevels = list(sbt.among = df_sbt$sbt))) 

# sbt.within
pred_sbt.within <- as.data.frame(Effect(c("sbt.within", "sbt.among"), 
                                        m5, 
                                        xlevels = list(sbt.within = seq(-1, 1, 0.1),
                                                       sbt.among = c(8.2, 10, 11.7)))) %>%
  mutate(sbt = sbt.among + sbt.within) 

### plot ----

#### sbt ----
ggplot() + 
  # sbt
  geom_line(data = pred_sbt, 
            aes(x = sbt, 
                y = exp(fit),
                color = "population average",
                linewidth = "population average",
                linetype = "population average")) +
  geom_ribbon(data = pred_sbt, 
              aes(x = sbt, 
                  ymin = exp(lower),
                  ymax = exp(upper)),
              fill = "grey",
              alpha = 0.3) +
  theme_classic() + 
  labs(x = "Temperature (°C)",
       y = "Predicted increment growth (μm)",
       color = "Temperature effect",
       linewidth = "Temperature effect",
       linetype = "Temperature effect") +
  scale_colour_manual(breaks = c("population average","within-individual","between-individual"),
                      values = c("grey", "black","black")) +
  scale_linewidth_manual(breaks = c("population average","within-individual","between-individual"),
                         values = c(1, 0.3, 1)) +
  scale_linetype_manual(breaks = c("population average","within-individual","between-individual"),
                        values = c("twodash", "solid", "solid")) +
  scale_y_continuous(limits = c(238, 500)) +
  scale_x_continuous(breaks = c(8,9,10,11,12))

ggsave(last_plot(), file = file.path(dir_presentation, "sbt_sbt.png"),
       width = 6, height = 4, 
       dpi = 600)

#### sbt.within ----

ggplot() + 
  # sbt
  geom_line(data = pred_sbt, 
            aes(x = sbt, 
                y = exp(fit),
                color = "population average",
                linewidth = "population average",
                linetype = "population average")) +
  geom_ribbon(data = pred_sbt, 
              aes(x = sbt, 
                  ymin = exp(lower),
                  ymax = exp(upper)),
              fill = "grey",
              alpha = 0.3) +
  # sbt.within
  geom_line(data = pred_sbt.within, 
            aes(x = sbt, 
                y = exp(fit), 
                group = sbt.among,
                color = "within-individual",
                linewidth = "within-individual",
                linetype = "within-individual")) +
  theme_classic() + 
  labs(x = "Temperature (°C)",
       y = "Predicted increment growth (μm)",
       color = "Temperature effect",
       linewidth = "Temperature effect",
       linetype = "Temperature effect") +
  scale_colour_manual(breaks = c("population average","within-individual","between-individual"),
                      values = c("grey", "black","black")) +
  scale_linewidth_manual(breaks = c("population average","within-individual","between-individual"),
                         values = c(1, 0.3, 1)) +
  scale_linetype_manual(breaks = c("population average","within-individual","between-individual"),
                        values = c("twodash", "solid", "solid")) +
  scale_y_continuous(limits = c(238, 500)) +
  scale_x_continuous(breaks = c(8,9,10,11,12))

ggsave(last_plot(), file = file.path(dir_presentation, "sbt_sbt.within.png"),
       width = 6, height = 4, 
       dpi = 600)

#### sbt.among ----
ggplot() + 
  # sbt
  geom_line(data = pred_sbt, 
            aes(x = sbt, 
                y = exp(fit),
                color = "population average",
                linewidth = "population average",
                linetype = "population average")) +
  geom_ribbon(data = pred_sbt, 
              aes(x = sbt, 
                  ymin = exp(lower),
                  ymax = exp(upper)),
              fill = "grey",
              alpha = 0.3) +
  # sbt.among
  geom_line(data = pred_sbt.among, 
            aes(x = sbt.among, 
                y = exp(fit),
                color = "between-individual",
                linewidth = "between-individual",
                linetype = "between-individual")) +
  # sbt.within
  geom_line(data = pred_sbt.within, 
            aes(x = sbt, 
                y = exp(fit), 
                group = sbt.among,
                color = "within-individual",
                linewidth = "within-individual",
                linetype = "within-individual")) +
  theme_classic() + 
  labs(x = "Temperature (°C)",
       y = "Predicted increment growth (μm)",
       color = "Temperature effect",
       linewidth = "Temperature effect",
       linetype = "Temperature effect") +
  scale_colour_manual(breaks = c("population average","within-individual","between-individual"),
                      values = c("grey", "black","black")) +
  scale_linewidth_manual(breaks = c("population average","within-individual","between-individual"),
                         values = c(1, 0.3, 1)) +
  scale_linetype_manual(breaks = c("population average","within-individual","between-individual"),
                        values = c("twodash", "solid", "solid")) +
  scale_y_continuous(limits = c(238, 500)) +
  scale_x_continuous(breaks = c(8,9,10,11,12))

  
ggsave(last_plot(), file = file.path(dir_presentation, "sbt_sbt.among.png"),
       width = 6, height = 4, 
       dpi = 600)
  

## sbt.within*age - width ----
### setup ----
# setup
df_age <- unique(select(data, age, c.log.age)) %>% round(2)
df_sbt.within <- unique(select(data, sbt.within)) %>% round(1) %>% unique()

df_legend <- tibble(sbt.within = c(-1, 0, 1),
                    legend = factor(c("-1°C cooler", "Lifetime average", "+1°C warmer"),
                                    levels = c("-1°C cooler", "Lifetime average", "+1°C warmer")))

pred <- as.data.frame(Effect(c("c.log.age", "sbt.within"), 
                             m5, 
                             xlevels = list(c.log.age = df_age$c.log.age,
                                            sbt.within = c(-1, 0, 1)))) %>%
  left_join(df_age) %>%
  left_join(df_legend)

## calculate otolith width 
pred <- pred %>% 
  group_by(sbt.within) %>%
  arrange(age) %>% 
  mutate(width = cumsum(exp(fit)),
         width_lower = cumsum(exp(lower)),
         width_upper = cumsum(exp(upper)))

### plot ----
#### warmer ----
ggplot() +
  geom_line(data = pred %>% filter(legend != "-1°C cooler"), 
            aes(x = age, 
                y = width, 
                color = legend),
            linewidth = 1) +
  geom_ribbon(data = pred %>% filter(legend != "-1°C cooler"), 
              aes(x = age, 
                  ymin = width_lower, 
                  ymax = width_upper,
                  fill = legend),
              alpha = 0.3) +
  scale_colour_manual(values = c("grey","#F8766D")) +
  scale_fill_manual(values = c("grey","#F8766D")) +
  theme_classic() +
  labs(x = "Age (years)",
       y = "Predicted width (μm)",
       color = expression('Temperature'['within-individual']),
       fill = expression('Temperature'['within-individual'])) +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(1490, 7100)) 

#### cooler ----
ggplot() +
  geom_line(data = pred, 
            aes(x = age, 
                y = width, 
                color = legend),
            linewidth = 1) +
  geom_ribbon(data = pred, 
              aes(x = age, 
                  ymin = width_lower, 
                  ymax = width_upper,
                  fill = legend),
              alpha = 0.3) +
  scale_colour_manual(values = c("#00BFC4","grey","#F8766D")) +
  scale_fill_manual(values = c("#00BFC4","grey","#F8766D")) +
  theme_classic() +
  labs(x = "Age (years)",
       y = "Predicted width (μm)",
       color = expression('Temperature'['within-individual']),
       fill = expression('Temperature'['within-individual'])) +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(1490, 7100)) 


