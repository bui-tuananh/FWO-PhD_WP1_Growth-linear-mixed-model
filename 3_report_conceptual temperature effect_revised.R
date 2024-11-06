library(tidyverse)
library(lme4)
library(effects)

# create a customised theme
theme_custom <- theme_classic() +
  theme(legend.background = element_rect(fill='transparent'),
        legend.key = element_rect(fill='transparent'),
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        legend.title = element_blank()) +
  theme(legend.position = c(0, 1),
        legend.justification = c(0.1, 0.75),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0.2),
        legend.text = element_text(size = 4, margin = margin(l = 1, unit = "pt")),
        legend.key.size = unit(6, "pt"),
        plot.margin = unit(c(1,1,1,1), "pt")) 

## temperature panel ----
### temp data (population) -----
#### setup
pop <- c("Population 1", "Population 2", "Population 3")
n_pop <- length(pop)
n_year <- 10

#### data
set.seed(17) #17 for 10 year, 24 for 20 year
data_temp <- tibble(year = rep(seq(1, n_year), n_pop),
                    temp = rep(rnorm(n_year, mean = 6, sd = 0.5), n_pop),
                    pop = rep(pop, each = n_year)) %>%
  mutate(temp = temp * rep(c(1, 1.2, 1.5), each = n_year)) 

#pop_ave vs pop_ano
data_temp <- data_temp %>% 
  group_by(pop) %>% 
  mutate(pop_ave = mean(temp),
         pop_ano = temp - pop_ave) 

# manual arrangement to have increasing trend in year 2-4 and 7-9
data_temp <- data_temp %>% mutate(year = c(1,4,3,2,7,6,5,9,8,10))

### fish data (individual) ----
#### setup
# each pop has 2 cohorts, each cohort has 3 fish, each fish has 3 ages 
cohort <- c(2,7) # cohort 2: low temp, cohort 7: high temp
n_cohort = length(cohort)
n_age = 3
n_fish = 3*n_pop*n_cohort # each pop and each cohort 3 fish

df_period <- tibble(cohort = cohort,
                    period_mean = c(3,8),
                    period_name = c("Cool period",
                                    "Warm period"))

data_fish <- tibble(
  fish = 1:n_fish,
  pop = rep(pop, each = n_fish/n_pop),
  cohort = rep(cohort, n_fish/n_cohort)) %>%
  crossing(age = seq(1, n_age)) %>%
  mutate(year = cohort + age - 1)

# calculate ind_temp
data_fish <- data_fish %>% 
  left_join(data_temp) %>%
  group_by(fish) %>% 
  mutate(ind_ave = mean(pop_ano),
         ind_ano = pop_ano - ind_ave) %>%
  ungroup() %>%
  left_join(df_period)

### plot ----
#### temp
ggplot(data = data_temp, aes(x = year, y = temp, color = pop)) +
  geom_point(size = 0.1) +
  geom_line(linewidth = 0.2) +
  #geom_point(alpha = 0) +
  labs(x = "Year",
       y = "Temperature (°C)",
       color = "Population") +
  scale_color_brewer(palette = "Dark2", direction = -1) +
  theme_custom

p_temp <- last_plot()

#### pop-average pop-anomaly 
# pop-average
ggplot(data = data_temp, aes(x = pop, y = pop_ave, color = pop)) + 
  geom_point(size = 0.1) +
  labs(x = "Population",
       y = expression('T'['population-average'] * ' (°C)'),
       color = "Population") +
  scale_color_brewer(palette = "Dark2", direction = -1) +
  theme_custom

p_temp_pop.ave <- last_plot()

#### pop-anomaly
ggplot(data = data_temp, aes(x = year, y = pop_ano, color = pop)) + 
  annotate("rect", xmin = 2, xmax = 4, 
           ymin = -Inf, 
           ymax = Inf,
           fill = "#313695",
           alpha = 0.2) +
  annotate("rect", xmin = 7, xmax = 9, 
           ymin = -Inf, 
           ymax = Inf,
           fill = "#a50026",
           alpha = 0.2) +
  geom_point(size = 0.1) +
  geom_line(linewidth = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.2) +
  labs(x = "Year",
       y = expression('T'['population-anomaly'] * ' (°C)'),
       color = "Population") +
  scale_color_brewer(palette = "Dark2", direction = -1) +
  annotate("text", 
           x = 3,
           y = min(data_temp$pop_ano), 
           label = "Cool period",
           size = 4/.pt) +
  annotate("text", 
           x = 8,
           y = min(data_temp$pop_ano), 
           label = "Warm period",
           size = 4/.pt) +
  theme_custom

p_temp_pop.ano <- last_plot()

#### ind-average
ggplot(data = data_fish, aes(x = period_mean, y = ind_ave, color = period_name)) + 
  geom_point(size = 0.1) +
  labs(x = "Period",
       y = expression('T'['individual-average'] * ' (°C)'),
       color = "Population") +
  scale_color_manual(values = c("Cool period" = "#313695", "Warm period" = "#a50026")) +
  xlim(1,10) +
  theme_custom

p_temp_ind.ave <- last_plot()

#### ind-anomaly
ggplot(data = data_fish, aes(x = age, y = ind_ano, group = fish, color = period_name)) + 
  geom_point(size = 0.1) +
  geom_line(linewidth = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.2) +
  labs(x = "Year",
       y = expression('T'['individual-anomaly'] * ' (°C)'),
       color = "Population") +
  scale_color_manual(values = c("Cool period" = "#313695", "Warm period" = "#a50026")) +
  theme_custom

p_temp_ind.ano <- last_plot()

## growth panel ----
### data ----
set.seed(1)
data <- data_fish %>%
  mutate(intercept_fishid = rep(rnorm(n_fish, 0, 0.15), each = n_age),
         slope_fishid = rep(rnorm(n_fish, 0, 0.1), each = n_age)) %>%
  mutate(
    log_age = log(age),
    
    # growth - temperature effect
    growth_fish = (intercept_fishid) + 0.3*pop_ave + 0.3*ind_ave + (0.3 + slope_fishid)*ind_ano,
    growth = growth_fish + rnorm(nrow(.), 0, 0.01),
    
    # growth - temperature + age effect
    log_growth_at_age_true = 3 - 1.2*log_age, 
    log_growth_at_age = log_growth_at_age_true + growth,
    growth_at_age = exp(log_growth_at_age),
    
    # growth-temperature true
    growth_pop_ave_true = exp(3 - 1.2*mean(log_age) + 0.3*pop_ave + 0.3*mean(pop_ano)),
    growth_pop_ano_true = exp(3 - 1.2*mean(log_age) + 0.3*mean(pop_ave) + 0.3*pop_ano),
    growth_ind_ave_true = exp(3 - 1.2*mean(log_age) + 0.3*mean(pop_ave) + 0.3*ind_ave + 0.3*mean(ind_ano)),
    growth_ind_ano_true = exp(3 - 1.2*mean(log_age) + 0.3*mean(pop_ave) + 0.3*mean(ind_ave) + 0.3*ind_ano),
    growth_ind_ano_fish = exp((3+intercept_fishid) - 1.2*mean(log_age) + 0.3*mean(pop_ave) + 0.3*mean(ind_ave) + (0.3 + slope_fishid)*ind_ano)
  )

### plot ----
#### growth 
ggplot(data = data, aes(x = age, y = growth_at_age, group = fish, color = pop)) +
  geom_point(size = 0.1) +
  geom_line(linewidth = 0.2) +
  labs(x = "Age",
       y = "Growth",
       color = "Population") +
  scale_color_brewer(palette = "Dark2", direction = -1) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(legend.background = element_rect(fill='transparent'),
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        legend.title = element_blank()) +
  theme_custom +
  theme(legend.position = c(0.6, 1))

p_growth <- last_plot()

#### pop-average
ggplot(data = data, aes(x = pop, y = growth, group = fish, color = pop)) +
  geom_point(size = 0.1, alpha = 0.5) +
  labs(x = "Population",
       y = "Growth",
       color = "Population") +
  scale_color_brewer(palette = "Dark2", direction = -1) +
  theme_custom

p_growth_pop.ave <- last_plot()

#### pop-anomaly
ggplot(data = data, aes(x = year, y = growth, group = fish, color = pop)) + 
  annotate("rect", xmin = 2, xmax = 4, 
           ymin = -Inf, 
           ymax = Inf,
           fill = "#313695",
           alpha = 0.2) +
  annotate("rect", xmin = 7, xmax = 9, 
           ymin = -Inf, 
           ymax = Inf,
           fill = "#a50026",
           alpha = 0.2) +
  geom_point(size = 0.1) +
  geom_line(linewidth = 0.2) +
  labs(x = "Year",
       y = "Growth",
       color = "Population") +
  scale_color_brewer(palette = "Dark2", direction = -1) +
  annotate("text", 
           x = 3,
           y = min(data$growth), 
           label = "Cool period",
           size = 4/.pt) +
  annotate("text", 
           x = 8,
           y = min(data$growth), 
           label = "Warm period",
           size = 4/.pt) +
  xlim(1,10) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme_custom

p_growth_pop.ano <- last_plot()

#### ind-average
ggplot(data = data, aes(x = period_mean, y = growth, group = fish, color = period_name)) +
  geom_point(size = 0.1, alpha = 0.5) +
  labs(x = "Period",
       y = "Growth") +
  scale_color_manual(values = c("Cool period" = "#313695", "Warm period" = "#a50026")) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  xlim(1,10) +
  theme_custom

p_growth_ind.ave <- last_plot()

# ind-anomaly
ggplot(data = data, aes(x = age, y = growth, group = fish, color = period_name)) +
  geom_point(size = 0.1) +
  geom_line(linewidth = 0.2) +
  labs(x = "Year",
       y = "Growth") +
  scale_color_manual(values = c("Cool period" = "#313695", "Warm period" = "#a50026")) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme_custom

p_growth_ind.ano <- last_plot()

## growth-temperature response panel ----
### note ----
# note: this panel can be visualised using either
# (1) - model results (pop and ind levels)
# (2) - simulation data (columns in section growth-temperature true of growth data)

# model
m_pop <- lmer(log_growth_at_age ~ 1 + log_age + pop_ano + pop_ave + (1 | fish), data = data)
m_ind <- lmer(log_growth_at_age ~ 1 + log_age + pop_ave + ind_ave + ind_ano + (1 + ind_ano | fish), data = data)

### plot ----
#### pop-average
# simulation data
ggplot(data = data, aes(x = pop_ave, y = growth_pop_ave_true)) +
  geom_line()

# use model results for smooth effect (if use simulation data, there are only 3 points)
pred <- as.data.frame(Effect(c("pop_ave"), m_pop, xlevels = 10))
ggplot(data = pred, aes(x = pop_ave, y = exp(fit))) +
  geom_line(linewidth = 0.3) +
  labs(x = expression('T'['population-average'] * ' (°C)'),
       y = "Growth") +
  theme_custom +
  ylim(59, 155) #from growth_pop_ave_true

p_response_pop.ave <- last_plot()

#### pop-anomaly
ggplot(data = data, aes(x = pop_ano, y = growth_pop_ano_true)) +
  geom_line(linewidth = 0.4) +
  labs(x = expression('T'['population-anomaly'] * ' (°C)'),
       y = "Growth") +
  theme_custom +
  ylim(59, 155)

p_response_pop.ano <- last_plot()

#### ind-average
ggplot(data = data, aes(x = ind_ave, y = growth_ind_ave_true)) +
  geom_line(linewidth = 0.4) +
  labs(x = expression('T'['individual-average'] * ' (°C)'),
       y = "Growth") +
  theme_custom +
  ylim(59, 155)# +
  #xlim(-0.5, 0.61) #ind_ano

p_response_ind.ave <- last_plot()

#### ind-anomaly
ggplot() +
  geom_line(data = data, 
            aes(x = ind_ano, y = growth_ind_ano_true),
            linewidth = 0.4) +
  geom_line(data = data, 
            aes(x = ind_ano, y = growth_ind_ano_fish, 
                group = fish), 
            alpha = 0.3,
            linewidth = 0.15) +
  labs(x = expression('T'['individual-anomaly'] * ' (°C)'),
       y = "Growth") +
  theme_custom +
  ylim(59, 155)# +
  #xlim(-0.5, 0.61)

p_response_ind.ano <- last_plot()

## save plot ----
plot_list <- list("p_temp" = p_temp,
                  "p_temp_pop.ave" = p_temp_pop.ave,
                  "p_temp_pop.ano" = p_temp_pop.ano,
                  "p_temp_ind.ave" = p_temp_ind.ave,
                  "p_temp_ind.ano" = p_temp_ind.ano,
                  "p_growth" = p_growth,
                  "p_growth_pop.ave" = p_growth_pop.ave,
                  "p_growth_pop.ano" = p_growth_pop.ano,
                  "p_growth_ind.ave" = p_growth_ind.ave,
                  "p_growth_ind.ano" = p_growth_ind.ano,
                  "p_response_pop.ave" = p_response_pop.ave,
                  "p_response_pop.ano" = p_response_pop.ano,
                  "p_response_ind.ave" = p_response_ind.ave,
                  "p_response_ind.ano" = p_response_ind.ano)

dir_report <- "./report"
for (i in 1:length(plot_list)) {
  ggsave(plot_list[[i]], file = file.path(dir_report, "fig3_conceptual temperature effect", paste0(names(plot_list[i]),".pdf") ),
         device = cairo_pdf,
         width = 3.3, height = 3.3,
         units = "cm")
  
  # ggsave(plot_list[[i]], file = file.path(dir_report, "fig3_conceptual temperature effect", paste0(names(plot_list[i]),".tiff") ),
  #        width = 3.3, height = 3.3,
  #        units = "cm",
  #        dpi = 1000,
  #        device = "tiff")
}
