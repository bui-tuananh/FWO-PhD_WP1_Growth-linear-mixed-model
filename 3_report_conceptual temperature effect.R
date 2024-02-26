library(tidyverse)

theme_set(theme_classic()) 


## temp data ----
#### setup
pop <- c("Population 1", "Population 2", "Population 3")
n_pop <- length(pop)
n_year <- 10

#### data
set.seed(17) #17 for 10 year, 24 for 20 years
data_temp <- tibble(year = rep(seq(1, n_year), n_pop),
                    temp = rep(rnorm(n_year, mean = 6, sd = 0.5), n_pop),
                    pop = rep(pop, each = n_year)) %>%
  mutate(temp = temp * rep(c(1, 1.2, 1.5), each = n_year)) 

#### plot
ggplot(data = data_temp, aes(x = year, y = temp, color = pop)) +
  geom_line() +
  #geom_point(alpha = 0) +
  labs(x = "Year",
       y = "Temperature (°C)",
       color = "Population") +
  scale_color_brewer(palette = "Dark2", direction = -1) +
  theme(legend.position = c(0, 1.05),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill='transparent'),
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        legend.title = element_blank()) #+
  #guides(colour = guide_legend(override.aes = list(alpha = 1)))

p_temp <- last_plot()

#### pop-average pop-anomaly 
#### data
data_temp <- data_temp %>% 
  group_by(pop) %>% 
  mutate(pop_ave = mean(temp),
         pop_ano = temp - pop_ave)

#### plot
# pop-average
ggplot(data = data_temp, aes(x = pop, y = pop_ave, color = pop)) + 
  geom_point() +
    labs(x = "Population",
         y = expression('Temperature'['population-average'] * ' (°C)'),
         color = "Population") +
  scale_color_brewer(palette = "Dark2", direction = -1) +
  theme(legend.position = c(0, 1.05),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill='transparent'),
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        legend.title = element_blank()) #+
  #theme(legend.position = "none")

p_temp_pop.ave <- last_plot()

# pop-anomaly
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
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Year",
       y = expression('Temperature'['population-anomaly'] * ' (°C)'),
       color = "Population") +
  scale_color_brewer(palette = "Dark2", direction = -1) +
  theme(legend.position = c(0, 1.05),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill='transparent'),
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        legend.title = element_blank()) +
  #theme(legend.position = "none") +
  annotate("text", 
           x = 3,
           y = min(data_temp$pop_ano), 
           label = "Cool period",
           size = 2.5) +
  annotate("text", 
           x = 8,
           y = min(data_temp$pop_ano), 
           label = "Warm period",
           size = 2.5) 

p_temp_pop.ano <- last_plot()

## fish data ----
#### setup
# each pop has 2 cohorts, each cohort has 3 fish, each fish has 3 ages 
cohort <- c(2,7) # cohort 2: low temp, cohort 7: high temp
n_cohort = length(cohort)
n_age = 3
n_fish = 3*n_pop*n_cohort # each pop and each cohort 3 fish

data_fish <- tibble(
  fish = 1:n_fish,
  pop = rep(pop, each = n_fish/n_pop),
  cohort = rep(cohort, n_fish/n_cohort)) %>%
  crossing(age = seq(1, n_age)) %>%
  mutate(year = cohort + age - 1)

### population-level growth response to temperature ----
#### data
data_pop <- data_fish %>% 
  left_join(data_temp)

set.seed(1)
data_pop <- data_pop %>%
  mutate(growth_pop_ave_true = 1 + 0.5*pop_ave,
         growth_pop_ave = 1 + 0.5*pop_ave + rnorm(nrow(.), 0, 0.1),
         growth_pop_ano_true = 1 + 0.5*pop_ano,
         growth_pop_ano = 1 + 0.5*pop_ano + rnorm(nrow(.), 0, 0.1)
         )

#### plot
# pop-average
ggplot(data = data_pop) +
  geom_point(aes(x = pop_ave, y = growth_pop_ave, color = pop), alpha = 0.5) +
  geom_line(aes(x = pop_ave, y = growth_pop_ave_true), linewidth = 1) +
  labs(x = expression('Temperature'['population-average'] * ' (°C)'),
       y = "Growth",
       color = "Population") +
  scale_color_brewer(palette = "Dark2", direction = -1) +
  theme(legend.position = c(0, 1.05),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill='transparent'),
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        legend.title = element_blank()) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) #+
  #theme(legend.position = "none")

p_growth_pop.ave <- last_plot()

# pop-anomaly
ggplot(data = data_pop) +
  #geom_point(aes(x = pop_ano, y = growth_pop_ano, color = pop), alpha = 0.5) +
  geom_point(aes(x = pop_ano, y = growth_pop_ano), alpha = 0.5) +
  geom_line(aes(x = pop_ano, y = growth_pop_ano_true), linewidth = 1) +
  labs(x = expression('Temperature'['population-anomaly'] * ' (°C)'),
       y = "Growth",
       color = "Population") +
  scale_color_brewer(palette = "Dark2", direction = -1) +
  theme(legend.position = c(0, 1.05),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill='transparent'),
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        legend.title = element_blank()) +
  guides(colour = guide_legend(override.aes = list(alpha = 1)))

p_growth_pop.ano <- last_plot()

### individual-level growth response to temperature: plasticity and between-individual ----
#### data
df_period <- tibble(cohort = cohort,
                    period_mean = c(3,8),
                    period_name = c("Cool period",
                                    "Warm period"))

data_ind <- data_fish %>% 
  left_join(data_temp) %>%
  group_by(fish) %>% 
  mutate(ind_ave = mean(pop_ano),
         ind_ano = pop_ano - ind_ave) %>%
  ungroup() %>%
  left_join(df_period)

set.seed(2) #
data_ind <- data_ind %>%
  mutate(intercept_fishid = rep(rnorm(n_fish, 0, 0.15), each = n_age),
         slope_fishid = rep(rnorm(n_fish, 0, 0.1), each = n_age)) %>%
  mutate(growth_ind_ave_true = 1 + 0.5*ind_ave,
         growth_ind_ave = 1 + 0.5*ind_ave + rnorm(nrow(.), 0, 0.1),
         growth_ind_ano_true = 1 + 0.5*ind_ano,
         growth_ind_ano_fish = (1 + intercept_fishid) + (0.5 + slope_fishid)*ind_ano,
         growth_ind_ano = growth_ind_ano_fish + rnorm(nrow(.), 0, 0.01)
  )

#### temperature plot ----
# ind-average
ggplot(data = data_ind, aes(x = period_mean, y = ind_ave, color = period_name)) + 
  geom_point() +
  #geom_line(alpha = 0) +
  labs(x = "Period",
       y = expression('Temperature'['individual-average'] * ' (°C)'),
       color = "Population") +
  scale_color_manual(values = c("Cool period" = "#313695", "Warm period" = "#a50026")) +
  theme(legend.position = c(0, 1.05),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill='transparent'),
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        legend.title = element_blank()) +
  #guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  xlim(1,10)

p_temp_ind.ave <- last_plot()

# ind-anomaly
ggplot(data = data_ind, aes(x = age, y = ind_ano, group = fish, color = period_name)) + 
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Year",
       y = expression('Temperature'['individual-anomaly'] * ' (°C)'),
       color = "Population") +
  scale_color_manual(values = c("Cool period" = "#313695", "Warm period" = "#a50026")) +
  theme(legend.position = c(0, 1.05),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill='transparent'),
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        legend.title = element_blank()) #+
  #theme(legend.position = "none")

p_temp_ind.ano <- last_plot()

#### growth response plot ----
# ind-average
ggplot(data = data_ind) +
  geom_point(aes(x = ind_ave, y = growth_ind_ave, color = period_name), alpha = 0.5) +
  geom_line(aes(x = ind_ave, y = growth_ind_ave_true), linewidth = 1) +
  labs(x = expression('Temperature'['individual-average'] * ' (°C)'),
       y = "Growth",
       color = "Population") +
  scale_color_manual(values = c("Cool period" = "#313695", "Warm period" = "#a50026")) +
  theme(legend.position = c(0, 1.05),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill='transparent'),
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        legend.title = element_blank()) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) #+
  #theme(legend.position = "none")

p_growth_ind.ave <- last_plot()

# ind-anomaly
ggplot(data = data_ind) +
  #geom_point(aes(x = ind_ano, y = growth_ind_ano, color = period_name), alpha = 0.5) +
  geom_point(aes(x = ind_ano, y = growth_ind_ano, group = fish), alpha = 0.5) +
  geom_line(aes(x = ind_ano, y = growth_ind_ano_fish, group = factor(fish)), alpha = 0.1) +
  geom_line(aes(x = ind_ano, y = growth_ind_ano_true), linewidth = 1) +
  labs(x = expression('Temperature'['individual-anomaly'] * ' (°C)'),
       y = "Growth",
       color = "Population") +
  scale_color_manual(values = c("Cool period" = "#313695", "Warm period" = "#a50026")) +
  theme(legend.position = c(0, 1.05),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill='transparent'),
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        legend.title = element_blank()) +
  guides(colour = guide_legend(override.aes = list(alpha = 1)))

p_growth_ind.ano <- last_plot()

## save plot ----
plot_list <- list("p_temp" = p_temp,
                  "p_temp_pop.ave" = p_temp_pop.ave,
                  "p_temp_pop.ano" = p_temp_pop.ano,
                  "p_temp_ind.ave" = p_temp_ind.ave,
                  "p_temp_ind.ano" = p_temp_ind.ano,
                  "p_growth_pop.ave" = p_growth_pop.ave,
                  "p_growth_pop.ano" = p_growth_pop.ano,
                  "p_growth_ind.ave" = p_growth_ind.ave,
                  "p_growth_ind.ano" = p_growth_ind.ano)

dir_report <- "./report"
for (i in 1:length(plot_list)) {
  ggsave(plot_list[[i]], file = file.path(dir_report, "fig3_conceptual temperature effect", paste0(names(plot_list[i]),".tiff") ),
         width = 3.5, height = 3.5,
         units = "cm",
         scaling = 0.4,
         dpi = 600,
         device = "tiff")
  
  ggsave(last_plot(), file = file.path(dir_report, "fig3_conceptual temperature effect", paste0(names(plot_list[i]),".eps")),
         width = 17, height = 5.75,
         units = "cm",
         dpi = 600, 
         device = "eps")
}



