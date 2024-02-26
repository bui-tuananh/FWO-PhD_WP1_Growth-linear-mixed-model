
# setup ----
#install.packages("icesSAG")
library(icesSAG)
library(tidyverse)
library(readxl)

# get stock assessment data ----
meta <- getListStocks(2023)
sol <- meta %>% filter(Purpose == "Advice",
                StockKeyLabel %in% c("sol.27.4", "sol.27.7a", "sol.27.8ab"))
# assessmentKey 
# North Sea     sol.27.4    18458
# Irish Sea     sol.27.7a   18059
# Bay of Biscay sol.27.8ab  18024

stock_4 <- getSummaryTable(assessmentKey = 18458)[[1]] %>%
  mutate(pop = "4bc") %>%
  rename(year = Year)

stock_7a <- getSummaryTable(assessmentKey = 18059)[[1]] %>%
  mutate(pop = "7a") %>%
  rename(year = Year)

stock_8ab <- getSummaryTable(assessmentKey = 18024)[[1]] %>%
  mutate(pop = "8ab") %>%
  rename(year = Year)

# merge all data
stock_all <- rbind(stock_4, stock_7a, stock_8ab)

# save data
write_rds(stock_all, file.path("./ICES/stock assessment 2023", "stock-assessment_2023.rds"))

# process abundance at age data ----
dir_ices_stock <- "./ICES/stock assessment 2023"

abun_4 <- read_excel(file.path(dir_ices_stock, "abundances at age (thousands)_4.xlsx"))
abun_4 <- abun_4 %>% 
  pivot_longer(
    cols = 2:ncol(abun_4),
    names_to = "age",
    values_to = "abundance"
    ) %>%
  mutate(
    unit = "thousdands",
    pop = "4bc"
  )

abun_7a <- read_excel(file.path(dir_ices_stock, "abundances at age (thousands)_7a.xlsx")) 
abun_7a <- abun_7a %>% 
  pivot_longer(
    cols = 2:ncol(abun_7a),
    names_to = "age",
    values_to = "abundance"
  ) %>%
  mutate(
    age = if_else(age == "8+", 8, as.numeric(age)),
    unit = "thousdands",
    pop = "7a"
  )

abun_8ab <- read_excel(file.path(dir_ices_stock, "abundances at age (thousands)_8ab.xlsx"))
abun_8ab <- abun_8ab %>% 
  pivot_longer(
    cols = 2:ncol(abun_8ab),
    names_to = "age",
    values_to = "abundance"
  ) %>%
  mutate(
    age = if_else(age == "8+", 8, as.numeric(age)),
    unit = "thousdands",
    pop = "8ab"
  )

# merge all data
abun_all <- rbind(abun_4, abun_7a, abun_8ab)

# save data
write_rds(abun_all, file.path("./ICES/stock assessment 2023", "abundance-at-age_2023.rds"))
