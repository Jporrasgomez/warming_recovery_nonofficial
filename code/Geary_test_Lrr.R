

rm(list = ls(all.names = TRUE))  #Se limpia el environment
pacman::p_unload(pacman::p_loaded(), character.only = TRUE) #


pacman::p_load(dplyr,reshape2,tidyverse, lubridate, ggplot2, ggpubr, gridExtra,
               car, ggsignif, dunn.test, rstatix, ggbreak, effsize, patchwork) 

source("code/palettes_labels.R")



arkaute <- read.csv("data/processed_data/arkaute.csv") %>% 
  mutate(
    year = as.factor(year),
    date = ymd(date),
    sampling = as.factor(sampling),
    plot = as.factor(plot),
    treatment = as.factor(treatment))  

arkaute_no0 <- arkaute %>% 
  filter(sampling != "0")

  
  
  ## Geary test to test suitability of data for LRR analysis
  ## We use a modification of Geary Test proposed by Lajeunesse 2015
  
  # 1. Geary test at Treatment level 
  geary_test_treatment0 <- arkaute_no0 %>%
    pivot_longer(
      cols = c(-date, -year, - date_label, -date_label_noyear, -sampling, -plot, -treatment,
               -OTC, -perturbation),
      values_to = "value", 
      names_to = "variable"
    ) %>% 
    mutate(value = ifelse(variable == "Y_zipf", value * -1, value)) %>%
    group_by(treatment, variable, plot) %>% 
    summarize(
      mean_plot = mean(value, na.rm = T), 
    )
  
  geary_test_treatment <- geary_test_treatment0 |> 
    group_by(treatment, variable) %>% 
    summarise(
      mean_variable = mean(mean_plot),
      sd_variable = sd(mean_plot),
      n = n(),
      .groups = "drop"
    ) |> 
    mutate(geary_test_value = (mean_variable/sd_variable) *((4 * n^1.5) / (1 + 4 * n))) %>% 
    mutate(geary_test_outcome = ifelse(geary_test_value >= 3, paste0("TRUE"), paste0("FALSE")))
  
  # 2. Geary test at sampling level 
  geary_test_sampling <- arkaute %>%
    pivot_longer(
      cols = c(-date, -year, - date_label, -date_label_noyear, -sampling, -plot, -treatment,
               -OTC, -perturbation),
      values_to = "value", 
      names_to = "variable"
    ) %>% 
    mutate(value = ifelse(variable == "Y_zipf", value * -1, value)) %>%
    group_by(treatment, sampling, variable) %>% 
    summarize(
      n = n(),
      mean_variable = mean(value, na.rm = T), 
      sd_variable = sd(value, na.rm = T)
    ) %>% 
    mutate(geary_test_value = (mean_variable/sd_variable) *((4 * n^1.5) / (1 + 4 * n))) %>% 
    mutate(geary_test_outcome = ifelse(geary_test_value >= 3, paste0("TRUE"), paste0("FALSE")))
  
  false_sampling <- geary_test_sampling %>% 
    filter(geary_test_outcome == "FALSE")
  unique(false_sampling$variable)
  
  nrow(false_sampling)/nrow(geary_test_sampling) *100
  
  false_cases_sampling <- false_sampling %>%
    ungroup() %>% 
    select(variable, treatment, sampling, geary_test_value, geary_test_outcome) %>% 
    group_by(variable, treatment) %>% 
    summarize(
      n = n(),
      samplings = list(as.numeric(sampling)),
      .groups = "drop"
    )
  
  
  
#View(geary_test_treatment)
#View(false_case_samplings)



