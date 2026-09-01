






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






arkaute_no0 |> 
  select(treatment, plot, sampling, richness, biomass_mice_lm) |> 
  #filter(treatment %in% c("p", "wp")) |> 
  ggplot(aes(x = richness, y = biomass_mice_lm, color = treatment)) +
  facet_wrap(~treatment, ncol = 2, scales = "free") + 
  geom_point() +
  geom_smooth(method = "lm") 



arkaute_no0 |> 
  select(treatment, plot, sampling, richness, biomass_mice_lm) |> 
  filter(treatment %in% c("p", "wp")) |> 
  ggplot(aes(x = richness, y = biomass_mice_lm, color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm") 





arkaute_no0 |> 
  select(treatment, plot, sampling, richness, biomass_mice_lm) |> 
  pivot_longer(
    cols = c("richness", "biomass_mice_lm"), 
    names_to = "variable", 
    values_to = "value"
  ) |> 
  group_by(treatment, sampling, variable) |> 
  summarize(
    mean_sampling = mean(value, na.rm = T)
  ) |> 
  pivot_wider(
    id_cols = c("treatment", "sampling"), 
    names_from = "variable", 
    values_from = "mean_sampling"
  ) |> 
  rename(
    biomass_mean_sampling = biomass_mice_lm, 
    richness_mean_sampling = richness
  ) |> 
  #filter(treatment %in% c("p", "wp")) |> 
  ggplot(aes(x = richness_mean_sampling, y = biomass_mean_sampling, color = treatment)) +
  facet_wrap(~treatment, ncol = 2, scales = "free") + 
  geom_point() +
  geom_smooth(method = "lm") 





arkaute_no0 |> 
  select(treatment, plot, sampling, richness, biomass_mice_lm) |> 
  pivot_longer(
    cols = c("richness", "biomass_mice_lm"), 
    names_to = "variable", 
    values_to = "value"
  ) |> 
  group_by(treatment, sampling, variable) |> 
  summarize(
    mean_sampling = mean(value, na.rm = T)
  ) |> 
  pivot_wider(
    id_cols = c("treatment", "sampling"), 
    names_from = "variable", 
    values_from = "mean_sampling"
  ) |> 
  rename(
    biomass_mean_sampling = biomass_mice_lm, 
    richness_mean_sampling = richness
  ) |> 
  filter(treatment %in% c("p", "wp")) |> 
  ggplot(aes(x = richness_mean_sampling, y = biomass_mean_sampling, color = treatment)) +
  #facet_wrap(~treatment, ncol = 2, scales = "free") + 
  geom_point() +
  geom_smooth(method = "lm") 

















