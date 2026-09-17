






rm(list = ls(all.names = TRUE))  #Se limpia el environment
pacman::p_unload(pacman::p_loaded(), character.only = TRUE) #


pacman::p_load(dplyr,reshape2,tidyverse, lubridate, ggplot2, ggpubr, gridExtra,
               car, ggsignif, dunn.test, rstatix, ggbreak, effsize, patchwork) 

source("code/palettes_labels.R")


theme_set(
  theme_bw() +
    theme(
      legend.position   = "right",
      panel.grid        = element_blank(),
      strip.background  = element_blank(),
      strip.text        = element_text(face = "bold"),
      text              = element_text(size = 11)
    ))

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

library(ggpmisc) # Paquete para añadir ecuaciones y métricas estadísticas al gráfico

arkaute_no0 |> 
  # Simplificación del cálculo de promedios sin necesidad de pivotar
  group_by(treatment, sampling) |> 
  summarize(
    richness_mean_sampling = mean(richness, na.rm = TRUE),
    biomass_mean_sampling = mean(biomass_mice_lm, na.rm = TRUE),
    .groups = "drop"
  ) |> 
  ggplot(aes(x = richness_mean_sampling, y = biomass_mean_sampling, color = treatment, fill = treatment)) +
  facet_wrap(~treatment, ncol = 4, scales = "free_x") +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  # Añade ecuación, R2 y p-valor automáticamente en cada facet
  stat_poly_eq(
    use_label("eq"),
    formula = y ~ x,
    label.x = "left", label.y = 0.95
  ) +
  # Capa 2: R2 (centro)
  stat_poly_eq(
    use_label("R2"),
    formula = y ~ x,
    label.x = "left", label.y = 0.90
  ) +
  # Capa 3: p-value (abajo)
  stat_poly_eq(
    use_label("p"),
    formula = y ~ x,
    label.x = "left", label.y = 0.85
  ) +
  scale_color_manual(values = palette_CB, labels = labels1) +
  scale_fill_manual(values = palette_CB, labels = labels1)+
  labs(x = "Mean richness at sampling level", y = "Mean biomass at sampling level") +
  theme(legend.position = "bottom")
  





arkaute_no0 |> 
  ggplot(aes(x = richness, y = biomass_mice_lm, color = treatment, fill = treatment)) +
  facet_wrap(~treatment, ncol = 4, scales = "free_x") +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  # Capa 1: Ecuación (arriba)
  stat_poly_eq(
    use_label("eq"),
    formula = y ~ x,
    label.x = "left", label.y = 0.95
  ) +
  # Capa 2: R2 (centro)
  stat_poly_eq(
    use_label("R2"),
    formula = y ~ x,
    label.x = "left", label.y = 0.90
  ) +
  # Capa 3: p-value (abajo)
  stat_poly_eq(
    use_label("p"),
    formula = y ~ x,
    label.x = "left", label.y = 0.85
  ) +
  scale_color_manual(values = palette_CB, labels = labels1) +
  scale_fill_manual(values = palette_CB, labels = labels1) +
  labs(x = "Richness at plot level", y = "Biomass at plot level") +
  theme(legend.position = "bottom")



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

















