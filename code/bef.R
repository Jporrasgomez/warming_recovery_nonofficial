






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





library(ggpmisc) 


BEF_sampling <- 
arkaute_no0 |> 
  group_by(treatment, sampling) |> 
  summarize(
    richness_mean_sampling = mean(richness, na.rm = TRUE),
    biomass_mean_sampling = mean(biomass_mice_lm, na.rm = TRUE),
    .groups = "drop"
  ) |> 
  ggplot(aes(x = richness_mean_sampling, y = biomass_mean_sampling, color = treatment, fill = treatment)) +
  facet_wrap(~treatment, ncol = 4, scales = "free_x",
             labeller = as_labeller(c(
               c     = "Control",
               p    = "Perturbed-only",
               w = "Warmed-only",
               wp = "Combined"))) +
  geom_point(size = 4, alpha = 0.8) +
  geom_smooth(method = "lm", formula = y ~ x) +
  stat_poly_eq(
    use_label("eq"),
    formula = y ~ x,
    label.x = "left", label.y = 0.95,
    color = "black"
  ) +
  stat_poly_eq(
    use_label("R2"),
    formula = y ~ x,
    label.x = "left", label.y = 0.90,
    color = "black"
  ) +
  stat_poly_eq(
    use_label("p"),
    formula = y ~ x,
    label.x = "left", label.y = 0.85,
    color = "black"
  ) +
  scale_color_manual(values = palette_CB, labels = labels1) +
  scale_fill_manual(values = palette_CB, labels = labels1)+
  labs(x = NULL, y = "Biomass") +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = "black"),
    text = element_text(size = 14),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 14),
    strip.placement = "outside",                
    strip.text.y.left = element_text(            
      angle = 90, face = "bold", size = 14
    ),
    axis.text.y = element_text(hjust = 0.5, face = "plain", size = 12),
    axis.text.x = element_text(face = "plain", size = 12),
    legend.position = " none",
    legend.text = element_text(size = 14, face = "plain"),
    legend.key = element_rect(fill = "white", colour = NA)
  )
  





BEF_plot <- 
arkaute_no0 |> 
  ggplot(aes(x = richness, y = biomass_mice_lm, color = treatment, fill = treatment)) +
  facet_wrap(~treatment, ncol = 4, scales = "free_x",
             labeller = as_labeller(c(
               c     = "Control",
               p    = "Perturbed-only",
               w = "Warmed-only",
               wp = "Combined"))) +
  geom_point(size = 2, alpha = 0.6) +
  geom_smooth(method = "lm", formula = y ~ x) +
  stat_poly_eq(
    use_label("eq"),
    formula = y ~ x,
    label.x = "left", label.y = 0.95,
    color = "black"
  ) +
  stat_poly_eq(
    use_label("R2"),
    formula = y ~ x,
    label.x = "left", label.y = 0.90,
    color = "black"
  ) +
  stat_poly_eq(
    use_label("p"),
    formula = y ~ x,
    label.x = "left", label.y = 0.85,
    color = "black"
  ) +
  scale_color_manual(values = palette_CB, labels = labels1) +
  scale_fill_manual(values = palette_CB, labels = labels1)+
  labs(x = "Richness", y = "Biomass") +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = "black"),
    text = element_text(size = 14),
    strip.background = element_blank(),
    strip.text = element_blank(),
    strip.placement = "outside",                
    strip.text.y.left = element_text(            
      angle = 90, face = "bold", size = 14
    ),
    axis.text.y = element_text(hjust = 0.5, face = "plain", size = 12),
    axis.text.x = element_text(face = "plain", size = 12),
    legend.position = " none",
    legend.text = element_text(size = 14, face = "plain"),
    legend.key = element_rect(fill = "white", colour = NA)
  )

#Checking
arkaute_no0 |> 
  filter(treatment == "w") |> 
  lm(biomass_mice_lm ~ richness, data = _) |> 
  summary()



#print(BEF_sampling) 
#print(BEF_plot)

gg_BEF <- 
BEF_sampling + 
  BEF_plot + theme (legend.position = "none") +
  plot_layout(guides = "collect", ncol = 1) +
  plot_annotation(theme = theme(legend.position = "bottom"),
                  tag_levels = 'A',)

print(gg_BEF)



ggsave("results/BEF.png", plot = gg_BEF, dpi = 600)












