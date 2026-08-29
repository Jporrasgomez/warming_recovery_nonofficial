


rm(list = ls(all.names = TRUE))
pacman::p_load(dplyr, tidyr, tidyverse, DT, viridis, ggrepel, codyn, vegan, eulerr, ggplot2, ggthemes, ggpubr, ggforce )#

source("code/palettes_labels.R")
sp_wide_plot <- read.csv("data/processed_data/relative_abudance_species_time.csv") 

spcomp <- sp_wide_plot %>% 
  select(-X) %>% 
  filter(sampling != "0") %>% 
  pivot_longer(5:44, values_to = "abundance", names_to = "code") %>% 
  mutate(abundance = na_if(abundance, 0)) %>%           # 0 -> NA
  group_by(treatment, code) %>% 
  summarise(mean_abundance = mean(abundance, na.rm = TRUE), .groups = "drop") %>% 
  filter(!is.nan(mean_abundance)) %>% 
  mutate(across(where(is.character), as.factor))


spcomp_c <- spcomp %>% 
  filter(treatment == "c") %>% 
  pull(code) %>% 
  droplevels() %>% 
  unique() %>% 
  print()

spcomp_w <- spcomp %>% 
  filter(treatment == "w") %>% 
  pull(code) %>% 
  droplevels() %>% 
  unique() %>% 
  print()

spcomp_p <- spcomp %>% 
  filter(treatment == "p") %>% 
  pull(code) %>% 
  droplevels() %>% 
  unique() %>% 
  print()

spcomp_wp <- spcomp %>% 
  filter(treatment == "wp") %>% 
  pull(code) %>% 
  droplevels() %>% 
  unique() %>% 
  print()



setdiff(spcomp_c, spcomp_w)
setdiff(spcomp_w, spcomp_c)

setdiff(spcomp_c, spcomp_p)
setdiff(spcomp_p, spcomp_c) # New species in perturbation treatment

setdiff(spcomp_c, spcomp_wp)
setdiff(spcomp_wp, spcomp_c) # New species in combined treatment

setdiff(spcomp_wp, spcomp_p)
setdiff(spcomp_p, spcomp_wp)

setdiff(spcomp_wp, spcomp_w)
setdiff(spcomp_w, spcomp_wp)




sp_colonizers <- 
sp_wide_plot %>% 
  select(-X) %>% 
  filter(treatment %in% c("p", "wp")) %>% 
  pivot_longer(5:44, values_to = "abundance", names_to = "code") %>% 
  mutate(abundance = na_if(abundance, 0)) %>% 
  filter(code %in% c("amsp", "casp", "chsp", "cisp", 
                     "kisp", "lapu", "mean", "cabu")) %>% 
  group_by(treatment, sampling, date, code) %>%
  summarize(mean_abundance = mean(abundance, na.rm = TRUE)) %>% 
  filter(!is.nan(mean_abundance)) %>% 
  mutate(across(where(is.character), as.factor)) %>% 
  mutate(code = fct_recode(code,
                           "Amaranthus sp."           = "amsp",
                           "Cardamine sp."            = "casp",
                           "Chenopodium sp."          = "chsp",
                           "Cirsium sp."              = "cisp",
                           "Kickxia spuria"           = "kisp",
                           "Lamium purpureum"         = "lapu",
                           "Mercurialis annua"        = "mean",
                           "Capsella bursa-pastoris"  = "cabu")) %>% 
  mutate(
    warming_excluded = as.factor(case_when(code %in% c("Amaranthus sp.", "Kickxia spuria") ~ "YES",
                                 .default = "NO"))
  ) %>% 
  as.data.frame()

  gg_colonizers <- 
  ggplot(sp_colonizers, aes(x = sampling, y = mean_abundance, color = code)) +
  facet_wrap(~ treatment,
             labeller = labeller(treatment = as_labeller(labels1)), ncol = 2) +
  geom_point(aes(shape = warming_excluded), size = 3, show.legend = FALSE) + 
  geom_path(aes(group = code), linewidth = 1) + 
  geom_vline(xintercept = 1.5, linetype = "dashed") +
  
   scale_color_brewer(
    palette = "Dark2",
    labels = c(
      "Amaranthus sp."          = expression(italic(Amaranthus~~sp.)),
      "Cardamine sp."           = expression(italic(Cardamine~~sp.)),
      "Chenopodium sp."         = expression(italic(Chenopodium~~sp.)),
      "Cirsium sp."             = expression(italic(Cirsium~~sp.)),
      "Kickxia spuria"          = expression(italic(Kickxia~~spuria)),
      "Lamium purpureum"        = expression(italic(Lamium~~purpureum)),
      "Mercurialis annua"       = expression(italic(Mercurialis~~annua)),
      "Capsella bursa-pastoris" = expression(italic(Capsella~~bursa-pastoris))
    )
  ) +
  
  labs(y = "Mean relative abundance (%)", 
       x = "Samplings", color = "Species", shape = NULL) +
  
  theme1 +
  theme(legend.position = "bottom")

print(gg_colonizers) # Supplementary Fig. 2



