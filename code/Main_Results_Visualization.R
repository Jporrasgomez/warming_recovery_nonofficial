

rm(list = ls(all.names = TRUE))  #Se limpia el environment
pacman::p_unload(pacman::p_loaded(), character.only = TRUE) #


pacman::p_load(dplyr,reshape2,tidyverse, lubridate, ggplot2, ggpubr, gridExtra,
               car, ggsignif, dunn.test, rstatix, ggbreak, effsize, patchwork) 

source("code/palettes_labels.R")


## JOINING GLMM RESULTS WITH LRR RESULTS ########################################################################


agg_glmm  <-  read.csv("results/GLMM_agg.csv") |> select(-X) |> 
  mutate(
    variable = fct_recode(variable,
                          "Y_zipf" = "evenness",
                          "biomass_mice_lm"  = "biomass"))

agg_LRR <- read.csv("results/LRR_results_agg.csv") |>  select(-X)


agg <- merge(agg_LRR, agg_glmm) |> 
  select(eff_descriptor, variable, eff_value, lower_limit, upper_limit,
         null_effect, scale, contrast,  estimate, SE, z.ratio,  p_value, AIC, estimate_type, effect_sign, 
         effect_significance) |> 
  rename(glmm_contrast = contrast,
         glmm_estimate = estimate, 
         glmm_SE = SE, 
         glmm_z.ratio = z.ratio,
         glmm_p_value = p_value, 
         glmm_AIC = AIC, 
         glmm_estimatetype = estimate_type, 
         glmm_effect_sign = effect_sign,
         glmm_effect_significance = effect_significance)




dyn_glmm  <-  read.csv("results/GLMM_dyn.csv") |> select(-X) |> 
  mutate(
    variable = fct_recode(variable,
                          "Y_zipf" = "evenness",
                          "biomass_mice_lm"  = "biomass"),
    sampling = as.factor(sampling))

dyn_LRR <- read.csv("results/LRR_results_dyn.csv") |>  select(-X) |> mutate(sampling = as.factor(sampling))

dyn <- left_join(dyn_LRR, dyn_glmm) |> 
  select(sampling, year, date, date_label_noyear,  eff_descriptor, variable, eff_value,
         lower_limit, upper_limit, null_effect, scale, contrast, estimate, SE, z.ratio, p_value, AIC, estimate_type,
         effect_sign, effect_significance) |> 
  rename(glmm_contrast = contrast,
         glmm_estimate = estimate, 
         glmm_SE = SE, 
         glmm_z.ratio = z.ratio,
         glmm_p_value = p_value, 
         glmm_AIC = AIC, 
         glmm_estimatetype = estimate_type, 
         glmm_effect_sign = effect_sign,
         glmm_effect_significance = effect_significance) |> 
  filter(eff_descriptor != "wp_vs_w") |> 
  mutate(
    year = as.factor(year),
    date = ymd(date),
    sampling = as.factor(sampling),
    date_label_noyear = factor(
      date_label_noyear,
      levels = unique(date_label_noyear[order(date)]), 
      ordered = TRUE
    ))




## 2. GENERATING PLOTS ##############################################################################################################

source("code/functions/gg_aggregated_function_2.R")   # Function for aggregated analysis plots
source("code/functions/gg_dynamics_function2.R")      # Function for temporal dynamics plots

k = 1  
# k = 1: To see all variable (but biomass raw and biomass LM)
# k = 2: biomass variables for sensitivity analysis
# k = 3: Richness, abundance and biomass for warming effects on recovery()
# K = 4: Evenness and functional traits for wp vs p

{
# When k = 1
labels_main_variables <- c("richness" = "Richness",            # 1
                           "abundance" = "Cover",              # 2     
                           "Y_zipf" = "Evenness",              # 3     
                           "SLA" = "SLA",                      # 4     
                           "LDMC" = "LDMC",                    # 5     
                           "leafN"= "LN",                      # 6     
                           "biomass_mice_lm" = "Biomass"       # 7
)    

# When k = 2
labels_biomass_variables <- c("biomass_raw" = "No imputation", "biomass_mice" = "MICE",
                              "biomass_mice_lm" = "MICE + LM")      

# When k = 3
labels_warming_recovery <- labels_main_variables[c(1, 2, 7)]

# When k = 4
labels_other_variables_wp <- labels_main_variables[c(2,3,4,5)]

labels_list <- list(labels_main_variables, labels_biomass_variables, labels_warming_recovery, labels_other_variables_wp)

lvls <- names(labels_list[[k]])
labs <- unname(labels_list[[k]])



# PERTURBATION AND WARMING TREATMENTS / CONTROL ################

comparissons <- c("p_vs_c", "w_vs_c")


# Aggregated analysis
gg_eff_agg_c2 <- agg %>% 
  filter(
    eff_descriptor %in% comparissons,
    variable %in% lvls
  ) %>% 
  mutate(
    eff_descriptor = factor(eff_descriptor, levels = comparissons),
    variable       = factor(variable, levels = lvls, labels = labs)
  ) %>% 
  ggagg2(
    palette   = palette_RR_CB,
    labels    = labels_RR2,
    colorline = "grey50",
    breaks_axix_y = 3
  )

# Temporal dynamics

pos_dod_c_dyn <- position_dodge2(width = 12, preserve = "single") 
gg_eff_dynamics_c2<- dyn %>% 
  filter(eff_descriptor %in% comparissons) %>% 
  filter(variable %in% lvls) %>%  
  mutate(
    variable = factor(variable, 
                      levels = lvls, 
                      labels = labs)) %>% 
  ggdyn2(palette_RR_CB,
         labels_RR2, 
         "grey50",
         position = position_dodge(width = 0.5),
         asterisk = 8, 
         caps = position_dodge(width = 0.5)$width,
         breaks_axix_y = 3)


gg_control <-
  (gg_eff_agg_c2 + 
     gg_eff_dynamics_c2 + theme (legend.position = "none") + 
     plot_layout(guides = "collect",
                 widths = c(1, 3))) +
  plot_annotation(theme = theme(legend.position = "bottom"))



#  COMBINED / PERTURBATION    ##################################

# Aggregated analysis
gg_eff_agg_wp2 <- agg %>% 
  filter(eff_descriptor == "wp_vs_p",
         variable %in% lvls) %>% 
  mutate(
    variable = factor(variable, levels = lvls, labels = labs)
  ) %>% 
  ggagg2(
    palette   = palette_RR_wp,
    labels    = labels_RR_wp,
    colorline = p_CB,
    breaks_axix_y = 3
  )


# Temporal dynamics
pos_dod_wp_dyn <- position_dodge2(width = 4, preserve = "single")

gg_eff_dynamics_wp2<- dyn %>% 
  filter(eff_descriptor %in% c("wp_vs_p")) %>% 
  filter(variable %in% lvls) %>%  
  mutate(
    variable = factor(variable, 
                      levels = lvls, 
                      labels = labs)) %>% 
  
  ggdyn2(palette_RR_wp,
         labels_RR_wp2, 
         p_CB,
         position = position_dodge(width = 0.5),
         asterisk = 8, 
         caps = position_dodge(width = 0.5)$width,
         breaks_axix_y = 3)


gg_wp <-
  (gg_eff_agg_wp2 + 
     gg_eff_dynamics_wp2 + theme (legend.position = "none") +
     plot_layout(guides = "collect",
                 widths = c(1, 3))) +
  plot_annotation(theme = theme(legend.position = "bottom"))

}



print(gg_control)
print(gg_wp) 


# All variables
#
#ggsave("results/Figure_2.png", plot = gg_control, dpi = 600)
#ggsave("results/Figure_2.svg", plot = gg_control, dpi = 600)
#ggsave("results/Figure_3.png", plot = gg_wp, dpi = 600)
#ggsave("results/Figure_3.svg", plot = gg_wp, dpi = 600)
#
#
#
## Biomass sensitivity
#
#ggsave("results/Supplementary_figure_biomass_c.png", plot = gg_control, dpi = 600)
#ggsave("results/Supplementary_figure_biomass_c.svg", plot = gg_control, dpi = 600)
#ggsave("results/Supplementary_figure_biomass_wp.png", plot = gg_wp, dpi = 600)
#ggsave("results/Supplementary_figure_biomass_wp.svg", plot = gg_wp, dpi = 600)
#
#
#
## Other variables WP vs P
#
#ggsave("results/Sup_fig_WPP.png", plot = gg_wp, dpi = 600)
#ggsave("results/Sup_fig_WPP.svg", plot = gg_wp, dpi = 600)
#


