






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
#View(false_cases_sampling)


#  Log Response Ratio


variables <- 
  (arkaute %>%  
     select(-date, -year, - date_label, -date_label_noyear, -sampling, -plot, -treatment,
            -OTC, -perturbation, -mean_temperature, -mean_vwc
     ) %>% 
     colnames()
  ) 

# Charging functions
source("code/functions/eff_size_LRR_function_delta.R")      # Function of LRR at aggregated level CORRECTED
source("code/functions/new_dynamics.R")               # Function of LRR at dynamics level
source("code/functions/gg_aggregated_function_2.R")   # Function for visualization of aggregated analysis
source("code/functions/gg_dynamics_function2.R")      # Function for visualization of dynamics analysis


k = 1   
# k = 1: To see all (but biomass raw and biomass LM)
# k = 2: biomass variables for sensitivity analysis
# k = 3: Richness, abundance and biomass for warming effects on recovery()
# K = 4: Evenness and functional traits for wp vs p

width_dynamics = 3 # size for plots (1:3) 1 for agg analysis, 3 for dynamics

{
  
  # 1. Log Response Analysis: 
  
  list_arkaute <- list(arkaute_no0, arkaute)
  list_agg <- list()
  list_dyn <- list()
  
  for(i in seq_along(variables)){ 
    
    
    LRR_agg(list_arkaute[[1]], variables[i])
    list_agg[[i]] <- effsize_data 
    
    LRR_dynamics(list_arkaute[[2]], variables[i])
    list_dyn[[i]] <- effsize_dynamics_data
    
  }
  
  agg <- do.call(rbind, list_agg) %>% 
    select(eff_descriptor, variable, scale, eff_value, lower_limit, upper_limit, null_effect)
  
 
  
  agg %>%   write.csv("results/effect_size_aggregated.csv")

  agg_glmm  <-  read.csv("results/GLMM_agg.csv") |> select(-X) |> 
    mutate(
      variable = fct_recode(variable,
                            "Y_zipf" = "evenness",
                            "biomass_mice_lm"  = "biomass"))
    
  
  agg <- merge(agg, agg_glmm) |> 
    select(eff_descriptor, variable, eff_value, lower_limit, upper_limit,
           null_effect, scale,  estimate, SE, p_value, AIC, estimate_type, effect_sign, 
           effect_significance) |> 
    rename(glmm_estimate = estimate, 
           glmm_SE = SE, 
           glmm_p_value = p_value, 
           glmm_AIC = AIC, 
           glmm_estimatetype = estimate_type, 
           glmm_effect_sign = effect_sign,
           glmm_effect_significance = effect_significance) |> 
    mutate(
      final_effect = case_when(
        eff_value > 0 & glmm_effect_sign == "positive"      ~ paste0("positive"), 
        eff_value > 0 & glmm_effect_sign == "negative"      ~ paste0("conflict"),
        eff_value < 0 & glmm_effect_sign == "negative"      ~ paste0("negative"),
        eff_value < 0 & glmm_effect_sign == "positive"      ~ paste0("conflict"),
      ),
      final_significance = case_when(
        null_effect == "NO"  & glmm_effect_significance == "significant"      ~ paste0("both"),
        null_effect == "NO"  & glmm_effect_significance == "marginal"         ~ paste0("both-marginal"),
        null_effect == "NO"  & glmm_effect_significance == "non-significant"  ~ paste0("lost"),
        null_effect == "YES" & glmm_effect_significance == "significant"      ~ paste0("conflict"), 
        null_effect == "YES" & glmm_effect_significance == "marginal"         ~ paste0("conflict-marginal"),
        null_effect == "YES" & glmm_effect_significance == "non-significant"  ~ paste0("none")
      )
    )
    
  
  dyn <- do.call(rbind, list_dyn) %>% 
    mutate(
      date_label_noyear = factor(
        date_label_noyear,
        levels = unique(date_label_noyear[order(date)]),
        ordered = TRUE
      )
    )
   
  dyn %>%   write.csv("results/effect_size_dynamics.csv")
  
  # We need samplings 0 and 1, for which GLMM model were not fitted. 
  dyn_01 <- dyn |> 
    select(sampling, year, date, date_label_noyear,  eff_descriptor, variable, eff_value,
           lower_limit, upper_limit, null_effect, scale) |> 
    filter(sampling %in% c("0", "1")) |> 
    mutate(
      glmm_estimate = NA, glmm_SE = NA, glmm_p_value = NA, glmm_AIC = NA, glmm_estimate_type = NA,
      glmm_effect_sign = NA, glmm_effect_significance = NA
    )
    
  
  dyn_glmm  <-  read.csv("results/GLMM_dyn.csv") |> select(-X) |> 
    mutate(
      variable = fct_recode(variable,
                            "Y_zipf" = "evenness",
                            "biomass_mice_lm"  = "biomass"),
      sampling = as.factor(sampling))
  
  dyn <- left_join(dyn, dyn_glmm) |> 
    select(sampling, year, date, date_label_noyear,  eff_descriptor, variable, eff_value,
           lower_limit, upper_limit, null_effect, scale, estimate, SE, p_value, AIC, estimate_type,
           effect_sign, effect_significance) |> 
    rename(glmm_estimate = estimate, 
           glmm_SE = SE, 
           glmm_p_value = p_value, 
           glmm_AIC = AIC, 
           glmm_estimatetype = estimate_type, 
           glmm_effect_sign = effect_sign,
           glmm_effect_significance = effect_significance) |> 
    mutate(
      final_effect = case_when(
        eff_value > 0 & glmm_effect_sign == "positive"      ~ paste0("positive"), 
        eff_value > 0 & glmm_effect_sign == "negative"      ~ paste0("conflict"),
        eff_value < 0 & glmm_effect_sign == "negative"      ~ paste0("negative"),
        eff_value < 0 & glmm_effect_sign == "positive"      ~ paste0("conflict"),
      ),
      final_significance = case_when(
        null_effect == "NO"  & glmm_effect_significance == "significant"      ~ paste0("both"),
        null_effect == "NO"  & glmm_effect_significance == "marginal"         ~ paste0("both-marginal"),
        null_effect == "NO"  & glmm_effect_significance == "non-significant"  ~ paste0("lost"),
        null_effect == "YES" & glmm_effect_significance == "significant"      ~ paste0("conflict"), 
        null_effect == "YES" & glmm_effect_significance == "marginal"         ~ paste0("conflict-marginal"),
        null_effect == "YES" & glmm_effect_significance == "non-significant"  ~ paste0("none")
      )
    ) |> 
    bind_rows(dyn_01) # Integration of samplings 0 and 1

  
  ## 2. GENERATING PLOTS ####
  
  
  limits_main_variables <- c("richness",                    # 1     
                             "abundance",                        # 2     
                             "Y_zipf",                           # 3     
                             "SLA",                              # 4    
                             "LDMC",                             # 5    
                             "leafN",                            # 6
                             "biomass_mice_lm"                   # 7
  )
  
  labels_main_variables <- c("richness" = "Richness",            # 1
                             "abundance" = "Cover",                   # 2
                             "Y_zipf" = "Evenness",                   # 3
                             "SLA" = "SLA",                           # 4
                             "LDMC" = "LDMC",                         # 5
                             "leafN"= "LN",                           # 6
                             "biomass_mice_lm" = "Biomass"    # 7
  )    
  
  
  limits_biomass_variables <- c(
    "biomass_raw",      # Biomass data without imputation of NA                   
    "biomass_mice",     # Biomass data with mice imputation of NA
    "biomass_mice_lm"   # Biomass data with mice imputation and LM imputation                  
  )
  
  labels_biomass_variables <- c(
    
    "biomass_raw" = "No imputation",
    "biomass_mice" = "MICE",
    "biomass_mice_lm" = "MICE + LM"
  )      
  
  
  
  limits_warming_recovery <- c("richness",                        
                               "abundance",                            
                               "biomass_mice_lm"                   
  )
  
  labels_warming_recovery <- c("richness" = "Richness",            
                             "abundance" = "Cover",                 
                             "biomass_mice_lm" = "Biomass"    
  )    
  
  
  
  limits_other_variables_wp <- c(    
                             "Y_zipf",                           # 3     
                             "SLA",                              # 4    
                             "LDMC",                             # 5    
                             "leafN"                           # 6
                             
  )
  
  labels_other_variables_wp <- c(
                             "Y_zipf" = "Evenness",                   # 3
                             "SLA" = "SLA",                           # 4
                             "LDMC" = "LDMC",                         # 5
                             "leafN"= "LN"                         # 6
                            
  )    
  
  limits_list <- list(limits_main_variables, limits_biomass_variables, limits_warming_recovery, limits_other_variables_wp)
  
  limits_variables <- limits_list[[k]]
  
  labels_list <- list(labels_main_variables, labels_biomass_variables, labels_warming_recovery, labels_other_variables_wp)
  
  labels_variables <- labels_list[[k]]
  
  
  lvls <- limits_variables
  labs <- unname(labels_variables[lvls])
  
  
  
  # PERTURBATION AND WARMING TREATMENTS / CONTROL #
  
  
  comparissons <- c("p_vs_c", "w_vs_c")
  
  pos_dod_c_agg <- position_dodge2(width = 0.3, preserve = "single")
  pos_dod_c_dyn <- position_dodge2(width = 12, preserve = "single")
  
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
      limitvar  = lvls,
      labelvar  = labels_variables[lvls], 
      breaks_axix_y = 3
    )
  
  
  
  gg_eff_dynamics_c2<- dyn %>% 
    filter(eff_descriptor %in% comparissons) %>% 
    filter(variable %in% limits_variables) %>%  
    mutate(
      variable = factor(variable, 
                        levels = limits_variables, 
                        labels = labels_variables)) %>% 
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
                   widths = c(1, width_dynamics))) +
    plot_annotation(theme = theme(legend.position = "bottom"))
  
  
  
  #  COMBINED / PERTURBATION    ###
  
  pos_dod_wp_agg <- position_dodge2(width = 0.1, preserve = "single")
  pos_dod_wp_dyn <- position_dodge2(width = 4, preserve = "single")
  
  
  
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
      limitvar  = lvls,
      labelvar  = labels_variables[lvls], 
      breaks_axix_y = 3
    )
  
  
  gg_eff_dynamics_wp2<- dyn %>% 
    filter(eff_descriptor %in% c("wp_vs_p")) %>% 
    filter(variable %in% limits_variables) %>%  
    mutate(
      variable = factor(variable, 
                        levels = limits_variables, 
                        labels = labels_variables)) %>% 
    
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
                   widths = c(1, width_dynamics))) +
    plot_annotation(theme = theme(legend.position = "bottom"))
  
}


print(gg_control)
print(gg_wp) 



ggsave("results/Figure_2.png", plot = gg_control, dpi = 600)
ggsave("results/Figure_2.svg", plot = gg_control, dpi = 600)
ggsave("results/Figure_3.png", plot = gg_wp, dpi = 600)
ggsave("results/Figure_3.svg", plot = gg_wp, dpi = 600)



# Biomass sensitivity

ggsave("results/Supplementary_figure_biomass_c.png", plot = gg_control, dpi = 600)
ggsave("results/Supplementary_figure_biomass_c.svg", plot = gg_control, dpi = 600)
ggsave("results/Supplementary_figure_biomass_wp.png", plot = gg_wp, dpi = 600)
ggsave("results/Supplementary_figure_biomass_wp.svg", plot = gg_wp, dpi = 600)


# Other variables WP vs P

ggsave("results/Sup_fig_WPP.png", plot = gg_wp, dpi = 600)
ggsave("results/Sup_fig_WPP.svg", plot = gg_wp, dpi = 600)



