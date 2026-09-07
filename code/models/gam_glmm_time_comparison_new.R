


rm(list = ls(all.names = TRUE))
pacman::p_load(
  dplyr, reshape2, tidyverse, lubridate, ggplot2, ggpubr, gridExtra, ggrepel,
  car, ggsignif, dunn.test, rstatix, performance, emmeans, DHARMa, mgcv, glmmTMB
)

source("code/palettes_labels.R")


#{
  
  theme_set(
    theme_bw() +
      theme(
        legend.position   = "right",
        panel.grid        = element_blank(),
        strip.background  = element_blank(),
        strip.text        = element_text(face = "bold"),
        text              = element_text(size = 11)
      )
  )
  
  arkaute <- read.csv("data/processed_data/arkaute.csv") %>%
    mutate(
      year      = factor(year),
      date      = ymd(date),
      sampling_num = as.numeric(as.character(sampling)), # Necesario para splines en GAM
      sampling = factor(sampling, levels = as.character(sort(unique(as.numeric(as.character(sampling)))))),
      sampling_f = numFactor(sampling_num),
      plot      = factor(plot),
      treatment = factor(treatment)
    ) %>%
    filter(sampling != "0",
           sampling != "1") |> 
    arrange(plot, sampling)  # Necessary for autocorrelation in models
  
  
  arkaute_richness  <- arkaute |> filter(!is.na(richness))
  arkaute_abundance <- arkaute |> filter(!is.na(abundance))
  arkaute_evenness  <- arkaute |> filter(!is.na(Y_zipf))
  arkaute_sla       <- arkaute |> filter(!is.na(SLA))
  arkaute_ldmc      <- arkaute |> filter(!is.na(LDMC))
  arkaute_leafN     <- arkaute |> filter(!is.na(leafN))
  arkaute_biomass   <- arkaute |> filter(!is.na(biomass_mice_lm))
  
  

  ## Defininf JOINING AND CLEANING functions
  common_function <- function(data){
    data |> 
      filter(!contrast %in% c("p - w", "p / w", "w - wp", "w / wp")) |> 
      rename(p_value = p.value) |> 
      mutate(
        effect_sign = case_when(
          estimate_type == "substract" & estimate < 0 ~ "positive",
          estimate_type == "substract" & estimate > 0 ~ "negative",
          estimate_type == "ratio" & estimate > 1     ~ "negative", 
          estimate_type == "ratio" & estimate < 1     ~ "positive"
        ),
        effect_significance = case_when(
          p_value < 0.05                     ~ "significant",
          p_value >= 0.05  & p_value < 0.10  ~ "marginal",
          TRUE                               ~ "non-significant"
        ),
        eff_descriptor = case_when(
          contrast %in% c("c - p", "c / p")    ~ "p_vs_c", 
          contrast %in% c("c - w" , "c / w")   ~ "w_vs_c", 
          contrast %in% c("c - wp", "c / wp")  ~ "wp_vs_c", 
          contrast %in% c("p - wp", "p / wp")  ~ "wp_vs_p"
        ))
  }
  
  lrr_reading <- function(data) {
    data |> 
      select(-scale, -X) |> 
      rename(estimate = eff_value) |> 
      mutate(model = paste0("LRR"),
             variable = as.factor(variable)) |> 
      filter(!variable %in% c("biomass_raw", "biomass_mice")) %>%
      mutate(
        variable = fct_recode(variable,
                              "evenness" = "Y_zipf",
                              "biomass"  = "biomass_mice_lm"
        ),
        variable = droplevels(variable),
        effect_sign = ifelse(estimate > 0 , "positive", "negative"),
        effect_significance = case_when(
          null_effect == "YES" ~ "non-significant",
          TRUE                 ~ "significant", 
        )
      ) |> 
      select(-upper_limit, -lower_limit, -null_effect)
  }
  
  ## Opening LOG RESPONSE RATIO results datasets
  lrr_table <- read.csv("results/effect_size_aggregated.csv") |> 
    lrr_reading()
  lrr_table_dyn <- read.csv("results/effect_size_dynamics.csv") |> 
    filter(sampling != "0") |> 
    mutate(sampling = as.factor(sampling)) |> 
    lrr_reading()
  
  
  ############# GLMM MODELS ##############################################################
  ## Richness ## 
  glmm_richness <- glmmTMB( richness ~ treatment * sampling + ar1(sampling + 0 | plot),
                            dispformula = ~ treatment, 
                            data = arkaute_richness,family = genpois())
  em_treat_richness <- emmeans(glmm_richness, ~ treatment, type = "response")
  em_time_richness <- emmeans(glmm_richness, ~ treatment | sampling, type = "response")
  
  ## ABUNDANCE ##
  # We use gaussian because it tolerates 0 (present in sampling 1 for treatments p and wp)
  glmm_abundance <- glmmTMB(abundance ~ treatment * sampling + ar1(sampling + 0 | plot),
                            dispformula = ~ treatment, # Allows variance/dispersion of treatments to be estimated independently
                            data = arkaute_abundance, family = lognormal(link = "log")
  )
  em_treat_abundance <- emmeans(glmm_abundance, ~ treatment, type = "response")
  em_time_abundance <- emmeans(glmm_abundance, ~ treatment | sampling, type = "response")
  
  ## EVENNESS ##
  glmm_evenness <- glmmTMB( Y_zipf ~ treatment * poly(sampling_num, 3) + ar1(sampling_f + 0 | plot),
                            data = arkaute_evenness, family = t_family(link = "identity"))
  em_treat_evenness <- emmeans(glmm_evenness, ~ treatment, type = "response")
  em_time_evenness <- emmeans(glmm_evenness, ~ treatment | sampling_num,
    at = list(sampling_num = unique(arkaute$sampling_num)),type = "response")
  
  ## SLA ##
  glmm_sla <- glmmTMB(SLA ~ treatment * poly(sampling_num, 3) + ar1(sampling_f + 0 | plot),
                      dispformula = ~ treatment,
                      data = arkaute_sla, family = lognormal(link = "log"))
  em_treat_sla <- emmeans(glmm_sla, ~ treatment, type = "response")
  em_time_sla <- emmeans(glmm_sla, ~ treatment | sampling_num,
                              at = list(sampling_num = unique(arkaute$sampling_num)),type = "response")
  
  ## LDMC ##
  glmm_LDMC <- glmmTMB( LDMC ~ treatment * poly(sampling_num, 3) + ar1(sampling_f + 0 | plot),
                        dispformula = ~treatment + poly(sampling_num, 2), 
                        data = arkaute_ldmc, family = lognormal(link = "log"))
  em_treat_ldmc<- emmeans(glmm_LDMC, ~ treatment, type = "response")
  em_time_ldmc <- emmeans(glmm_LDMC, ~ treatment | sampling_num,
                         at = list(sampling_num = unique(arkaute$sampling_num)),type = "response")
  
  ## Leaf nitrogen ##
  glmm_leafN <-  glmmTMB(leafN ~ treatment * poly(sampling_num, 3) + ar1(sampling_f + 0 | plot),
                         dispformula = ~treatment, 
                         data = arkaute_leafN, family = gaussian(link = "identity"))
  em_treat_leafN <- emmeans(glmm_leafN, ~ treatment, type = "response")
  em_time_leafN <- emmeans(glmm_leafN, ~ treatment | sampling_num,
                          at = list(sampling_num = unique(arkaute$sampling_num)),type = "response")
  
  ### BIOMASS ###
  glmm_biomass <- glmmTMB(biomass_mice_lm ~ treatment * sampling + ar1(sampling + 0 | plot),
                          dispformula = ~treatment,
                          data = arkaute_biomass,  family = Gamma(link = "log"))
  em_treat_biomass <- emmeans(glmm_biomass, ~ treatment, type = "response")
  em_time_biomass <- emmeans(glmm_biomass, ~ treatment | sampling, type = "response")
  ###################################
  
  #### Joining GLMM results ####
  
  glmm_em_treat_list <- list(em_treat_richness, em_treat_abundance, em_treat_evenness, em_treat_sla, 
                        em_treat_ldmc, em_treat_leafN, em_treat_biomass)
  
  glmm_em_time_list <- list(em_time_richness, em_time_abundance, em_time_evenness, em_time_sla, 
                        em_time_ldmc, em_time_leafN, em_time_biomass)
  
  glmm_em_list <- list(glmm_em_treat_list, glmm_em_time_list)
  
  glmm_results_treat_list <- list()
  glmm_results_time_list <- list()
  glmm_results <- list(glmm_results_treat_list, glmm_results_time_list)
  
for(i in seq_along(em_list)){
  
  # 1 is glmm at treatment level
  # 2 is glmm at sampling level
  
  glmm_results[[i]][[1]] <- as.data.frame(pairs(glmm_em_list[[i]][[1]], adjust = "tukey")) |>
    mutate(variable = paste0("richness"), AIC = AIC(glmm_richness), estimate_type = "ratio") |> 
    rename(estimate = ratio)|> 
    select(-null)
  
  glmm_results[[i]][[2]]  <- as.data.frame(pairs(glmm_em_list[[i]][[2]], adjust = "tukey")) |> 
    mutate(variable = paste0("abundance"), AIC = AIC(glmm_abundance), estimate_type = "ratio")|> 
    rename(estimate = ratio)|> 
    select(-null)
  
  glmm_results[[i]][[3]]  <- as.data.frame(pairs(glmm_em_list[[i]][[3]], adjust = "tukey")) |>
    mutate(variable = paste0("evenness"), AIC = AIC(glmm_evenness), estimate_type = "substract")
  
  glmm_results[[i]][[4]]  <- as.data.frame(pairs(glmm_em_list[[i]][[4]], adjust = "tukey")) |> 
    mutate(variable = paste0("SLA"), AIC = AIC(glmm_sla), estimate_type = "ratio")|> 
    rename(estimate = ratio) |> 
    select(-null)
  
  glmm_results[[i]][[5]]  <-  as.data.frame(pairs(glmm_em_list[[i]][[5]], adjust = "tukey")) |>
    mutate(variable = paste0("LDMC"), AIC = AIC(glmm_LDMC), estimate_type = "ratio")|> 
    rename(estimate = ratio)|> 
    select(-null)
  
  glmm_results[[i]][[6]]  <- as.data.frame(pairs(glmm_em_list[[i]][[6]], adjust = "tukey")) |> 
    mutate(variable = paste0("leafN"), AIC = AIC(glmm_leafN), estimate_type = "substract")
  
  glmm_results[[i]][[7]]  <- as.data.frame(pairs(glmm_em_list[[i]][[7]], adjust = "tukey")) |>
    mutate(variable = paste0("biomass"), AIC = AIC(glmm_biomass), estimate_type = "ratio")|> 
    rename(estimate = ratio) |> 
    select(-null)
}
  
  glmm_treatment <- do.call(rbind, glmm_results[[1]]) |> 
    common_function() |> 
    mutate(model = paste0("GLMM")) |>  
    select(-contrast)
  
  # For those models where we used poly(sampling_num, X) we need to change the name of the variable
  glmm_results[[2]][[3]] <- glmm_results[[2]][[3]] |> rename(sampling = sampling_num)
  glmm_results[[2]][[4]] <- glmm_results[[2]][[4]] |> rename(sampling = sampling_num)
  glmm_results[[2]][[5]] <- glmm_results[[2]][[5]] |> rename(sampling = sampling_num)
  glmm_results[[2]][[6]] <- glmm_results[[2]][[6]] |> rename(sampling = sampling_num)
  
  glmm_dynamics <- do.call(rbind, glmm_results[[2]]) |> 
    common_function() |> 
    mutate(model = paste0("GLMM")) |>  
    select(-contrast)
  
  
  
  
  
  ############# GAM MODELS ###########################################################
  samplings_eval <- sort(unique(arkaute_evenness$sampling_num))
  
  ## Richness ## 
  gam_richness <-  gam(
    richness ~ treatment + s(sampling_num, by = treatment, k = 10) + s(plot, bs = "re"),
    data = arkaute_richness, family = gaussian())
  gam_em_treat_richness <- emmeans(gam_richness, ~ treatment)
  gam_em_time_richness  <- emmeans(gam_richness, ~ treatment | sampling_num, at = list(sampling_num = samplings_eval))
  
  ## Abundance ##
  gam_abundance <- gam(abundance ~ treatment + s(sampling_num, by = treatment, k = 18) + 
      s(plot, bs = "re"), data = arkaute_abundance,  family = gaussian(link = "identity"))
  gam_em_treat_abundance <- emmeans(gam_abundance, ~ treatment)
  gam_em_time_abundance  <- emmeans(gam_abundance, ~ treatment | sampling_num, at = list(sampling_num = samplings_eval))
  
  ## Evenness ##
  gam_evenness <- gam( Y_zipf ~ treatment + s(sampling_num, by = treatment, k = 15) + s(plot, bs = "re"),
                       data = arkaute_evenness, family = scat())
  gam_em_treat_evenness <- emmeans(gam_evenness, ~ treatment)
  gam_em_time_evenness  <- emmeans(gam_evenness, ~ treatment | sampling_num, at = list(sampling_num = samplings_eval))
  
  ## SLA ##
  gam_sla <- gam(SLA ~ treatment + s(sampling_num, by = treatment, k = 10) + s(plot, bs = "re"),
                 data = arkaute_sla, family = gaussian(link = "log"))
  gam_em_treat_sla <- emmeans(gam_sla, ~ treatment, type = "response")
  gam_em_time_sla  <- emmeans(gam_sla, ~ treatment | sampling_num, at = list(sampling_num = samplings_eval),
                          type = "response")
  
  ## LDMC ##
  gam_ldmc <- gam(LDMC ~ treatment + s(sampling_num, by = treatment, k = 18) + s(plot, bs = "re"),
                  data = arkaute_ldmc, family = gaussian())
  gam_em_treat_ldmc <- emmeans(gam_ldmc, ~ treatment)
  gam_em_time_ldmc  <- emmeans(gam_ldmc, ~ treatment | sampling_num, at = list(sampling_num = samplings_eval))
  
  ## Leaf Nitrogen ##
  gam_leafN <- gam(leafN ~ treatment + s(sampling_num, by = treatment, k = 10) + s(plot, bs = "re"),
                   data = arkaute_leafN, family = gaussian(link = "identity"))
  gam_em_treat_leafN <- emmeans(gam_leafN, ~ treatment)
  gam_em_time_leafN  <- emmeans(gam_leafN, ~ treatment | sampling_num, at = list(sampling_num = samplings_eval))
  
  ## Biomass ##
  gam_biomass <- gam(biomass_mice_lm ~ treatment + s(sampling_num, by = treatment, k = 15) + s(plot, bs = "re"),
                     data = arkaute_biomass, family = tw(link = "log"))
  gam_em_treat_biomass <-  emmeans(gam_biomass, ~ treatment, type = "response")
  gam_em_time_biomass  <- emmeans(gam_biomass, ~ treatment | sampling_num, at = list(sampling_num = samplings_eval),
                              type = "response")
  
  ###############
  
  gam_em_treat_list <- list(gam_em_treat_richness, gam_em_treat_abundance, gam_em_treat_evenness, gam_em_treat_sla, 
                             gam_em_treat_ldmc, gam_em_treat_leafN, gam_em_treat_biomass)
  
  gam_em_time_list <- list(gam_em_time_richness, gam_em_time_abundance, gam_em_time_evenness, gam_em_time_sla, 
                            gam_em_time_ldmc, gam_em_time_leafN, gam_em_time_biomass)
  
  gam_em_list <- list(gam_em_treat_list, gam_em_time_list)
  
  gam_results_treat_list <- list()
  gam_results_time_list <- list()
  gam_results <- list(gam_results_treat_list, gam_results_time_list)
  
  
  
  for(i in seq_along(gam_em_list)){
  
  
  gam_results[[i]][[1]] <- as.data.frame(pairs(gam_em_list[[i]][[1]], adjust = "tukey")) |>
    mutate(variable = "richness", AIC = AIC(gam_richness), estimate_type = "substract")
  
  gam_results[[i]][[2]] <- as.data.frame(pairs(gam_em_list[[i]][[2]], adjust = "tukey")) |>
    mutate(variable = "abundance", AIC = AIC(gam_abundance), estimate_type = "substract")
  
  gam_results[[i]][[3]] <- as.data.frame(pairs(gam_em_list[[i]][[3]], adjust = "tukey")) |>
    mutate(variable = "evenness", AIC = AIC(gam_evenness), estimate_type = "substract")
  
  gam_results[[i]][[4]] <- as.data.frame(pairs(gam_em_list[[i]][[4]], adjust = "tukey")) |>
    mutate(variable = "SLA", AIC = AIC(gam_sla), estimate_type = "ratio") |>
    rename(estimate = ratio) |> 
    select(-null)
  
  gam_results[[i]][[5]] <- as.data.frame(pairs(gam_em_list[[i]][[5]], adjust = "tukey")) |>
    mutate(variable = "LDMC", AIC = AIC(gam_ldmc), estimate_type = "substract")
  
  gam_results[[i]][[6]] <- as.data.frame(pairs(gam_em_list[[i]][[6]], adjust = "tukey")) |>
    mutate(variable = "leafN", AIC = AIC(gam_leafN), estimate_type = "substract")
  
  gam_results[[i]][[7]] <- as.data.frame(pairs(gam_em_list[[i]][[7]], adjust = "tukey")) |>
    mutate(variable = "biomass", AIC = AIC(gam_biomass), estimate_type = "ratio") |>
    rename(estimate = ratio) |> 
    select(-null)
  }
  
  gam_treatment <- do.call(rbind,  gam_results[[1]]) |>
    common_function() |> 
    mutate(model = "GAM") |>
    select(-contrast)
  

  gam_dynamics <- do.call(rbind, gam_results[[2]]) |>
    common_function() |> 
    mutate(model = "GAM") |>
    rename(sampling = sampling_num) |> 
    select(-contrast)
  
  
  
  ####### Joining GLMM, GAM and LRR ############
  
  
  model_agg <- full_join(glmm_treatment, gam_treatment) |> 
    full_join(lrr_table) |> 
    mutate(variable.bis = variable) |> 
    select(variable, eff_descriptor, model, AIC, p_value, effect_significance, estimate,
           effect_sign, variable.bis)
  
  
  models_dynamics <- full_join(glmm_dynamics, gam_dynamics) |>
    full_join(lrr_table_dyn) |> 
    mutate(
      variable.bis = variable,
      variable_model = paste0(variable, "-", model),
      sampling = fct_reorder(sampling, as.numeric(sampling))
    ) |>
    select(
      variable, sampling, eff_descriptor, model,p_value, effect_significance, estimate,
      effect_sign, variable.bis, variable_model)
  
  
#}


### VISUALIZATION ###


gg_wpp <- 
  model_agg |> 
  filter(eff_descriptor == "wp_vs_p") |> 
  ggplot(aes(x = model, y = variable, color = effect_significance, shape = effect_sign,
             label = round(p_value, 2))) +
  geom_point(size = 10) +
  geom_label_repel(show.legend = FALSE) +
  scale_color_manual(values = palette_significance) +
  scale_shape_manual(values = shape_significance) +
  labs(title = "Warming effect on recovery (wp vs p)",
       subtitle = "Model: variable ~ treatment * sampling + ar(sampling) + (1 | plot)",
       x = "Statistical analysis", y = NULL,
       shape = "Sign of effect", color = "Effect stat. significance")
print(gg_wpp)

gg_pc <- 
  model_agg |> 
  filter(eff_descriptor == "p_vs_c") |> 
  ggplot(aes(x = model, y = variable, color = effect_significance, shape = effect_sign,
             label = round(p_value, 2))) +
  geom_point(size = 10) +
  geom_label_repel(show.legend = FALSE) +
  scale_color_manual(values = palette_significance) +
  scale_shape_manual(values = shape_significance) +
  labs(title = "Perturbation effect (p vs c)",
       subtitle = "Model: variable ~ treatment * sampling + ar(sampling) + (1 | plot)",
       x = "Statistical analysis", y = NULL,
       shape = "Sign of effect", color = "Effect stat. significance")
print(gg_pc)


gg_wc <- 
  model_agg |> 
  filter(eff_descriptor == "w_vs_c") |> 
  ggplot(aes(x = model, y = variable, color = effect_significance, shape = effect_sign,
             label = round(p_value, 2))) +
  geom_point(size = 10) +
  geom_label_repel(show.legend = FALSE) +
  scale_color_manual(values = palette_significance) +
  scale_shape_manual(values = shape_significance) +
  labs(title = "Warming effect on assembly (w vs c)",,
       subtitle = "Model: variable ~ treatment * sampling + ar(sampling) + (1 | plot)",
       x = "Statistical analysis", y = NULL,
       shape = "Sign of effect", color = "Effect stat. significance")
print(gg_wc)

gg_wpc <-
  model_agg |> 
  filter(eff_descriptor == "wp_vs_c") |> 
  ggplot(aes(x = model, y = variable, color = effect_significance, shape = effect_sign,
             label = round(p_value, 2))) +
  geom_point(size = 10) +
  geom_label_repel(show.legend = FALSE) +
  scale_color_manual(values = palette_significance) +
  scale_shape_manual(values = shape_significance) +
  labs(title = "Combined effect (wp vs c)",
       subtitle = "Model: variable ~ treatment * sampling + ar(sampling) + (1 | plot)",
       x = "Statistical analysis", y = NULL,
       shape = "Sign of effect", color = "Effect stat. significance")
print(gg_wpc)


gg_wpp_dyn <- 
  models_dynamics |> 
  filter(eff_descriptor == "wp_vs_p") |> 
  ggplot(aes(x = sampling, y = variable_model, color = effect_significance, shape = effect_sign)) +
  geom_point(size = 8) +
  geom_hline(yintercept = seq(3.5, 18.5, by = 3), linetype = "solid", linewidth = 0.5) +
  scale_color_manual(values = palette_significance) +
  scale_shape_manual(values = shape_significance) +
  labs(title = "Warming effect on recovery (wp/p)",
       subtitle = "Model: variable ~ treatment * sampling + ar(sampling) + (1 | plot)",
       x = "Sampling", 
       y = "Variable and model",
       color = "Effect significance",
       shape = "Sign of the effect")
print(gg_wpp_dyn)




gg_pc_dyn <- 
  models_dynamics |> 
  filter(eff_descriptor == "p_vs_c") |> 
  ggplot(aes(x = sampling, y = variable_model, color = effect_significance, shape = effect_sign)) +
  geom_point(size = 8) +
  geom_hline(yintercept = seq(3.5, 18.5, by = 3), linetype = "solid", linewidth = 0.5) +
  scale_color_manual(values = palette_significance) +
  scale_shape_manual(values = shape_significance) + 
  labs(title = "Recovery (p/c)",
       subtitle = "Model: variable ~ treatment * sampling + ar(sampling) + (1 | plot)",
       x = "Sampling", 
       y = "Variable and model")
print(gg_pc_dyn)



gg_wc_dyn <- 
  models_dynamics |> 
  filter(eff_descriptor == "w_vs_c") |> 
  ggplot(aes(x = sampling, y = variable_model, color = effect_significance, shape = effect_sign)) +
  geom_point(size = 8) +
  geom_hline(yintercept = seq(3.5, 18.5, by = 3), linetype = "solid", linewidth = 0.5) +
  scale_color_manual(values = palette_significance) +
  scale_shape_manual(values = shape_significance) +
  labs(title = "Warming effect(w/c)",
       subtitle = "Model: variable ~ treatment * sampling + ar(sampling) + (1 | plot)",
       x = "Sampling", 
       y = "Variable and model")
print(gg_wc_dyn)




gg_wpc_dyn <- 
  models_dynamics |> 
  filter(eff_descriptor == "wp_vs_c") |> 
  ggplot(aes(x = sampling, y = variable_model, color = effect_significance, shape = effect_sign)) +
  geom_point(size = 8) +
  geom_hline(yintercept = seq(3.5, 18.5, by = 3), linetype = "solid", linewidth = 0.5) +
  scale_color_manual(values = palette_significance) +
  scale_shape_manual(values = shape_significance) +
  labs(title = "Combined effect(wp/c)", 
       subtitle = "Model: variable ~ treatment * sampling + ar(sampling) + (1 | plot)",
       x = "Sampling", 
       y = "Variable and model")
print(gg_wpc_dyn)




ggsave("results/model_comparison_time_treatment_wp_vs_p.png", plot = gg_wpp, dpi = 600)
ggsave("results/model_comparison_time_treatment_p_vs_c.png", plot = gg_pc, dpi = 600)
ggsave("results/model_comparison_time_treatment_w_vs_c.png", plot = gg_wc, dpi = 600)
ggsave("results/model_comparison_time_treatment_wp_vs_c.png", plot = gg_wpc, dpi = 600)

ggsave("results/model_comparison_time_sampling_wp_vs_p.png", plot = gg_wpp_dyn, dpi = 600)
ggsave("results/model_comparison_time_sampling_w_vs_c.png", plot = gg_pc_dyn, dpi = 600)
ggsave("results/model_comparison_time_sampling_p_vs_c.png", plot = gg_wc_dyn, dpi = 600)
ggsave("results/model_comparison_time_sampling_wp_vs_c.png", plot = gg_wpc_dyn, dpi = 600)
