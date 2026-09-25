


rm(list = ls(all.names = TRUE))
pacman::p_load(
  dplyr, reshape2, tidyverse, lubridate, ggplot2, ggpubr, gridExtra, ggrepel,
  car, ggsignif, dunn.test, rstatix, performance, emmeans, DHARMa, mgcv, glmmTMB
)

source("code/palettes_labels.R")


{
  
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
  arkaute_biomass_mice <- arkaute |> filter(!is.na(biomass_mice))
  arkaute_biomass_raw <- arkaute |> filter(!is.na(biomass_raw), biomass_raw != 0 )
  
  
  
  
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
  
  
  
  
  extract_glmm_summary <- function(model, var_name) {
    vc <- VarCorr(model)
    sd_plot <- attr(vc$cond$plot, "stddev")["(Intercept)"]

    disp_val <- suppressWarnings(sigma(model))
    if (is.na(disp_val)) {
      disp_val <- mean(predict(model, type = "disp"))
    }
    
    glmm_summary_result <<- 
      as.data.frame(car::Anova(model, type = "III")) |> 
      tibble::rownames_to_column(var = "Effect") |> 
      mutate(
        SD_re_plot           = as.numeric(sd_plot),
        Var_re_plot          = as.numeric(sd_plot^2),
        dispersion_parameter = disp_val,
        AIC                  = AIC(model),
        Family               = family(model)$family,
        Link                 = family(model)$link,
        Variable             = var_name
      )
    
    return(glmm_summary_result)
    
  }

  
  ############# GLMM MODELS ##############################################################
  ## Richness ## 
  glmm_richness <- glmmTMB( richness ~ treatment * sampling + (1 | plot),
                            data = arkaute_richness,family = genpois())
  em_treat_richness <- emmeans(glmm_richness, ~ treatment, type = "response")
  em_time_richness <- emmeans(glmm_richness, ~ treatment | sampling, type = "response")
  
  ## ABUNDANCE ##
  # We use gaussian because it tolerates 0 (present in sampling 1 for treatments p and wp)
  glmm_abundance <- glmmTMB(abundance ~ treatment * sampling + (1 | plot),
                            data = arkaute_abundance, family = gaussian())
  em_treat_abundance <- emmeans(glmm_abundance, ~ treatment, type = "response")
  em_time_abundance <- emmeans(glmm_abundance, ~ treatment | sampling, type = "response")
  
  ## EVENNESS ##
  glmm_evenness <- glmmTMB(Y_zipf ~ treatment * sampling + (1 | plot),
                           dispformula = ~treatment + sampling, 
                           data = arkaute_evenness, family = gaussian())
  em_treat_evenness <- emmeans(glmm_evenness, ~ treatment, type = "response")
  em_time_evenness <- emmeans(glmm_evenness, ~ treatment | sampling, type = "response")
  
  ## SLA ##
  glmm_sla <- glmmTMB(SLA ~ treatment * sampling + (1 | plot), 
                      dispformula = ~ treatment, 
                      data = arkaute_sla, family = gaussian())
  em_treat_sla <- emmeans(glmm_sla, ~ treatment, type = "response")
  em_time_sla <- emmeans(glmm_sla, ~ treatment | sampling, type = "response")
  
  ## LDMC ##
  glmm_LDMC <- glmmTMB(LDMC ~ treatment * sampling + (1 | plot), 
                       dispformula = ~ treatment, 
                       data = arkaute_ldmc, family = gaussian())
  em_treat_ldmc<- emmeans(glmm_LDMC, ~ treatment, type = "response")
  em_time_ldmc <- emmeans(glmm_LDMC, ~ treatment | sampling, type = "response")
  
  ## Leaf nitrogen ##
  glmm_leafN <-  glmmTMB(leafN ~ treatment * sampling + (1 | plot),
                         dispformula = ~ treatment,
                         data = arkaute_leafN, family = gaussian())
  em_treat_leafN <- emmeans(glmm_leafN, ~ treatment, type = "response")
  em_time_leafN <- emmeans(glmm_leafN, ~ treatment | sampling, type = "response")
  
  ### BIOMASS ###
  glmm_biomass <- glmmTMB(biomass_mice_lm ~ treatment * sampling + mean_vwc + (1 | plot),
                          dispformula = ~treatment + sampling, 
                          data = arkaute_biomass, family = Gamma(link = "log"))
  em_treat_biomass <- emmeans(glmm_biomass, ~ treatment, type = "response")
  em_time_biomass <- emmeans(glmm_biomass, ~ treatment | sampling, type = "response")
  
  
  ## BIOMASS - just MICE ####
  
  glmm_biomass_mice <- glmmTMB(biomass_mice ~ treatment * sampling + mean_vwc +  (1 | plot),
                                dispformula = ~treatment + sampling, 
                                data = arkaute_biomass_mice,  family = Gamma(link = "log"))
  em_treat_biomass_mice <- emmeans(glmm_biomass_mice, ~ treatment, type = "response")
  em_time_biomass_mice <- emmeans(glmm_biomass_mice, ~ treatment | sampling, type = "response")
  
  
  ## BIOMASS raw ##
  
  glmm_biomass_raw <- glmmTMB(biomass_raw ~ treatment * sampling + (1 | plot),
                              dispformula = ~treatment + sampling, 
                              data = arkaute_biomass_raw,  family = Gamma(link = "log"))
  em_treat_biomass_raw <- emmeans(glmm_biomass_raw, ~ treatment, type = "response")
  em_time_biomass_raw <- emmeans(glmm_biomass_raw, ~ treatment | sampling, type = "response")
  
  
  
  ################################### GLMM RESULTS ########################################
  
  
  
  ### Joining GLMM summary results #############
  
  
  
  glmm_list <- list(glmm_richness, glmm_abundance, glmm_evenness, glmm_sla,
                    glmm_LDMC, glmm_leafN, glmm_biomass, glmm_biomass_mice, glmm_biomass_raw)
  variables <- c("richness", "abundance", "Y_zipf", "SLA", "LDMC", "leafN", "biomass_mice_lm", 
                 "biomass_mice", "biomass_raw")
  result_list <- list()
  
  

  
  for(i in seq_along(glmm_list)){
    extract_glmm_summary(glmm_list[[i]], variables[i])
    result_list[[i]] <- glmm_summary_result
  }
  
  glmm_summary_result <- do.call(rbind, result_list)
  
 
  
  
  #### Joining GLMM post hoc comparison results ####
  
  glmm_em_treat_list <- list(em_treat_richness, em_treat_abundance, em_treat_evenness, em_treat_sla, 
                             em_treat_ldmc, em_treat_leafN, em_treat_biomass, em_treat_biomass_mice, em_treat_biomass_raw)
  
  glmm_em_time_list <- list(em_time_richness, em_time_abundance, em_time_evenness, em_time_sla, 
                            em_time_ldmc, em_time_leafN, em_time_biomass, em_time_biomass_mice, em_time_biomass_raw)
  
  glmm_em_list <- list(glmm_em_treat_list, glmm_em_time_list)
  
  glmm_results_treat_list <- list()
  glmm_results_time_list <- list()
  glmm_results <- list(glmm_results_treat_list, glmm_results_time_list)
  
  for(i in seq_along(glmm_em_list)){
    

    glmm_results[[i]][[1]] <- as.data.frame(pairs(glmm_em_list[[i]][[1]], adjust = "tukey")) |>
      mutate(variable = paste0("richness"), AIC = AIC(glmm_richness), estimate_type = "ratio") |> 
      rename(estimate = ratio)|> 
      select(-null)
    
    glmm_results[[i]][[2]]  <- as.data.frame(pairs(glmm_em_list[[i]][[2]], adjust = "tukey")) |> 
      mutate(variable = paste0("abundance"), AIC = AIC(glmm_abundance), estimate_type = "substract")

    
    glmm_results[[i]][[3]]  <- as.data.frame(pairs(glmm_em_list[[i]][[3]], adjust = "tukey")) |>
      mutate(variable = paste0("evenness"), AIC = AIC(glmm_evenness), estimate_type = "substract")
    
    glmm_results[[i]][[4]]  <- as.data.frame(pairs(glmm_em_list[[i]][[4]], adjust = "tukey")) |> 
      mutate(variable = paste0("SLA"), AIC = AIC(glmm_sla), estimate_type = "substract")
    
    glmm_results[[i]][[5]]  <-  as.data.frame(pairs(glmm_em_list[[i]][[5]], adjust = "tukey")) |>
      mutate(variable = paste0("LDMC"), AIC = AIC(glmm_LDMC), estimate_type = "substract")
    
    glmm_results[[i]][[6]]  <- as.data.frame(pairs(glmm_em_list[[i]][[6]], adjust = "tukey")) |> 
      mutate(variable = paste0("leafN"), AIC = AIC(glmm_leafN), estimate_type = "substract")
    
    
    glmm_results[[i]][[7]]  <- as.data.frame(pairs(glmm_em_list[[i]][[7]], adjust = "tukey")) |>
      mutate(variable = paste0("biomass"), AIC = AIC(glmm_biomass), estimate_type = "ratio") |> 
      rename(estimate = ratio)|> 
      select(-null)
    
    
    glmm_results[[i]][[8]]  <- as.data.frame(pairs(glmm_em_list[[i]][[8]], adjust = "tukey")) |> 
      mutate(variable = paste0("biomass_mice"), AIC = AIC(glmm_biomass_mice), estimate_type = "ratio") |> 
      rename(estimate = ratio)|> 
      select(-null)
    
    glmm_results[[i]][[9]]  <- as.data.frame(pairs(glmm_em_list[[i]][[9]], adjust = "tukey")) |> 
      mutate(variable = paste0("biomass_raw"), AIC = AIC(glmm_biomass_raw), estimate_type = "ratio") |> 
      rename(estimate = ratio)|> 
      select(-null)
    
  }
  
  glmm_treatment <- do.call(rbind, glmm_results[[1]]) |> 
    common_function() |> 
    mutate(model = paste0("GLMM"))
  
  glmm_dynamics <- do.call(rbind, glmm_results[[2]]) |> 
    common_function() |> 
    mutate(model = paste0("GLMM")) 
  
  
  
  glmm_treatment |>  write.csv("results/GLMM_agg.csv")
  glmm_dynamics  |>  write.csv("results/GLMM_dyn.csv")
  
  glmm_summary_result |> write.csv("results/GLMM_summary_results.csv")
  
}


