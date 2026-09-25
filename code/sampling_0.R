



rm(list = ls(all.names = TRUE))
pacman::p_load(
  dplyr, reshape2, tidyverse, lubridate, ggplot2, ggpubr, gridExtra, ggrepel,
  car, ggsignif, dunn.test, rstatix, performance, emmeans, DHARMa, mgcv, glmmTMB
)

source("code/palettes_labels.R")


{
  

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
    filter(
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
  
  diagnose_glmm <- function(model, data = NULL, group_var = NULL) {
    
    print(summary(model))
    
    # Simulate DHARMa residuals
    sim <- DHARMa::simulateResiduals(fittedModel = model)
    plot(sim)
    
    # --- DHARMa DIAGNOSTIC TESTS ---
    cat("       DHARMa DIAGNOSTIC TESTS          \n")
    cat("\n1. Uniformity Test (KS Test):\n")
    cat("   -> Evaluates overall distributional fit. A significant p-value (p < 0.05) indicates model mis-specification.\n")
    print(DHARMa::testUniformity(sim))
    
    cat("\n2. Dispersion Test:\n")
    cat("   -> Evaluates residual variance. Detects if data suffer from overdispersion (ratio > 1) or underdispersion (ratio < 1).\n")
    print(DHARMa::testDispersion(sim))
    
    cat("\n3. Outliers Test:\n")
    cat("   -> Evaluates extreme values. Tests whether the frequency of 0 or 1 scaled residuals exceeds expectation.\n")
    print(DHARMa::testOutliers(sim))
    
    cat("\n4. Quantiles Test:\n")
    cat("   -> Evaluates residual patterns across predictions. Detects non-linearities or heteroscedasticity along fitted values.\n")
    print(DHARMa::testQuantiles(sim))
    
    # --- GROUP-LEVEL HETEROSCEDASTICITY TEST ---
    if (!is.null(data) && !is.null(group_var)) {
      grp  <- data[[group_var]]
      resu <- sim$scaledResiduals
      cat("\n----------------------------------------\n")
      cat("5. Levene Test on DHARMa residuals by", group_var, ":\n")
      cat("   -> Evaluates variance homogeneity. Tests if residual spread differs significantly across levels of the grouping factor.\n")
      print(car::leveneTest(resu ~ grp))
      DHARMa::plotResiduals(sim, form = grp)
    }
    
    invisible(sim)
  }
  
  
  
  lrr_reading <- function(data) {
    data 
  }
  
  
  
  ## Opening LOG RESPONSE RATIO results datasets


  
  ############# GLMM MODELS ##############################################################
  # I diagnose the models again because I am including the sampling 0 
  
  ## Richness ## 
  glmm_richness <- glmmTMB( richness ~ treatment * sampling + (1 | plot),
                            data = arkaute_richness,family = genpois())
  #diagnose_glmm(glmm_richness)
  em_time_richness <- emmeans(glmm_richness, ~ treatment | sampling, type = "response")
  
  ## ABUNDANCE ##
  # We use gaussian because it tolerates 0 (present in sampling 1 for treatments p and wp)
  glmm_abundance <- glmmTMB(abundance ~ treatment * sampling + (1 | plot),
                            data = arkaute_abundance, family = gaussian())
  #diagnose_glmm(glmm_abundance)
  em_time_abundance <- emmeans(glmm_abundance, ~ treatment | sampling, type = "response")
  
  ## EVENNESS ##
  glmm_evenness <- glmmTMB(Y_zipf ~ treatment * sampling + (1 | plot),
                           dispformula = ~treatment + sampling, 
                           data = arkaute_evenness, family = gaussian())
  #diagnose_glmm(glmm_evenness)
  em_time_evenness <- emmeans(glmm_evenness, ~ treatment | sampling, type = "response")
  
  ## SLA ##
  glmm_sla <- glmmTMB(SLA ~ treatment * sampling +  (1 | plot), 
                      dispformula = ~ treatment + sampling, 
                      data = arkaute_sla, family = gaussian())
  #AIC(glmm_sla)
  #diagnose_glmm(glmm_sla) # There is a very slight overdispersion
  em_time_sla <- emmeans(glmm_sla, ~ treatment | sampling, type = "response")
  
  ## LDMC ##
  glmm_LDMC <- glmmTMB(LDMC ~ treatment * sampling + (1 | plot), 
                       dispformula = ~ treatment, 
                       data = arkaute_ldmc, family = gaussian())
  #diagnose_glmm(glmm_LDMC)
  em_time_ldmc <- emmeans(glmm_LDMC, ~ treatment | sampling, type = "response")
  
  ## Leaf nitrogen ##
  glmm_leafN <-  glmmTMB(leafN ~ treatment * sampling + (1 | plot),
                         dispformula = ~ treatment + sampling,
                         data = arkaute_leafN, family = gaussian())
  #diagnose_glmm(glmm_leafN)
  em_time_leafN <- emmeans(glmm_leafN, ~ treatment | sampling, type = "response")
  
  ### BIOMASS ###
  glmm_biomass <- glmmTMB(biomass_mice_lm ~ treatment * sampling + mean_vwc + (1 | plot),
                          dispformula = ~treatment + sampling, 
                          data = arkaute_biomass, family = Gamma(link = "log"))
  #diagnose_glmm(glmm_biomass)
  em_time_biomass <- emmeans(glmm_biomass, ~ treatment | sampling, type = "response")
  
  
  ## BIOMASS - just MICE ####
  
  glmm_biomass_mice <- glmmTMB(biomass_mice ~ treatment * sampling + mean_vwc +  (1 | plot),
                               dispformula = ~treatment + sampling, 
                               data = arkaute_biomass_mice,  family = Gamma(link = "log"))
  #diagnose_glmm(glmm_biomass_mice)
  em_time_biomass_mice <- emmeans(glmm_biomass_mice, ~ treatment | sampling, type = "response")
  
  
  ## BIOMASS raw ##
  
  glmm_biomass_raw <- glmmTMB(biomass_raw ~ treatment * sampling + (1 | plot),
                              dispformula = ~treatment + sampling, 
                              data = arkaute_biomass_raw,  family = Gamma(link = "log"))
  #diagnose_glmm(glmm_biomass_raw)
  em_time_biomass_raw <- emmeans(glmm_biomass_raw, ~ treatment | sampling, type = "response")
  
  
  
  


  
  glmm0_list <- list()
  
  glmm0_list[[1]] <- as.data.frame(pairs(em_time_richness, adjust = "tukey")) |>
    mutate(variable = paste0("richness"), AIC = AIC(glmm_richness), estimate_type = "ratio") |> 
    rename(estimate = ratio)|> 
    select(-null)
  
  glmm0_list[[2]]  <- as.data.frame(pairs(em_time_abundance, adjust = "tukey")) |> 
    mutate(variable = paste0("abundance"), AIC = AIC(glmm_abundance), estimate_type = "substract")
  
  glmm0_list[[3]]  <- as.data.frame(pairs(em_time_evenness, adjust = "tukey")) |>
    mutate(variable = paste0("evenness"), AIC = AIC(glmm_evenness), estimate_type = "substract")
  
  glmm0_list[[4]]  <- as.data.frame(pairs(em_time_sla, adjust = "tukey")) |> 
    mutate(variable = paste0("SLA"), AIC = AIC(glmm_sla), estimate_type = "substract")
  
  glmm0_list[[5]]  <-  as.data.frame(pairs(em_time_ldmc, adjust = "tukey")) |>
    mutate(variable = paste0("LDMC"), AIC = AIC(glmm_LDMC), estimate_type = "substract")
  
  glmm0_list[[6]]  <- as.data.frame(pairs(em_time_leafN, adjust = "tukey")) |> 
    mutate(variable = paste0("leafN"), AIC = AIC(glmm_leafN), estimate_type = "substract")
  
  glmm0_list[[7]]  <- as.data.frame(pairs(em_time_biomass, adjust = "tukey")) |>
    mutate(variable = paste0("biomass"), AIC = AIC(glmm_biomass), estimate_type = "ratio") |> 
    rename(estimate = ratio)|> 
    select(-null)
  
  glmm0_list[[8]]  <- as.data.frame(pairs(em_time_biomass_mice, adjust = "tukey")) |> 
    mutate(variable = paste0("biomass_mice"), AIC = AIC(glmm_biomass_mice), estimate_type = "ratio") |> 
    rename(estimate = ratio)|> 
    select(-null)
  
  glmm0_list[[9]]  <- as.data.frame(pairs(em_time_biomass_raw, adjust = "tukey")) |> 
    mutate(variable = paste0("biomass_raw"), AIC = AIC(glmm_biomass_raw), estimate_type = "ratio") |> 
    rename(estimate = ratio)|> 
    select(-null)

  glmm_results0 <- do.call(rbind, glmm0_list) |> 
    filter(sampling == "0")|> 
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
        contrast %in% c("p - wp", "p / wp")  ~ "wp_vs_p",
      )) |> 
    mutate(
      variable = fct_recode(variable,
                            "Y_zipf" = "evenness",
                            "biomass_mice_lm"  = "biomass"),
      sampling = as.factor(sampling))
  
  
  
  lrr_sampling_0 <- read.csv("results/effect_size_dynamics.csv") |> 
    filter(sampling == "0",
           eff_descriptor != "wp_vs_w") |> 
    mutate(sampling = as.factor(sampling)) |> 
    select(-scale, -X) |> 
    mutate(
      variable = as.factor(variable),
      variable = droplevels(variable),
    ) |> 
    select(-null_effect)
  

  
  RESULT0 <- full_join(lrr_sampling_0, glmm_results0) |> 
    filter(eff_descriptor != "wp_vs_p") |> 
    select(sampling, year, date, date_label_noyear,  eff_descriptor, variable, eff_value,
           lower_limit, upper_limit, contrast, estimate, SE, z.ratio, p_value, AIC, estimate_type,
           effect_sign, effect_significance) |> 
    rename(glmm_contrast = contrast,
           glmm_estimate = estimate, 
           glmm_SE = SE, 
           glmm_z.ratio = z.ratio,
           glmm_p_value = p_value, 
           glmm_AIC = AIC, 
           glmm_estimatetype = estimate_type, 
           glmm_effect_sign = effect_sign,
           glmm_effect_significance = effect_significance)

  
  

   
  


  labs_variable <- c(
    "richness"        = "Richness",            
    "abundance"       = "Cover",                   
    "Y_zipf"          = "Evenness",                   
    "SLA"             = "SLA",                           
    "LDMC"            = "LDMC",                         
    "leafN"           = "LN",                           
    "biomass_mice_lm" = "Biomass"
  )
  
  RESULT0 <- RESULT0 |> 
    mutate(
      eff_descriptor = factor(eff_descriptor),
      x_jit          = as.integer(eff_descriptor) - 2,
      variable       = factor(
        variable, 
        levels = names(labs_variable), 
        labels = unname(labs_variable)
      )
    )

  
  gg_s0 <- 
  ggplot(RESULT0, aes(
    x = eff_descriptor,                 # centrado en 0 + pequeño desplazamiento
    y = eff_value,
    color = eff_descriptor
  )) +
    #facet_grid(rows = vars(variable), scales = "free_y", switch = "y") +
    facet_wrap(~ variable, scales = "free_y", strip.position = "top", nrow = 2) +

    
    geom_hline(yintercept = 0, linetype = "dashed",
               color = "grey50", linewidth = 0.5) +
    
    
    geom_linerange(aes(ymin = lower_limit, ymax = upper_limit),
                   linewidth = 1, alpha = 1) +
    
    geom_point(size = 6) +
    
    geom_text(
      aes(
        y = ifelse(eff_value < 0, lower_limit - 0.1, upper_limit + 0.2),
        label = case_when(
          glmm_effect_significance == "significant"     ~ "*",
          glmm_effect_significance == "marginal"        ~ "·",
          glmm_effect_significance == "non-significant" ~ NA_character_
        )
      ),
      
      #vjust = 0.7,          # Ajuste vertical para centrar el '*' dentro de la figura
      show.legend = FALSE,
      size = 10
    ) +
    
    scale_color_manual(values = palette_RR_CB, labels = c("w_vs_c" = "Warmed-only", "p_vs_c" = "Perturbed-only", "wp_vs_c" = "Combined")) +
    scale_y_continuous( breaks = scales::pretty_breaks(n = 3), 
                        expand = expansion(mult = c(0.05, 0.1)) ) +
    
    
    labs(x = NULL, y = NULL, color = NULL) +
    
    #gg_RR_theme +
    theme(
      strip.background   = element_blank(),
      strip.placement    = "outside",
      strip.text         = element_text(face = "bold", size = 12),
      strip.text.y.left  = element_text(face = "bold"),
      axis.text.y        = element_text(angle = 90, hjust = 0.5, face = "plain", size = 12),
      axis.text.x        = element_blank(),
      axis.ticks.x       = element_blank(),
      legend.position    = "bottom",
      legend.text        = element_text(size = 14, face = "plain")
    )
  
  print(gg_s0)
  
}


ggsave("results/SAMPLING_0.png", plot = gg_s0, dpi = 600)
#ggsave("results/SAMPLING_=.svg", plot = gg_s0, dpi = 600)  
  