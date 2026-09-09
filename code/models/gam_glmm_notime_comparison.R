



rm(list = ls(all.names = TRUE))

# Cargar paquetes
pacman::p_load(
  dplyr, reshape2, tidyverse, lubridate, ggplot2, ggpubr, gridExtra, ggrepel,
  car, ggsignif, dunn.test, rstatix, performance, emmeans, DHARMa, mgcv, glmmTMB
)

# Custom palettes and labels
source("code/palettes_labels.R")


# Global ggplot2 theme
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

# Carga y limpieza de datos
arkaute <- read.csv("data/processed_data/arkaute.csv") %>%
  mutate(
    year      = factor(year),
    date      = ymd(date),
    sampling = factor(sampling, levels = as.character(sort(unique(as.numeric(as.character(sampling)))))),
    plot      = factor(plot),
    treatment = factor(treatment)
  ) %>%
  filter(sampling != "0") |>
  filter(sampling != "1") 
  #filter(!(sampling == "1" & treatment %in% c("p", "wp")))

{
# databases for models:
arkaute_richness <- arkaute |>  filter(!is.na(richness))
arkaute_abundance <- arkaute |>  filter(!is.na(abundance))
arkaute_evenness <- arkaute |> filter(!is.na(Y_zipf))
arkaute_sla <- arkaute |> filter(!is.na(SLA))
arkaute_ldmc <- arkaute |> filter(!is.na(LDMC))
arkaute_leafN <- arkaute |> filter(!is.na(leafN))
arkaute_biomass <- arkaute |> filter(!is.na(biomass_mice_lm))

# Diagnosis functions




  ############# GLMM ##############
  ## Richness ## 
  glmm_richness <- glmmTMB(richness ~ treatment + (1 | plot), data = arkaute, family = gaussian)
  em_treat_richness <- emmeans(glmm_richness, ~ treatment, type = "response")
  
  ##### ABUNDANCE #####
  # The model for abundance present difficulties to be adjusted without problems of residuals dispersion
  glmm_abundance <- glmmTMB(abundance ~ treatment + (1 | plot),data = arkaute, family = gaussian(link = "identity"))
  em_treat_abundance <- emmeans(glmm_abundance, ~ treatment, type = "response")
  
  ##### EVENNESS #####
  glmm_evenness <- glmmTMB(Y_zipf ~ treatment + (1 | plot),
                              dispformula = ~ treatment, 
                              data = arkaute_evenness, family = gaussian())
  em_treat_evenness <- emmeans(glmm_evenness, ~ treatment, type = "response")
  
  ### SLA ###
  glmm_sla <- glmmTMB(SLA ~ treatment + (1 | plot), 
                      dispformula = ~ treatment, 
                      data = arkaute_sla, family = gaussian(link = "log"))
  em_treat_sla <- emmeans(glmm_sla, ~ treatment, type = "response")

  ### LDMC ###
  glmm_LDMC <- glmmTMB(LDMC ~ treatment + (1 | plot), 
                          dispformula = ~ treatment, 
                          data = arkaute_ldmc, family = gaussian())
  em_treat_ldmc<- emmeans(glmm_LDMC, ~ treatment, type = "response")

  ### Leaf nitrogen ###
  glmm_leafN <- glmmTMB(leafN ~ treatment + (1 | plot),
                        dispformula = ~ treatment, 
                        data = arkaute_leafN, family = gaussian(link = "identity"))
  em_treat_leafN <- emmeans(glmm_leafN, ~ treatment, type = "response")

  ### BIOMASS ###
  #glmm_biomass <- glmmTMB(biomass_mice_lm ~ treatment + (1 | plot), data = arkaute_biomass, family = Gamma(link = "log"))
  glmm_biomass <- glmmTMB(biomass_mice_lm ~ treatment + (1 | plot),
                          dispformula = ~ treatment, 
                          data = arkaute_biomass, family = Gamma(link = "log"))
  em_treat_biomass <- emmeans(glmm_biomass, ~ treatment, type = "response")

  
  
  
  ############# GAM MODELS ##############
  ## Richness ## 
  gam_richness <- gam( richness ~ treatment + s(plot, bs = "re"), data = arkaute_richness, family = gaussian())
  em_treat_richness_gam <- emmeans(gam_richness, ~ treatment)
  
  ## Abundance ##
  gam_abundance <- gam(abundance ~ treatment + s(plot, bs = "re"), data = arkaute_abundance, family = gaussian(link = "identity"))
  em_treat_abundance_gam <- emmeans(gam_abundance, ~ treatment)
 
  ## Evenness ##
  gam_evenness <- gam(Y_zipf ~ treatment + s(plot, bs = "re"), data = arkaute_evenness, family = gaussian())
  em_treat_evenness_gam <- emmeans(gam_evenness, ~ treatment)
  
  ## SLA ##
  gam_sla <- gam(SLA ~ treatment + s(plot, bs = "re"), data = arkaute_sla, family = gaussian(link = "log"))
  em_treat_sla_gam <- emmeans(gam_sla, ~ treatment, type = "response")
  
  ## LDMC ##
  gam_ldmc <- gam(LDMC ~ treatment + s(plot, bs = "re"), data = arkaute_ldmc, family = gaussian())
  em_treat_ldmc_gam <- emmeans(gam_ldmc, ~ treatment)
 
  ## Leaf Nitrogen ##
  gam_leafN <- gam(leafN ~ treatment + s(plot, bs = "re"), data = arkaute_leafN, family = gaussian(link = "identity"))
  em_treat_leafN_gam <- emmeans(gam_leafN, ~ treatment)
 
  ## Biomass ##
  gam_biomass <- gam(biomass_mice_lm ~ treatment + s(plot, bs = "re"), data = arkaute_biomass, family = Gamma(link = "log"))
  em_treat_biomass_gam <-  emmeans(gam_biomass, ~ treatment, type = "response")
  
  
  
  
  #### JOINING TABLES ####
  
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
  
  glmm_list <- list()
  
  glmm_list[[1]] <- as.data.frame(pairs(em_treat_richness, adjust = "tukey")) |>
    mutate(variable = paste0("richness"), AIC = AIC(glmm_richness), estimate_type = "substract")
  
  glmm_list[[2]] <- as.data.frame(pairs(em_treat_abundance, adjust = "tukey")) |>
    mutate(variable = paste0("abundance"), AIC = AIC(glmm_abundance), estimate_type = "substract")
  
  glmm_list[[3]]  <- as.data.frame(pairs(em_treat_evenness, adjust = "tukey"))|>
    mutate(variable = paste0("evenness"), AIC = AIC(glmm_evenness), estimate_type = "substract")
  
  glmm_list[[4]] <- as.data.frame(pairs(em_treat_sla, adjust = "tukey")) |>
    mutate(variable = paste0("SLA"), AIC = AIC(glmm_sla), estimate_type = "ratio"
           )|> rename(estimate = ratio) |>  select(-null)
  
  glmm_list[[5]]  <- as.data.frame(pairs(em_treat_ldmc, adjust = "tukey"))|>
    mutate(variable = paste0("LDMC"), AIC = AIC(glmm_LDMC), estimate_type = "substract")
  
  glmm_list[[6]]  <- as.data.frame(pairs(em_treat_leafN, adjust = "tukey"))|> 
    mutate(variable = paste0("leafN"), AIC = AIC(glmm_leafN), estimate_type = "substract")
  
  glmm_list[[7]]  <- as.data.frame(pairs(em_treat_biomass, adjust = "tukey"))|> 
    mutate(variable = paste0("biomass"), AIC = AIC(glmm_biomass), estimate_type = "ratio"
           )|> rename(estimate = ratio) |> select(-null)
  
  
  glmm_results <- do.call(rbind, glmm_list) |> 
    common_function() |> 
    mutate(model = paste0("GLMM"))|>  
    select(-contrast, -z.ratio)
  
  
  

  gam_list <- list()
  
  gam_list[[1]] <- as.data.frame(pairs(em_treat_richness_gam, adjust = "tukey")) |>
    mutate(variable = "richness", AIC = AIC(gam_richness), estimate_type = "substract")
  
  gam_list[[2]] <- as.data.frame(pairs(em_treat_abundance_gam, adjust = "tukey")) |>
    mutate(variable = "abundance", AIC = AIC(gam_abundance), estimate_type = "substract")
  
  gam_list[[3]] <- as.data.frame(pairs(em_treat_evenness_gam, adjust = "tukey")) |>
    mutate(variable = "evenness", AIC = AIC(gam_evenness), estimate_type = "substract")
  
  gam_list[[4]] <- as.data.frame(pairs(em_treat_sla_gam, adjust = "tukey")) |>
    mutate(variable = "SLA", AIC = AIC(gam_sla), estimate_type = "ratio") |>
    rename(estimate = ratio) |> 
    select(-null)
  
  gam_list[[5]] <- as.data.frame(pairs(em_treat_ldmc_gam, adjust = "tukey")) |>
    mutate(variable = "LDMC", AIC = AIC(gam_ldmc), estimate_type = "substract")
  
  gam_list[[6]] <- as.data.frame(pairs(em_treat_leafN_gam, adjust = "tukey")) |>
    mutate(variable = "leafN", AIC = AIC(gam_leafN), estimate_type = "substract")
  
  gam_list[[7]] <- as.data.frame(pairs(em_treat_biomass_gam, adjust = "tukey")) |>
    mutate(variable = "biomass", AIC = AIC(gam_biomass), estimate_type = "ratio") |>
    rename(estimate = ratio) |> 
    select(-null)
  
  
  gam_results <- do.call(rbind, gam_list) |>
    common_function() |> 
   mutate(model = "GAM")|>
    select(-contrast, -t.ratio)
  
  
 models <- rbind(glmm_results, gam_results) 
  
  
  
  lrr_table <- read.csv("results/effect_size_aggregated.csv") |> 
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
  
  
  
 model_comparison <- full_join(models, lrr_table) |> 
    mutate(variable.bis = variable) |> 
    select(variable, eff_descriptor, model, AIC, 
           p_value, effect_significance, estimate,
           effect_sign, variable.bis
    )
 
 