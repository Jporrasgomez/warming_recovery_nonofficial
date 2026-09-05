



rm(list = ls(all.names = TRUE))


pacman::p_load(
  dplyr, reshape2, tidyverse, lubridate, ggplot2, ggpubr, gridExtra, ggrepel,
  car, ggsignif, dunn.test, rstatix, performance, emmeans, DHARMa, mgcv, glmmTMB
)

source("code/palettes_labels.R")

arkaute <- read.csv("data/processed_data/arkaute.csv") %>%
  mutate(
    year      = factor(year),
    date      = ymd(date),
    sampling = factor(sampling, levels = as.character(sort(unique(as.numeric(as.character(sampling)))))),
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

# Defining MODEL DIAGNOSIS FUNCTION
# 4.2 GLMM diagnostics via DHARMa (including heteroscedasticity)
diagnose_glmm <- function(model, data = NULL, group_var = NULL) {
  print(summary(model))
  sim <- DHARMa::simulateResiduals(fittedModel = model)
  plot(sim)
  print(DHARMa::testDispersion(sim))
  
  if (!is.null(data) && !is.null(group_var)) {
    grp  <- data[[group_var]]
    resu <- sim$scaledResiduals
    cat("Levene test on DHARMa residuals by", group_var, ":\n")
    print(car::leveneTest(resu ~ grp))
    plotResiduals(sim, form = grp)
  }
  invisible(sim)
}

diagnose_gam <- function(model, data = NULL, group_var = NULL) {
  print(summary(model))
  gam.check(model)
  sim <- DHARMa::simulateResiduals(fittedModel = model)
  plot(sim)
  print(DHARMa::testDispersion(sim))
  
  if (!is.null(data) && !is.null(group_var)) {
    grp  <- data[[group_var]]
    resu <- sim$scaledResiduals
    cat("Levene test on DHARMa residuals by", group_var, ":\n")
    print(car::leveneTest(resu ~ grp))
    plotResiduals(sim, form = grp)
  }
  invisible(sim)
}


############# GLMM ##########################
## Richness ## 
glmm_richness <- glmmTMB( richness ~ treatment * sampling + ar1(sampling + 0 | plot),
  data = arkaute_richness,family = gaussian)
#diagnose_glmm(glmm_richness)
em_treat_richness <- emmeans(glmm_richness, ~ treatment, type = "response")
em_time_richness <- emmeans(glmm_richness, ~ treatment | sampling, type = "response")

##### ABUNDANCE #####
# We use gaussian because it tolerates 0 (present in sampling 1 for treatments p and wp)
glmm_abundance <- glmmTMB(abundance ~ treatment * sampling + ar1(sampling + 0 | plot),
  data = arkaute_abundance, family = gaussian(link = "identity"))
#diagnose_glmm(glmm_abundance)
em_treat_abundance <- emmeans(glmm_abundance, ~ treatment, type = "response")
em_time_abundance <- emmeans(glmm_abundance, ~ treatment | sampling, type = "response")

##### EVENNESS #####
# We choose gaussian() family because Y_zipf is a continuos variable of real numbers, and gaussian accepts negative values. Besides,  Y_zipf present a (kind of) symmetrical distribution
glmm_evenness <- glmmTMB(Y_zipf ~ treatment * sampling + ar1(sampling + 0 | plot),
  data = arkaute_evenness, family = gaussian())
#diagnose_glmm(glmm_evenness)
em_treat_evenness <- emmeans(glmm_evenness, ~ treatment, type = "response")
em_time_evenness <- emmeans(glmm_evenness, ~ treatment | sampling, type = "response")

### SLA ###
glmm_sla <- glmmTMB(SLA ~ treatment * sampling + ar1(sampling + 0 | plot),
  data = arkaute_sla, family = gaussian(link = "log"))
#diagnose_glmm(glmm_sla)
em_treat_sla <- emmeans(glmm_sla, ~ treatment, type = "response")
em_time_sla <- emmeans(glmm_sla, ~ treatment | sampling, type = "response")

### LDMC ###
glmm_LDMC <- glmmTMB( LDMC ~ treatment * sampling + ar1(sampling + 0 | plot),
                      data = arkaute_ldmc, family = gaussian())
#diagnose_glmm(glmm_LDMC)
em_treat_ldmc<- emmeans(glmm_LDMC, ~ treatment, type = "response")
em_time_ldmc <- emmeans(glmm_LDMC, ~ treatment | sampling, type = "response")

### Leaf nitrogen ###
glmm_leafN <- glmmTMB(leafN ~ treatment * sampling + ar1(sampling + 0 | plot),
  data = arkaute_leafN, family = gaussian(link = "identity"))
#diagnose_glmm(glmm_leafN)
em_treat_leafN <- emmeans(glmm_leafN, ~ treatment, type = "response")
em_time_leafN <- emmeans(glmm_leafN, ~ treatment | sampling, type = "response")

### BIOMASS ###
glmm_biomass <- glmmTMB(biomass_mice_lm ~ treatment * sampling + ar1(sampling + 0 | plot),
  data = arkaute_biomass,  family = Gamma(link = "log"))
#diagnose_glmm(glmm_biomass)
em_treat_biomass <- emmeans(glmm_biomass, ~ treatment, type = "response")
em_time_biomass <- emmeans(glmm_biomass, ~ treatment | sampling, type = "response")
###################################



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

lrr_reading <- function(data){
  data |> 
    select(-scale, -X) |> 
    rename(diff_value = eff_value) |> 
    mutate(model = paste0("LRR"),
           variable = as.factor(variable)) |> 
    filter(!variable %in% c("biomass_raw", "biomass_mice"),
           sampling != "0") %>%
    mutate(
      variable = fct_recode(variable,
                            "evenness" = "Y_zipf",
                            "biomass"  = "biomass_mice_lm"
      ),
      variable = droplevels(variable),
      effect_sign = ifelse(diff_value > 0 , "positive", "negative"),
      effect_significance = case_when(
        null_effect == "YES" ~ "non-significant",
        TRUE                 ~ "significant", 
      ),
      sampling = as.factor(sampling)
    ) |> 
    select(-upper_limit, -lower_limit, -null_effect)
}

## Aggregating model results at treatment level ##
glmm_list <- list()

  glmm_list[[1]] <- as.data.frame(pairs(em_treat_richness, adjust = "tukey")) |>
    mutate(variable = paste0("richness"), AIC = AIC(glmm_richness), estimate_type = "substract") |> 
    rename(estimate = estimate)

  glmm_list[[2]] <- as.data.frame(pairs(em_treat_abundance, adjust = "tukey")) |> 
    mutate(variable = paste0("abundance"), AIC = AIC(glmm_abundance), estimate_type = "substract")|> 
    rename(estimate = estimate)

  glmm_list[[3]]  <- as.data.frame(pairs(em_treat_evenness, adjust = "tukey")) |>
    mutate(variable = paste0("evenness"), AIC = AIC(glmm_evenness), estimate_type = "substract")|> 
    rename(estimate = estimate)

  glmm_list[[4]] <- as.data.frame(pairs(em_treat_sla, adjust = "tukey")) |> 
    mutate(variable = paste0("SLA"), AIC = AIC(glmm_sla), estimate_type = "ratio")|> 
    rename(estimate = ratio) |> 
    select(-null)

  glmm_list[[5]]  <-  as.data.frame(pairs(em_treat_ldmc, adjust = "tukey")) |>
    mutate(variable = paste0("LDMC"), AIC = AIC(glmm_LDMC), estimate_type = "substract")|> 
    rename(estimate = estimate)

  glmm_list[[6]]  <- as.data.frame(pairs(em_treat_leafN, adjust = "tukey")) |> 
    mutate(variable = paste0("leafN"), AIC = AIC(glmm_leafN), estimate_type = "substract")|> 
    rename(estimate = estimate)

  glmm_list[[7]]  <- as.data.frame(pairs(em_treat_biomass, adjust = "tukey")) |>
    mutate(variable = paste0("biomass"), AIC = AIC(glmm_biomass), estimate_type = "ratio")|> 
    rename(estimate = ratio) |> 
    select(-null)

  
  
glmm_results <- do.call(rbind, glmm_list) |> 
  common_function() |> 
  mutate(model = paste0("GLMM")) |>  
  select(-contrast)

lrr_table <- read.csv("results/effect_size_aggregated.csv") |> 
 lrr_reading()

glmm_agg <- full_join(glmm_result, lrr_table) |> 
  mutate(variable.bis = variable) |> 
  select(variable, eff_descriptor, model, AIC, 
          p_value, effect_significance, estimate,
         effect_sign, variable.bis
  )


## Aggregating model results at sampling  level ##
glmm_time_list <- list()

  glmm_time_list[[1]] <- as.data.frame(pairs(em_time_richness, adjust = "tukey")) |>
    mutate(variable = paste0("richness"), AIC = AIC(glmm_richness), estimate_type = "substract") |> 
    rename(estimate = estimate)
  
  glmm_time_list[[2]] <- as.data.frame(pairs(em_time_abundance, adjust = "tukey")) |> 
    mutate(variable = paste0("abundance"), AIC = AIC(glmm_abundance), estimate_type = "substract")|> 
    rename(estimate = estimate)
  
  glmm_time_list[[3]]  <- as.data.frame(pairs(em_time_evenness, adjust = "tukey")) |>
    mutate(variable = paste0("evenness"), AIC = AIC(glmm_evenness), estimate_type = "substract")|> 
    rename(estimate = estimate)
  
  glmm_time_list[[4]] <- as.data.frame(pairs(em_time_sla, adjust = "tukey")) |> 
    mutate(variable = paste0("SLA"), AIC = AIC(glmm_sla), estimate_type = "ratio")|> 
    rename(estimate = ratio) |> 
    select(-null)
  
  glmm_time_list[[5]]  <-  as.data.frame(pairs(em_time_ldmc, adjust = "tukey")) |>
    mutate(variable = paste0("LDMC"), AIC = AIC(glmm_LDMC), estimate_type = "substract")|> 
    rename(estimate = estimate)
  
  glmm_time_list[[6]]  <- as.data.frame(pairs(em_time_leafN, adjust = "tukey")) |> 
    mutate(variable = paste0("leafN"), AIC = AIC(glmm_leafN), estimate_type = "substract")|> 
    rename(estimate = estimate)
  
  glmm_time_list[[7]]  <- as.data.frame(pairs(em_time_biomass, adjust = "tukey")) |>
    mutate(variable = paste0("biomass"), AIC = AIC(glmm_biomass), estimate_type = "ratio")|> 
    rename(estimate = ratio) |> 
    select(-null)
  
glmm_time <- do.call(rbind, glmm_time_list) |> 
  common_function() |> 
  mutate(model = paste0("GLMM")) |>  
  select(-contrast)

lrr_table_dyn <- read.csv("results/effect_size_dynamics.csv") |> 
  lrr_reading()

glmm_dynamics <- full_join(glmm_time, lrr_table_dyn) |> 
  mutate(variable.bis = variable) |> 
  select(variable, sampling, eff_descriptor, model, 
         p_value, effect_significance, diff_value,
         effect_sign, variable.bis) |> 
  mutate(
    variable_model = paste0(variable, "-", model)
  )









glmm_agg |> 
  filter(eff_descriptor == "wp_vs_p") |> 
  ggplot(aes(x = model, y = variable, color = effect_significance, shape = effect_sign,
             label = round(p_value, 2))) +
  geom_point(size = 10) +
  #geom_label_repel(show.legend = FALSE) +
  scale_color_manual(values = palette_significance) +
  scale_shape_manual(values = shape_significance) +
  labs(title = "Warming effect on recovery (wp vs p)", x = "Statistical analysis", y = NULL,
       shape = "Sign of effect", color = "Effect stat. significance")

glmm_agg |> 
  filter(eff_descriptor == "p_vs_c") |> 
  ggplot(aes(x = model, y = variable, color = effect_significance, shape = effect_sign,
             label = round(p_value, 2))) +
  geom_point(size = 10) +
  geom_label_repel(show.legend = FALSE) +
  scale_color_manual(values = palette_significance) +
  scale_shape_manual(values = shape_significance) +
  labs(title = "Perturbation effect (p vs c)", x = "Statistical analysis", y = NULL,
       shape = "Sign of effect", color = "Effect stat. significance")


glmm_agg |> 
  filter(eff_descriptor == "wp_vs_c") |> 
  ggplot(aes(x = model, y = variable, color = effect_significance, shape = effect_sign,
             label = round(p_value, 2))) +
  geom_point(size = 10) +
  #geom_label_repel(show.legend = FALSE) +
  scale_color_manual(values = palette_significance) +
  scale_shape_manual(values = shape_significance) +
  labs(title = "Combined effect (wp vs c)", x = "Statistical analysis", y = NULL,
       shape = "Sign of effect", color = "Effect stat. significance")


glmm_agg |> 
  filter(eff_descriptor == "w_vs_c") |> 
  ggplot(aes(x = model, y = variable, color = effect_significance, shape = effect_sign,
             label = round(p_value, 2))) +
  geom_point(size = 10) +
  #geom_label_repel(show.legend = FALSE) +
  scale_color_manual(values = palette_significance) +
  scale_shape_manual(values = shape_significance) +
  labs(title = "Warming effect on assembly (w vs c)", x = "Statistical analysis", y = NULL,
       shape = "Sign of effect", color = "Effect stat. significance")


# Juntar en una columna variable-model para usar en el eje Y y así verlo todo junto
glmm_dynamics |> 
  filter(eff_descriptor == "wp_vs_p") |> 
  ggplot(aes(x = sampling, y = variable_model, color = effect_significance, shape = effect_sign)) +
  #facet_wrap(~model) + 
  geom_point(size = 8) +
  scale_color_manual(values = palette_sig) +
  scale_shape_manual(values = palette_shape) +
  labs(title = "Warming effect on recovery (wp/p)", 
       x = "Sampling", 
       y = "Variable and model",
       color = "Effect significance",
       shape = "Sign of the effect")

glmm_dynamics |> 
  filter(eff_descriptor == "p_vs_c") |> 
  ggplot(aes(x = sampling, y = variable_model, color = effect_significance, shape = effect_sign)) +
  #facet_wrap(~model) + 
  geom_point(size = 8) +
  scale_color_manual(values = palette_sig) +
  scale_shape_manual(values = palette_shape) +
  labs(title = "Recovery(p/c)", 
       x = "Sampling", 
       y = "Variable and model")

glmm_dynamics |> 
  filter(eff_descriptor == "w_vs_c") |> 
  ggplot(aes(x = sampling, y = variable_model, color = effect_significance, shape = effect_sign)) +
  #facet_wrap(~model) + 
  geom_point(size = 8) +
  scale_color_manual(values = palette_sig) +
  scale_shape_manual(values = palette_shape) +
  labs(title = "Warming effect(w/c)", 
       x = "Sampling", 
       y = "Variable and model")


glmm_dynamics |> 
  filter(eff_descriptor == "wp_vs_c") |> 
  ggplot(aes(x = sampling, y = variable_model, color = effect_significance, shape = effect_sign)) +
  #facet_wrap(~model) + 
  geom_point(size = 8) +
  scale_color_manual(values = palette_sig) +
  scale_shape_manual(values = palette_shape) +
  labs(title = "Combined effect(wp/c)", 
       x = "Sampling", 
       y = "Variable and model")

