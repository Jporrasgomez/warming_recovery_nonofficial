

rm(list = ls(all.names = TRUE))

# Cargar paquetes
pacman::p_load(
  dplyr, reshape2, tidyverse, lubridate, ggplot2, ggpubr, gridExtra, ggrepel,
  car, ggsignif, dunn.test, rstatix, performance, emmeans, DHARMa, mgcv
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
    year         = factor(year),
    date         = ymd(date),
    sampling_num = as.numeric(as.character(sampling)), # Necesario para splines en GAM
    sampling     = factor(sampling, levels = as.character(sort(unique(as.numeric(as.character(sampling)))))),
    plot         = factor(plot),
    treatment    = factor(treatment)
  ) %>%
  filter(sampling != "0") |> 
  arrange(plot, sampling_num)

# Subconjuntos filtrando muestreo 1 y NAs
arkaute_richness  <- arkaute |> filter(sampling != "1")
arkaute_abundance <- arkaute |> filter(sampling != "1")
arkaute_evenness  <- arkaute |> filter(sampling != "1", !is.na(Y_zipf))
arkaute_sla       <- arkaute |> filter(sampling != "1", !is.na(SLA))
arkaute_ldmc      <- arkaute |> filter(sampling != "1", !is.na(LDMC))
arkaute_leafN     <- arkaute |> filter(sampling != "1", !is.na(leafN))
arkaute_biomass   <- arkaute |> filter(sampling != "1")

# Función de diagnóstico para GAM via DHARMa y mgcv
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

# Muestreos a evaluar en dinámicas temporales
samplings_eval <- sort(unique(arkaute_evenness$sampling_num))

{
  ############# GAM MODELS ##############
  
  ## Richness ## 
  gam_richness <- gam(
    richness ~ treatment + s(sampling_num, by = treatment, k = 10) + s(plot, bs = "re"),
    data = arkaute_richness,
    family = gaussian()
  )
  diagnose_gam(gam_richness)
  
  em_treat_richness <- emmeans(gam_richness, ~ treatment)
  em_time_richness  <- emmeans(gam_richness, ~ treatment | sampling_num, at = list(sampling_num = samplings_eval))
  
  ## Abundance ##
  gam_abundance <- gam(
    abundance ~ treatment + s(sampling_num, by = treatment, k = 10) + s(plot, bs = "re"),
    data = arkaute_abundance,
    family = gaussian(link = "identity")
  )
  diagnose_gam(gam_abundance)
  
  em_treat_abundance <- emmeans(gam_abundance, ~ treatment)
  em_time_abundance  <- emmeans(gam_abundance, ~ treatment | sampling_num, at = list(sampling_num = samplings_eval))
  
  ## Evenness ##
  gam_evenness <- gam(
    Y_zipf ~ treatment + s(sampling_num, by = treatment, k = 10) + s(plot, bs = "re"),
    data = arkaute_evenness,
    family = gaussian()
  )
  diagnose_gam(gam_evenness)
  
  em_treat_evenness <- emmeans(gam_evenness, ~ treatment)
  em_time_evenness  <- emmeans(gam_evenness, ~ treatment | sampling_num, at = list(sampling_num = samplings_eval))
  
  ## SLA ##
  gam_sla <- gam(
    SLA ~ treatment + s(sampling_num, by = treatment, k = 10) + s(plot, bs = "re"),
    data = arkaute_sla,
    family = gaussian(link = "log")
  )
  diagnose_gam(gam_sla)
  
  em_treat_sla <- emmeans(gam_sla, ~ treatment, type = "response")
  em_time_sla  <- emmeans(gam_sla, ~ treatment | sampling_num, at = list(sampling_num = samplings_eval), type = "response")
  
  ## LDMC ##
  gam_ldmc <- gam(
    LDMC ~ treatment + s(sampling_num, by = treatment, k = 10) + s(plot, bs = "re"),
    data = arkaute_ldmc,
    family = gaussian()
  )
  diagnose_gam(gam_ldmc)
  
  em_treat_ldmc <- emmeans(gam_ldmc, ~ treatment)
  em_time_ldmc  <- emmeans(gam_ldmc, ~ treatment | sampling_num, at = list(sampling_num = samplings_eval))
  
  ## Leaf Nitrogen ##
  gam_leafN <- gam(
    leafN ~ treatment + s(sampling_num, by = treatment, k = 10) + s(plot, bs = "re"),
    data = arkaute_leafN,
    family = gaussian(link = "identity")
  )
  diagnose_gam(gam_leafN)
  
  em_treat_leafN <- emmeans(gam_leafN, ~ treatment)
  em_time_leafN  <- emmeans(gam_leafN, ~ treatment | sampling_num, at = list(sampling_num = samplings_eval))
  
  ## Biomass ##
  gam_biomass <- gam(
    biomass_mice_lm ~ treatment + 
      s(sampling_num, by = treatment, k = 10) + s(plot, bs = "re"),
    data = arkaute_biomass,
    family = tw(link = "log")
  )
  
diagnose_gam(gam_biomass)
em_treat_biomass <-  emmeans(gam_biomass, ~ treatment, type = "response")
em_time_biomass  <- emmeans(gam_biomass, ~ treatment | sampling_num, at = list(sampling_num = samplings_eval), type = "response")
#pairs(em_treat_biomass, adjust = "tukey")
  
}

# --- Extracción y procesamiento de resultados agregados GAM ---
{
  model_list <- list()
  
  model_list[[1]] <- as.data.frame(pairs(em_treat_richness, adjust = "tukey")) |>
    mutate(variable = "richness", AIC = AIC(gam_richness), estimate_type = "substract")
  
  model_list[[2]] <- as.data.frame(pairs(em_treat_abundance, adjust = "tukey")) |>
    mutate(variable = "abundance", AIC = AIC(gam_abundance), estimate_type = "substract")
  
  model_list[[3]] <- as.data.frame(pairs(em_treat_evenness, adjust = "tukey")) |>
    mutate(variable = "evenness", AIC = AIC(gam_evenness), estimate_type = "substract")
  
  model_list[[4]] <- as.data.frame(pairs(em_treat_sla, adjust = "tukey")) |>
    mutate(variable = "SLA", AIC = AIC(gam_sla), estimate_type = "ratio") |>
    rename(estimate = ratio) |> 
    select(-null)
  
  model_list[[5]] <- as.data.frame(pairs(em_treat_ldmc, adjust = "tukey")) |>
    mutate(variable = "LDMC", AIC = AIC(gam_ldmc), estimate_type = "substract")
  
  model_list[[6]] <- as.data.frame(pairs(em_treat_leafN, adjust = "tukey")) |>
    mutate(variable = "leafN", AIC = AIC(gam_leafN), estimate_type = "substract")
  
  model_list[[7]] <- as.data.frame(pairs(em_treat_biomass, adjust = "tukey")) |>
    mutate(variable = "biomass", AIC = AIC(gam_biomass), estimate_type = "ratio") |>
    rename(estimate = ratio) |> 
    select(-null)
  
  model_result <- do.call(rbind, model_list) |>
    filter(!contrast %in% c("p - w", "p / w", "w - wp", "w / wp")) |>
    rename(p_value = p.value) |>
    mutate(
      effect_sign = case_when(
        estimate_type == "substract" & estimate < 0     ~ "positive",
        estimate_type == "substract" & estimate > 0     ~ "negative",
        estimate_type == "ratio"     & estimate > 1     ~ "negative",
        estimate_type == "ratio"     & estimate < 1     ~ "positive"
      ),
      effect_significance = case_when(
        p_value < 0.05                   ~ "significant",
        p_value >= 0.05 & p_value < 0.10 ~ "marginal",
        TRUE                             ~ "non-significant"
      ),
      eff_descriptor = case_when(
        contrast %in% c("c - p", "c / p")   ~ "p_vs_c",
        contrast %in% c("c - w", "c / w")   ~ "w_vs_c",
        contrast %in% c("c - wp", "c / wp") ~ "wp_vs_c",
        contrast %in% c("p - wp", "p / wp") ~ "wp_vs_p"
      ),
      model = "GAM"
    ) |>
    select(-contrast)
  
  # Cargar tabla agregada LRR
  lrr_table <- read.csv("results/effect_size_aggregated.csv") |>
    select(-scale, -X) |>
    rename(estimate = eff_value) |>
    mutate(
      model = "LRR",
      variable = as.factor(variable)
    ) |>
    filter(!variable %in% c("biomass_raw", "biomass_mice")) %>%
    mutate(
      variable = fct_recode(variable,
                            "evenness" = "Y_zipf",
                            "biomass"  = "biomass_mice_lm"
      ),
      variable = droplevels(variable),
      effect_sign = ifelse(estimate > 0, "positive", "negative"),
      effect_significance = case_when(
        null_effect == "YES" ~ "non-significant",
        TRUE                 ~ "significant"
      )
    ) |>
    select(-upper_limit, -lower_limit, -null_effect)
  
  models <- full_join(model_result, lrr_table) |>
    mutate(variable.bis = variable) |>
    select(
      variable, eff_descriptor, model, AIC,
      p_value, effect_significance, estimate,
      effect_sign, variable.bis
    )
  
  palette_sig <- c("significant" = "blue", "marginal" = "orange", "non-significant" = "grey")
  palette_shape <- c("positive" = "+", "negative" = "-")
}

# --- Gráficos Comparativos Agregados (GAM vs LRR) ---
models |>
  filter(eff_descriptor == "wp_vs_p") |>
  ggplot(aes(x = model, y = variable, color = effect_significance, shape = effect_sign)) +
  geom_point(size = 10) +
  scale_color_manual(values = palette_sig) +
  scale_shape_manual(values = palette_shape) +
  labs(title = "Warming effect on recovery (wp vs p)", x = "Statistical analysis", y = NULL,
       shape = "Sign of effect", color = "Effect stat. significance")

models |>
  filter(eff_descriptor == "p_vs_c") |>
  ggplot(aes(x = model, y = variable, color = effect_significance, shape = effect_sign)) +
  geom_point(size = 10) +
  scale_color_manual(values = palette_sig) +
  scale_shape_manual(values = palette_shape) +
  labs(title = "Perturbation effect (p vs c)", x = "Statistical analysis", y = NULL,
       shape = "Sign of effect", color = "Effect stat. significance")

models |>
  filter(eff_descriptor == "wp_vs_c") |>
  ggplot(aes(x = model, y = variable, color = effect_significance, shape = effect_sign)) +
  geom_point(size = 10) +
  scale_color_manual(values = palette_sig) +
  scale_shape_manual(values = palette_shape) +
  labs(title = "Combined effect (wp vs c)", x = "Statistical analysis", y = NULL,
       shape = "Sign of effect", color = "Effect stat. significance")

models |>
  filter(eff_descriptor == "w_vs_c") |>
  ggplot(aes(x = model, y = variable, color = effect_significance, shape = effect_sign)) +
  geom_point(size = 10) +
  scale_color_manual(values = palette_sig) +
  scale_shape_manual(values = palette_shape) +
  labs(title = "Warming effect on assembly (w vs c)", x = "Statistical analysis", y = NULL,
       shape = "Sign of effect", color = "Effect stat. significance")

# --- Extracción y procesamiento de dinámicas temporales GAM ---
{
  model_time_result <- list()
  
  model_time_result[[1]] <- as.data.frame(pairs(em_time_richness, adjust = "tukey")) |>
    mutate(variable = "richness", AIC = AIC(gam_richness), estimate_type = "substract")
  
  model_time_result[[2]] <- as.data.frame(pairs(em_time_abundance, adjust = "tukey")) |>
    mutate(variable = "abundance", AIC = AIC(gam_abundance), estimate_type = "substract")
  
  model_time_result[[3]] <- as.data.frame(pairs(em_time_evenness, adjust = "tukey")) |>
    mutate(variable = "evenness", AIC = AIC(gam_evenness), estimate_type = "substract")
  
  model_time_result[[4]] <- as.data.frame(pairs(em_time_sla, adjust = "tukey")) |>
    mutate(variable = "SLA", AIC = AIC(gam_sla), estimate_type = "ratio")|>
    rename(estimate = ratio) |> 
    select(-null)
  
  model_time_result[[5]] <- as.data.frame(pairs(em_time_ldmc, adjust = "tukey")) |>
    mutate(variable = "LDMC", AIC = AIC(gam_ldmc), estimate_type = "substract")
  
  model_time_result[[6]] <- as.data.frame(pairs(em_time_leafN, adjust = "tukey")) |>
    mutate(variable = "leafN", AIC = AIC(gam_leafN), estimate_type = "substract")
  
  model_time_result[[7]] <- as.data.frame(pairs(em_time_biomass, adjust = "tukey")) |>
    mutate(variable = "biomass", AIC = AIC(gam_biomass), estimate_type = "ratio") |>
    rename(estimate = ratio) |> 
    select(-null)
  
  model_time <- do.call(rbind, model_time_result) |>
    filter(!contrast %in% c("p - w", "p / w", "w - wp", "w / wp")) |>
    rename(p_value = p.value, sampling = sampling_num) |>
    mutate(
      sampling = as.factor(sampling),
      effect_sign = case_when(
        estimate_type == "substract" & estimate < 0     ~ "positive",
        estimate_type == "substract" & estimate > 0     ~ "negative",
        estimate_type == "ratio"     & estimate > 1     ~ "negative",
        estimate_type == "ratio"     & estimate < 1     ~ "positive"
      ),
      effect_significance = case_when(
        p_value < 0.05                   ~ "significant",
        p_value >= 0.05 & p_value < 0.10 ~ "marginal",
        TRUE                             ~ "non-significant"
      ),
      eff_descriptor = case_when(
        contrast %in% c("c - p", "c / p")   ~ "p_vs_c",
        contrast %in% c("c - w", "c / w")   ~ "w_vs_c",
        contrast %in% c("c - wp", "c / wp") ~ "wp_vs_c",
        contrast %in% c("p - wp", "p / wp") ~ "wp_vs_p"
      ),
      model = "GAM"
    ) |>
    select(-contrast)
  
  # Cargar dinámicas LRR
  lrr_table_dyn <- read.csv("results/effect_size_dynamics.csv") |>
    select(-scale, -X) |>
    rename(estimate = eff_value) |>
    mutate(
      model = "LRR",
      variable = as.factor(variable)
    ) |>
    filter(
      !variable %in% c("biomass_raw", "biomass_mice"),
      sampling != "0"
    ) %>%
    mutate(
      variable = fct_recode(variable,
                            "evenness" = "Y_zipf",
                            "biomass"  = "biomass_mice_lm"
      ),
      variable = droplevels(variable),
      effect_sign = ifelse(estimate > 0, "positive", "negative"),
      effect_significance = case_when(
        null_effect == "YES" ~ "non-significant",
        TRUE                 ~ "significant"
      ),
      sampling = as.factor(sampling)
    ) |>
    select(-upper_limit, -lower_limit, -null_effect)
  
  models_dynamics <- full_join(model_time, lrr_table_dyn) |>
    mutate(
      variable.bis = variable,
      variable_model = paste0(variable, "-", model)
    ) |>
    select(
      variable, sampling, eff_descriptor, model,
      p_value, effect_significance, estimate,
      effect_sign, variable.bis, variable_model
    )
  }

# --- Gráficos Dinámicos Temporales (GAM vs LRR) ---
models_dynamics |>
  filter(eff_descriptor == "wp_vs_p") |>
  ggplot(aes(x = sampling, y = variable_model, color = effect_significance, shape = effect_sign)) +
  geom_point(size = 8) +
  scale_color_manual(values = palette_sig) +
  scale_shape_manual(values = palette_shape) +
  labs(title = "Warming effect on recovery (wp/p)", x = "Sampling", y = "Variable and model",
       color = "Effect significance", shape = "Sign of the effect")

models_dynamics |>
  filter(eff_descriptor == "p_vs_c") |>
  ggplot(aes(x = sampling, y = variable_model, color = effect_significance, shape = effect_sign)) +
  geom_point(size = 8) +
  scale_color_manual(values = palette_sig) +
  scale_shape_manual(values = palette_shape) +
  labs(title = "Recovery (p/c)", x = "Sampling", y = "Variable and model")

models_dynamics |>
  filter(eff_descriptor == "w_vs_c") |>
  ggplot(aes(x = sampling, y = variable_model, color = effect_significance, shape = effect_sign)) +
  geom_point(size = 8) +
  scale_color_manual(values = palette_sig) +
  scale_shape_manual(values = palette_shape) +
  labs(title = "Warming effect (w/c)", x = "Sampling", y = "Variable and model")

models_dynamics |>
  filter(eff_descriptor == "wp_vs_c") |>
  ggplot(aes(x = sampling, y = variable_model, color = effect_significance, shape = effect_sign)) +
  geom_point(size = 8) +
  scale_color_manual(values = palette_sig) +
  scale_shape_manual(values = palette_shape) +
  labs(title = "Combined effect (wp/c)", x = "Sampling", y = "Variable and model")



