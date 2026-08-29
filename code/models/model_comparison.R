



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
    sampling  = factor(sampling),
    plot      = factor(plot),
    treatment = factor(treatment)
  ) %>%
  filter(sampling != "0") %>%
  filter(!(sampling == "1" & treatment %in% c("p", "wp")))





# Diagnosis functions


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
  sim <- DHARMa::simulateResiduals(fittedModel = model)
  plot(sim)
  print(DHARMa::testDispersion(sim))
  gam.check(model)
  
  if (!is.null(data) && !is.null(group_var)) {
    grp  <- data[[group_var]]
    resu <- sim$scaledResiduals
    cat("Levene test on DHARMa residuals by", group_var, ":\n")
    print(car::leveneTest(resu ~ grp))
    plotResiduals(sim, form = grp)
  }
  invisible(sim)
}




# =================================================================================
# Función de Extracción Unificada con el orden de columnas ajustado
# =================================================================================

extract_model_diagnostics_and_contrasts <- 
  function(model, var_name, model_type = "GAM", data = NULL, group_var = "treatment") {
  
  # 1. Simulación de residuos con DHARMa
  sim <- DHARMa::simulateResiduals(fittedModel = model, plot = FALSE)
  
  # Test de Dispersión (DHARMa)
  disp_test <- DHARMa::testDispersion(sim, plot = FALSE)
  disp_status <- if (disp_test$p.value >= 0.05) "NS" else "problem"
  
  # Test Uniformidad / KS (DHARMa)
  ks_test <- DHARMa::testUniformity(sim, plot = FALSE)
  ks_status <- if (ks_test$p.value >= 0.05) "NS" else "problem"
  
  # Test de Levene (Homocedasticidad entre grupos)
  levene_status <- "NS"
  if (!is.null(data) && !is.null(group_var)) {
    grp  <- data[[group_var]]
    resu <- sim$scaledResiduals
    lev  <- car::leveneTest(resu ~ grp)
    p_lev <- lev$`Pr(>F)`[1]
    levene_status <- if (!is.na(p_lev) && p_lev >= 0.05) "NS" else "problem"
  }
  
  # 2. Emmeans y Contrastes
  emm <- emmeans(model, pairwise ~ treatment)
  contrasts_df <- as.data.frame(emm$contrasts)
  
  ratio_col <- if ("t.ratio" %in% colnames(contrasts_df)) "t.ratio" else "z.ratio"
  target_contrasts <- c("c - p", "c - w", "c - wp", "p - wp")
  
  # Extraer metadatos
  if (model_type == "GAM") {
    fam_name  <- model$family$family
    link_name <- model$family$link
  } else {
    fam_info  <- family(model)
    fam_name  <- fam_info$family
    link_name <- fam_info$link
  }
  
  # Construir dataframe final con el orden exacto especificado
  res_df <- contrasts_df %>%
    filter(contrast %in% target_contrasts) %>%
    mutate(
      variable     = var_name,
      model        = model_type,
      family       = fam_name,
      link         = link_name,
      dispersion   = disp_status,
      dharma_ks    = ks_status,
      levene_test  = levene_status,
      contrast     = gsub(" ", "", contrast),
      ratio        = .data[[ratio_col]],
      p_value      = p.value
    ) %>%
    select(variable, model, family, link, dispersion, dharma_ks,
           levene_test, contrast, ratio, p_value)
  
  return(res_df)
}

# =================================================================================
# Preparación de datasets filtrados
# =================================================================================



arkaute_abundance <- arkaute |>  filter (!is.na(abundance))
arkaute_Yzipf   <- arkaute %>% filter(!is.na(Y_zipf), Y_zipf < -1e-5
                                     ) %>% mutate(Y_zipf = Y_zipf * -1)
arkaute_biomass <- arkaute %>% filter(!is.na(biomass_mice_lm))
arkaute_SLA     <- arkaute %>% filter(!is.na(SLA))
arkaute_LDMC    <- arkaute %>% filter(!is.na(LDMC))
arkaute_leafN   <- arkaute %>% filter(!is.na(leafN))

# =================================================================================
#  GAM
# =================================================================================

results_gam_list <- list()

# Richness
gam_richness <- 
  gam(richness ~ treatment + s(sampling, bs = "re"),
      data = arkaute,
      family = gaussian(link = "identity"))
diagnose_gam(gam_richness)
emmeans(gam_richness, pairwise ~ treatment)
results_gam_list[["richness"]] <-
  extract_model_diagnostics_and_contrasts(gam_richness, "richness", "GAM", arkaute)

# Abundance

gam_abundance <- gam(abundance ~ treatment + s(sampling, bs = "re"),
                     data = abundance_noOutliers,
                     family = gaussian(link = "log"))
diagnose_gam(gam_abundance)
emmeans(gam_abundance, pairwise ~ treatment)
results_gam_list[["abundance"]] <-
  extract_model_diagnostics_and_contrasts(gam_abundance, "abundance", "GAM", abundance_noOutliers)


gam_evenness <- gam(Y_zipf ~ treatment + s(sampling, bs = "re"),
                    data = arkaute_Yzipf,
                    family = Gamma(link = "log"))
diagnose_gam(gam_evenness)
emmeans(gam_abundance, pairwise ~ treatment)
results_gam_list[["evenness"]] <-
  extract_model_diagnostics_and_contrasts(gam_evenness, "evenness", "GAM", arkaute_Yzipf)


gam_biomass <- gam(biomass_mice_lm ~ treatment + s(sampling, bs = "re"),
                   data = arkaute_biomass,
                   family = gaussian(link = "identity"))
diagnose_gam(gam_biomass)
emmeans(gam_biomass, pairwise ~ treatment)
results_gam_list[["biomass"]] <-
  extract_model_diagnostics_and_contrasts(gam_biomass, "biomass", "GAM", arkaute_biomass)

gam_SLA <- gam(SLA ~ treatment + s(sampling, bs = "re"),
               data = arkaute_SLA,
               family = gaussian(link = "log"))
diagnose_gam(gam_SLA)
emmeans(gam_SLA, pairwise ~ treatment)
results_gam_list[["SLA"]] <-
  extract_model_diagnostics_and_contrasts(gam_SLA, "SLA", "GAM", arkaute_SLA)

gam_LDMC <- gam(LDMC ~ treatment + s(sampling, bs = "re"),
                data = arkaute_LDMC,
                family = gaussian(link = "identity"))
diagnose_gam(gam_LDMC)
emmeans(gam_LDMC, pairwise ~ treatment)
results_gam_list[["LDMC"]] <-
  extract_model_diagnostics_and_contrasts(gam_LDMC, "LDMC", "GAM", arkaute_LDMC)


gam_leafN <- gam(leafN ~ treatment + s(sampling, bs = "re"),
                 data = arkaute_leafN,
                 family = gaussian(link = "identity"))
diagnose_gam(gam_leafN)
emmeans(gam_leafN, pairwise ~ treatment)
results_gam_list[["leafN"]] <-
  extract_model_diagnostics_and_contrasts(gam_leafN, "leafN", "GAM", arkaute_leafN)

emmeans_gam_summary_table <- bind_rows(results_gam_list)


# =================================================================================
#  GLMM
# =================================================================================

results_glmm_list <- list()

glmm_richness <- glmmTMB(richness ~ treatment + (1 | sampling),
                         data = arkaute,
                         family = genpois(link = "log"))
diagnose_glmm(glmm_richness)
emmeans(glmm_richness, pairwise ~ treatment)
emmeans(glmm_richness, trt.vs.ctrl ~ treatment, ref = "p", type = "response")
results_glmm_list[["richness"]] <-
  extract_model_diagnostics_and_contrasts(glmm_richness, "richness", "GLMM", arkaute)


glmm_abundance <- glmmTMB(abundance ~ treatment + (1 | sampling),
                          data = arkaute_abundance,
                          family = tweedie(link = "log"))
diagnose_glmm(glmm_abundance)
emmeans(glmm_abundance, pairwise ~ treatment)
emmeans(glmm_abundance, trt.vs.ctrl ~ treatment, ref = "p", type = "response")
results_glmm_list[["abundance"]] <-
  extract_model_diagnostics_and_contrasts(glmm_abundance, "abundance", "GLMM", abundance_noOutliers)


glmm_evenness <- glmmTMB(Y_zipf ~ treatment + (1 | sampling),
                         data = arkaute_Yzipf,
                         family = Gamma(link = "log"))
diagnose_glmm(glmm_evenness)
emmeans(glmm_evenness, pairwise ~ treatment)
emmeans(glmm_evenness, trt.vs.ctrl ~ treatment, ref = "p", type = "response")
results_glmm_list[["evenness"]] <-
  extract_model_diagnostics_and_contrasts(glmm_evenness, "evenness", "GLMM", arkaute_Yzipf)


glmm_biomass <- glmmTMB(biomass_mice_lm ~ treatment + (1 | sampling),
                        data = arkaute_biomass,
                        family = Gamma(link = "log"))
diagnose_glmm(glmm_biomass)
emmeans(glmm_biomass, pairwise ~ treatment)
emmeans(glmm_biomass, trt.vs.ctrl ~ treatment, ref = "p", type = "response")
results_glmm_list[["biomass"]] <-
  extract_model_diagnostics_and_contrasts(glmm_biomass, "biomass", "GLMM", arkaute_biomass)


glmm_SLA <- glmmTMB(SLA ~ treatment + (1 | sampling),
                    data = arkaute_SLA,
                    family = Gamma(link = "log"))
diagnose_glmm(glmm_SLA)
emmeans(glmm_SLA, pairwise ~ treatment)
emmeans(glmm_SLA, trt.vs.ctrl ~ treatment, ref = "p", type = "response")
results_glmm_list[["SLA"]] <-
  extract_model_diagnostics_and_contrasts(glmm_SLA, "SLA", "GLMM", arkaute_SLA)


glmm_LDMC <- glmmTMB(LDMC ~ treatment + (1 | sampling),
                     data = arkaute_LDMC,
                     family = gaussian(link = "identity"))
diagnose_glmm(glmm_LDMC)
emmeans(glmm_LDMC, pairwise ~ treatment)
emmeans(glmm_LDMC, trt.vs.ctrl ~ treatment, ref = "p", type = "response")
results_glmm_list[["LDMC"]] <- extract_model_diagnostics_and_contrasts(glmm_LDMC, "LDMC", "GLMM", arkaute_LDMC)


glmm_leafN <- glmmTMB(leafN ~ treatment + (1 | sampling),
                      data = arkaute_leafN,
                      family = gaussian(link = "identity"))
diagnose_glmm(glmm_leafN)
emmeans(glmm_leafN, pairwise ~ treatment)
emmeans(glmm_leafN, trt.vs.ctrl ~ treatment, ref = "p", type = "response")
results_glmm_list[["leafN"]] <- extract_model_diagnostics_and_contrasts(glmm_leafN, "leafN", "GLMM", arkaute_leafN)

emmeans_glmm_summary_table <-
  bind_rows(results_glmm_list)


# =================================================================================
# TABLA COMPARATIVA CONSOLIDADA
# =================================================================================

emmeans_table <- bind_rows(emmeans_gam_summary_table, emmeans_glmm_summary_table) %>%
  arrange(variable, contrast, model) |> 
  mutate(
    effect_sign = ifelse(ratio > 0 , "negative", "positive"), 
    effect_significance = case_when(
      p_value < 0.05                     ~ "significant",
      p_value >= 0.05  & p_value < 0.10  ~ "marginal",
      TRUE                               ~ "non-significant"
    ),
    eff_descriptor = case_when(
      contrast == "c-p"  ~ "p_vs_c", 
      contrast == "c-w"  ~ "w_vs_c", 
      contrast == "c-wp" ~ "wp_vs_c", 
      contrast == "p-wp" ~ "wp_vs_p"
    ), 
    variable = fct_recode(variable,
                            "dominance" = "evenness",
      )
  ) |> 
  rename(diff_value = ratio) |> 
  select(-contrast)



# Visualizar la tabla resultante
#print(emmeans_table)


lrr_table <- read.csv("results/effect_size_aggregated.csv") |> 
  select(-scale, -X) |> 
  rename(diff_value = eff_value) |> 
  mutate(model = paste0("LRR"),
         variable = as.factor(variable)) |> 
  filter(!variable %in% c("biomass_raw", "biomass_mice")) %>%
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
    )
  ) |> 
  select(-upper_limit, -lower_limit, -null_effect)




models <- full_join(emmeans_table, lrr_table) |> 
  mutate(variable.bis = variable) |> 
  select(variable, eff_descriptor, model, family, link, dispersion, 
         dharma_ks, levene_test, p_value, effect_significance, diff_value,
          effect_sign, variable.bis
         )

#models |>  write.csv("results/model_comparison.csv")

models_p_vs_c  <- models |> filter (eff_descriptor == "p_vs_c")
models_w_vs_c  <- models |> filter (eff_descriptor == "w_vs_c")
models_wp_vs_c <- models |> filter (eff_descriptor == "wp_vs_c")
models_wp_vs_p <- models |> filter (eff_descriptor == "wp_vs_p")

palette_sig <- 
  c("significant" = "blue", "marginal" = "orange", "non-significant" = "grey")
palette_shape <- c(
  "positive" = "+",  # o pch 43 / 3
  "negative" = "-"   # o pch 45
)

models |> 
  filter(eff_descriptor == "wp_vs_p") |> 
  ggplot(aes(x = model, y = variable, color = effect_significance, shape = effect_sign,
             label = round(p_value, 2))) +
  geom_point(size = 10) +
  geom_label_repel(show.legend = FALSE) +
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

