

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
                          data = arkaute_richness,family = gaussian)
#diagnose_glmm(glmm_richness)
em_treat_richness <- emmeans(glmm_richness, ~ treatment, type = "response")
em_time_richness <- emmeans(glmm_richness, ~ treatment | sampling, type = "response")

## ABUNDANCE ##
# We use gaussian because it tolerates 0 (present in sampling 1 for treatments p and wp)
glmm_abundance <- glmmTMB(abundance ~ treatment * sampling + ar1(sampling + 0 | plot),
                          data = arkaute_abundance, family = gaussian(link = "identity"))
#diagnose_glmm(glmm_abundance)
em_treat_abundance <- emmeans(glmm_abundance, ~ treatment, type = "response")
em_time_abundance <- emmeans(glmm_abundance, ~ treatment | sampling, type = "response")

## EVENNESS ##
# We choose gaussian() family because Y_zipf is a continuos variable of real numbers, and gaussian accepts negative values. Besides,  Y_zipf present a (kind of) symmetrical distribution
glmm_evenness <- glmmTMB(Y_zipf ~ treatment * sampling + ar1(sampling + 0 | plot),
                         data = arkaute_evenness, family = gaussian())
#diagnose_glmm(glmm_evenness)
em_treat_evenness <- emmeans(glmm_evenness, ~ treatment, type = "response")
em_time_evenness <- emmeans(glmm_evenness, ~ treatment | sampling, type = "response")

## SLA ##
glmm_sla <- glmmTMB(SLA ~ treatment * sampling + ar1(sampling + 0 | plot),
                    data = arkaute_sla, family = gaussian(link = "log"))
#diagnose_glmm(glmm_sla)
em_treat_sla <- emmeans(glmm_sla, ~ treatment, type = "response")
em_time_sla <- emmeans(glmm_sla, ~ treatment | sampling, type = "response")

## LDMC ##
glmm_LDMC <- glmmTMB( LDMC ~ treatment * sampling + ar1(sampling + 0 | plot),
                      data = arkaute_ldmc, family = gaussian())
#diagnose_glmm(glmm_LDMC)
em_treat_ldmc<- emmeans(glmm_LDMC, ~ treatment, type = "response")
em_time_ldmc <- emmeans(glmm_LDMC, ~ treatment | sampling, type = "response")

## Leaf nitrogen ##
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

glmm_dynamics <- do.call(rbind, glmm_time_list) |> 
  common_function() |> 
  mutate(model = paste0("GLMM")) |>  
  select(-contrast)



############# GAM MODELS ###########################################################
samplings_eval <- sort(unique(arkaute_evenness$sampling_num))

## Richness ## 
gam_richness <- gam(
  richness ~ treatment + s(sampling_num, by = treatment, k = 10) + s(plot, bs = "re"),
  data = arkaute_richness, family = gaussian())
#diagnose_gam(gam_richness)
em_treat_richness <- emmeans(gam_richness, ~ treatment)
em_time_richness  <- emmeans(gam_richness, ~ treatment | sampling_num, at = list(sampling_num = samplings_eval))

## Abundance ##
gam_abundance <- gam(abundance ~ treatment + s(sampling_num, by = treatment, k = 10) + s(plot, bs = "re"),
                     data = arkaute_abundance, family = gaussian(link = "identity"))
#diagnose_gam(gam_abundance)
em_treat_abundance <- emmeans(gam_abundance, ~ treatment)
em_time_abundance  <- emmeans(gam_abundance, ~ treatment | sampling_num, at = list(sampling_num = samplings_eval))

## Evenness ##
gam_evenness <- gam(Y_zipf ~ treatment + s(sampling_num, by = treatment, k = 10) + s(plot, bs = "re"),
                    data = arkaute_evenness, family = gaussian())
#diagnose_gam(gam_evenness)
em_treat_evenness <- emmeans(gam_evenness, ~ treatment)
em_time_evenness  <- emmeans(gam_evenness, ~ treatment | sampling_num, at = list(sampling_num = samplings_eval))

## SLA ##
gam_sla <- gam(SLA ~ treatment + s(sampling_num, by = treatment, k = 10) + s(plot, bs = "re"),
               data = arkaute_sla, family = gaussian(link = "log"))
#diagnose_gam(gam_sla)
em_treat_sla <- emmeans(gam_sla, ~ treatment, type = "response")
em_time_sla  <- emmeans(gam_sla, ~ treatment | sampling_num, at = list(sampling_num = samplings_eval),
                        type = "response")

## LDMC ##
gam_ldmc <- gam(LDMC ~ treatment + s(sampling_num, by = treatment, k = 10) + s(plot, bs = "re"),
                data = arkaute_ldmc, family = gaussian())
#diagnose_gam(gam_ldmc)
em_treat_ldmc <- emmeans(gam_ldmc, ~ treatment)
em_time_ldmc  <- emmeans(gam_ldmc, ~ treatment | sampling_num, at = list(sampling_num = samplings_eval))

## Leaf Nitrogen ##
gam_leafN <- gam(leafN ~ treatment + s(sampling_num, by = treatment, k = 10) + s(plot, bs = "re"),
                 data = arkaute_leafN, family = gaussian(link = "identity"))
#diagnose_gam(gam_leafN)
em_treat_leafN <- emmeans(gam_leafN, ~ treatment)
em_time_leafN  <- emmeans(gam_leafN, ~ treatment | sampling_num, at = list(sampling_num = samplings_eval))

## Biomass ##
gam_biomass <- gam(biomass_mice_lm ~ treatment + s(sampling_num, by = treatment, k = 10) + s(plot, bs = "re"),
                   data = arkaute_biomass, family = tw(link = "log"))
#diagnose_gam(gam_biomass)
em_treat_biomass <-  emmeans(gam_biomass, ~ treatment, type = "response")
em_time_biomass  <- emmeans(gam_biomass, ~ treatment | sampling_num, at = list(sampling_num = samplings_eval),
                            type = "response")

###############


gam_list <- list()

gam_list[[1]] <- as.data.frame(pairs(em_treat_richness, adjust = "tukey")) |>
  mutate(variable = "richness", AIC = AIC(gam_richness), estimate_type = "substract")

gam_list[[2]] <- as.data.frame(pairs(em_treat_abundance, adjust = "tukey")) |>
  mutate(variable = "abundance", AIC = AIC(gam_abundance), estimate_type = "substract")

gam_list[[3]] <- as.data.frame(pairs(em_treat_evenness, adjust = "tukey")) |>
  mutate(variable = "evenness", AIC = AIC(gam_evenness), estimate_type = "substract")

gam_list[[4]] <- as.data.frame(pairs(em_treat_sla, adjust = "tukey")) |>
  mutate(variable = "SLA", AIC = AIC(gam_sla), estimate_type = "ratio") |>
  rename(estimate = ratio) |> 
  select(-null)

gam_list[[5]] <- as.data.frame(pairs(em_treat_ldmc, adjust = "tukey")) |>
  mutate(variable = "LDMC", AIC = AIC(gam_ldmc), estimate_type = "substract")

gam_list[[6]] <- as.data.frame(pairs(em_treat_leafN, adjust = "tukey")) |>
  mutate(variable = "leafN", AIC = AIC(gam_leafN), estimate_type = "substract")

gam_list[[7]] <- as.data.frame(pairs(em_treat_biomass, adjust = "tukey")) |>
  mutate(variable = "biomass", AIC = AIC(gam_biomass), estimate_type = "ratio") |>
  rename(estimate = ratio) |> 
  select(-null)

gam_results <- do.call(rbind, gam_list) |>
  common_function() |> 
  mutate(model = "GAM") |>
  select(-contrast)


gam_time_result <- list()

gam_time_result[[1]] <- as.data.frame(pairs(em_time_richness, adjust = "tukey")) |>
  mutate(variable = "richness", AIC = AIC(gam_richness), estimate_type = "substract")

gam_time_result[[2]] <- as.data.frame(pairs(em_time_abundance, adjust = "tukey")) |>
  mutate(variable = "abundance", AIC = AIC(gam_abundance), estimate_type = "substract")

gam_time_result[[3]] <- as.data.frame(pairs(em_time_evenness, adjust = "tukey")) |>
  mutate(variable = "evenness", AIC = AIC(gam_evenness), estimate_type = "substract")

gam_time_result[[4]] <- as.data.frame(pairs(em_time_sla, adjust = "tukey")) |>
  mutate(variable = "SLA", AIC = AIC(gam_sla), estimate_type = "ratio")|>
  rename(estimate = ratio) |> 
  select(-null)

gam_time_result[[5]] <- as.data.frame(pairs(em_time_ldmc, adjust = "tukey")) |>
  mutate(variable = "LDMC", AIC = AIC(gam_ldmc), estimate_type = "substract")

gam_time_result[[6]] <- as.data.frame(pairs(em_time_leafN, adjust = "tukey")) |>
  mutate(variable = "leafN", AIC = AIC(gam_leafN), estimate_type = "substract")

gam_time_result[[7]] <- as.data.frame(pairs(em_time_biomass, adjust = "tukey")) |>
  mutate(variable = "biomass", AIC = AIC(gam_biomass), estimate_type = "ratio") |>
  rename(estimate = ratio) |> 
  select(-null)

gam_dynamics <- do.call(rbind, gam_time_result) |>
  common_function() |> 
  mutate(model = "GAM") |>
  rename(sampling = sampling_num) |> 
  select(-contrast)




####### Joining GLMM, GAM and LRR ############


model_agg <- full_join(glmm_results, gam_results) |> 
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


}


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
  labs(title = "Recovery(p/c)",
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
