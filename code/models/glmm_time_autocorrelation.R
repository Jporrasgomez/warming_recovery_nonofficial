



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
  filter(sampling != "0")


# databases for models:

arkaute_evenness <- arkaute |> filter(sampling != "1", !is.na(Y_zipf))
arkaute_sla <- arkaute |> filter(sampling != "1", !is.na(SLA))
arkaute_ldmc <- arkaute |> filter(sampling != "1", !is.na(LDMC))
arkaute_leafN <- arkaute |> filter(sampling != "1", !is.na(leafN))

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



{
############# GLMM ##############


## Richness ## 

#arkaute_richness <- arkaute |> filter(sampling != "1")
glmm_richness <- glmmTMB(
  richness ~ treatment * sampling + ar1(sampling + 0 | plot),
  data = arkaute,
  family = gaussian
)

diagnose_glmm(glmm_richness)

#General effect of the treatment
em_treat_richness <- emmeans(glmm_richness, ~ treatment, type = "response")
pairs(em_treat_richness, adjust = "tukey")
#Specific contrast wp-p (direct contrast)
contrast(em_treat_richness, method = list("p vs wp" = c(c = 0, p = 1, w = 0, wp = -1)))
as.data.frame(pairs(em_treat_richness, adjust = "tukey"))
# Effect of treatment*samplings on the variable
Anova(glmm_richness, type = "II")

# Pairwise comparisons by sampling event
em_time_richness <- emmeans(glmm_richness, ~ treatment | sampling, type = "response")
time_pairs_richness <- pairs(em_time_richness, adjust = "tukey")
time_pairs_richness

#Specific pairwise compaisons p-wp across samplings
subset(as.data.frame(time_pairs_richness), contrast == "p - wp")

# Specific contrast p-wp 
contrast(
  em_time_richness, 
  method = list("p_vs_wp" = c(c = 0, p = 1, w = 0, wp = -1)), 
  by = "sampling",
  adjust = "none" # Use "none" or "fdr" across the 20 sampling events
)


##### ABUNDANCE #####
# We use gaussian because it tolerates 0 (present in sampling 1 for treatments p and wp)
hist(arkaute$abundance)
arkaute_abundance <- arkaute |> filter(sampling != "1")
glmm_abundance <- glmmTMB(
  abundance ~ treatment * sampling + ar1(sampling + 0 | plot),
  data = arkaute,
  family = gaussian(link = "identity")
)

diagnose_glmm(glmm_abundance)

#General effect of the treatment
em_treat_abundance <- emmeans(glmm_abundance, ~ treatment, type = "response")
pairs(em_treat_abundance, adjust = "tukey")
#Specific contrast wp-p (direct contrast)
contrast(em_treat_abundance, method = list("p vs wp" = c(c = 0, p = 1, w = 0, wp = -1)))

# Effect of treatment*samplings on the variable
Anova(glmm_abundance, type = "II")

# Pairwise comparisons by sampling event
em_time_abundance <- emmeans(glmm_abundance, ~ treatment | sampling, type = "response")
time_pairs_abundance <- pairs(em_time_abundance, adjust = "tukey")
time_pairs_abundance

#Specific pairwise compaisons p-wp across samplings
subset(as.data.frame(time_pairs_abundance), contrast == "p - wp")

# Specific contrast p-wp 
contrast(
  em_time_abundance, 
  method = list("p_vs_wp" = c(c = 0, p = 1, w = 0, wp = -1)), 
  by = "sampling",
  adjust = "none" # Use "none" or "fdr" across the 20 sampling events
)





##### EVENNESS #####
# We choose gaussian() family because Y_zipf is a continuos variable
# of real numbers, and gaussian accepts negative values. Besides, 
# Y_zipf present a (kind of) symmetrical distribution
hist(arkaute$Y_zipf)

glmm_evenness <- glmmTMB(
  Y_zipf ~ treatment * sampling + ar1(sampling + 0 | plot),
  data = arkaute_evenness,
  family = gaussian()
)

diagnose_glmm(glmm_evenness)

#General effect of the treatment
em_treat_evenness <- emmeans(glmm_evenness, ~ treatment, type = "response")
pairs(em_treat_evenness, adjust = "tukey")
#Specific contrast wp-p (direct contrast)
contrast(em_treat_evenness, method = list("p - wp" = c(c = 0, p = 1, w = 0, wp = -1)))


# Effect of treatment*samplings on the variable
Anova(glmm_evenness, type = "II")

# Pairwise comparisons by sampling event
em_time_evenness <- emmeans(glmm_evenness, ~ treatment | sampling, type = "response")
time_pairs_evenness <- pairs(em_time_evenness, adjust = "tukey")
time_pairs_evenness

#Specific pairwise compaisons p-wp across samplings
subset(as.data.frame(time_pairs_evenness), contrast == "p - wp")

# Specific contrast p-wp 
contrast(
  em_time_evenness, 
  method = list("p_vs_wp" = c(c = 0, p = 1, w = 0, wp = -1)), 
  by = "sampling",
  adjust = "none" # Use "none" or "fdr" across the 20 sampling events
)







### SLA ###


glmm_sla <- glmmTMB(
  SLA ~ treatment * sampling + ar1(sampling + 0 | plot),
  data = arkaute_sla,
  family = gaussian(link = "log")
)

diagnose_glmm(glmm_sla)

#General effect of the treatment
em_treat_sla <- emmeans(glmm_sla, ~ treatment, type = "response")
pairs(em_treat_sla, adjust = "tukey")
#Specific contrast wp-p (direct contrast)
contrast(em_treat_sla, method = list("p vs wp" = c(c = 0, p = 1, w = 0, wp = -1)))


# Effect of treatment*samplings on the variable
Anova(glmm_sla, type = "II")

# Pairwise comparisons by sampling event
em_time_sla <- emmeans(glmm_sla, ~ treatment | sampling, type = "response")
time_pairs_sla <- pairs(em_time_sla, adjust = "tukey")
time_pairs_sla

#Specific pairwise compaisons p-wp across samplings
subset(as.data.frame(time_pairs_sla), contrast == "p / wp")

# Specific contrast p-wp 
contrast(
  em_time_sla, 
  method = list("p_vs_wp" = c(c = 0, p = 1, w = 0, wp = -1)), 
  by = "sampling",
  adjust = "none" # Use "none" or "fdr" across the 20 sampling events
)







### LDMC ###

hist(arkaute_ldmc$LDMC)
glmm_LDMC <- glmmTMB(
  LDMC ~ treatment * sampling + ar1(sampling + 0 | plot),
  data = arkaute_ldmc,
  family = gaussian()
)

diagnose_glmm(glmm_LDMC)

#General effect of the treatment
em_treat_ldmc<- emmeans(glmm_LDMC, ~ treatment, type = "response")
pairs(em_treat_ldmc, adjust = "tukey")
#Specific contrast wp-p (direct contrast)
contrast(em_treat_ldmc, method = list("p vs wp" = c(c = 0, p = 1, w = 0, wp = -1)))


# Effect of treatment*samplings on the variable
Anova(glmm_LDMC, type = "II")

# Pairwise comparisons by sampling event
em_time_ldmc <- emmeans(glmm_LDMC, ~ treatment | sampling, type = "response")
res_pairs_ldmc <- pairs(em_time_ldmc, adjust = "tukey")
res_pairs_ldmc

#Specific pairwise compaisons p-wp across samplings
subset(as.data.frame(res_pairs_ldmc), contrast == "p - wp")

# Specific contrast p-wp 
contrast(
  em_time_ldmc, 
  method = list("p_vs_wp" = c(c = 0, p = 1, w = 0, wp = -1)), 
  by = "sampling",
  adjust = "none" # Use "none" or "fdr" across the 20 sampling events
)






### Leaf nitrogen ###

hist(arkaute_leafN$leafN)
glmm_leafN <- glmmTMB(
  leafN ~ treatment * sampling + ar1(sampling + 0 | plot),
  data = arkaute_leafN,
  family = gaussian(link = "identity")
)

diagnose_glmm(glmm_leafN)

#General effect of the treatment
em_treat_leafN <- emmeans(glmm_leafN, ~ treatment, type = "response")
pairs(em_treat_leafN, adjust = "tukey")
#Specific contrast wp-p (direct contrast)
contrast(em_treat_leafN, method = list("p vs wp" = c(c = 0, p = 1, w = 0, wp = -1)))


# Effect of treatment*samplings on the variable
Anova(glmm_leafN, type = "II")

# Pairwise comparisons by sampling event
em_time_leafN <- emmeans(glmm_leafN, ~ treatment | sampling, type = "response")
time_pairs_leafN <- pairs(em_time_leafN, adjust = "tukey")
time_pairs_leafN

#Specific pairwise compaisons p-wp across samplings
subset(as.data.frame(time_pairs_leafN), contrast == "p - wp")

# Specific contrast p-wp 
contrast(
  em_time_leafN, 
  method = list("p_vs_wp" = c(c = 0, p = 1, w = 0, wp = -1)), 
  by = "sampling",
  adjust = "none" # Use "none" or "fdr" across the 20 sampling events
)







### BIOMASS ###

arkaute_biomass <- arkaute |> filter(sampling != "1")
hist(arkaute$biomass_mice_lm)
quantile(arkaute$biomass_mice_lm)
min(arkaute$biomass_mice_lm)


arkaute$sampling <- relevel(arkaute$sampling, ref = "2")
glmm_biomass <- glmmTMB(
  biomass_mice_lm ~ treatment * sampling + ar1(sampling + 0 | plot),
  data = arkaute,
  family = tweedie(link = "log")
)

diagnose_glmm(glmm_biomass)

#General effect of the treatment
em_treat_biomass <- emmeans(glmm_biomass, ~ treatment, type = "response")
pairs(em_treat_biomass, adjust = "tukey")
#Specific contrast wp-p (direct contrast)
contrast(em_treat_biomass, method = list("p vs wp" = c(c = 0, p = 1, w = 0, wp = -1)))


# Effect of treatment*samplings on the variable
Anova(glmm_biomass, type = "II")

# Pairwise comparisons by sampling event
em_time_biomass <- emmeans(glmm_biomass, ~ treatment | sampling, type = "response")
time_pairs_biomass <- pairs(em_time_biomass, adjust = "tukey")
time_pairs_biomass

#Specific pairwise compaisons p-wp across samplings
subset(as.data.frame(time_pairs_biomass), contrast == "p - wp")

# Specific contrast p-wp 
contrast(
  em_time_biomass, 
  method = list("p_vs_wp" = c(c = 0, p = 1, w = 0, wp = -1)), 
  by = "sampling",
  adjust = "none" # Use "none" or "fdr" across the 20 sampling events
)


}

## Cambiar nombres de contrastes (p_vs_c, por ejemplo).
## Añadir familias
## Añadir problemas de dispersion y tal
## Juntar con resultados actualiados de log response ratio. 
{
model_list <- list()

model_list[[1]] <- 
  as.data.frame(pairs(em_treat_richness, adjust = "tukey")
                ) |>
  mutate(variable = paste0("richness"),
         AIC = AIC(glmm_richness)
                ) |> 
  rename(estimate_ratio = estimate)

model_list[[2]] <- 
  as.data.frame(pairs(em_treat_abundance, adjust = "tukey")
                ) |> 
  mutate(variable = paste0("abundance"),
         AIC = AIC(glmm_abundance)
         )|> 
  rename(estimate_ratio = estimate)

model_list[[3]]  <- 
  as.data.frame(pairs(em_treat_evenness, adjust = "tukey")
  ) |> mutate(
    variable = paste0("evenness"),
    AIC = AIC(glmm_evenness)
    )|> 
  rename(estimate_ratio = estimate)

model_list[[4]] <- 
  as.data.frame(pairs(em_treat_sla, adjust = "tukey")
  ) |> mutate(
    variable = paste0("SLA"),
    AIC = AIC(glmm_sla)
    )|> 
  rename(estimate_ratio = ratio) |> 
  select(-null)

model_list[[5]]  <- 
  as.data.frame(pairs(em_treat_ldmc, adjust = "tukey")
  ) |> mutate(
    variable = paste0("LDMC"),
    AIC = AIC(glmm_LDMC)
    )|> 
  rename(estimate_ratio = estimate)

model_list[[6]]  <- 
  as.data.frame(pairs(em_treat_leafN, adjust = "tukey")
  ) |> 
  mutate(variable = paste0("leafN"),
         AIC = AIC(glmm_leafN)
         )|> 
  rename(estimate_ratio = estimate)

model_list[[7]]  <- 
  as.data.frame(pairs(em_treat_biomass, adjust = "tukey")
  ) |> mutate(
    variable = paste0("biomass"),
    AIC = AIC(glmm_biomass)
    )|> 
  rename(estimate_ratio = ratio) |> 
  select(-null)


model_result <- do.call(rbind, model_list) |> 
  filter(!contrast %in% c("p - w", "p / w", "w - wp", "w /wp")) |> 
  rename(diff_value = estimate_ratio,
         p_value = p.value) |> 
  mutate(
    effect_sign = ifelse(diff_value > 0 , "negative", "positive"), 
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
    ),
    model = paste0("GLMM")
  ) |>  
  select(-contrast)

}

# Actualizar 
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




models <- full_join(model_result, lrr_table) |> 
  mutate(variable.bis = variable) |> 
  select(variable, eff_descriptor, model, AIC, 
          p_value, effect_significance, diff_value,
         effect_sign, variable.bis
  )


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
  #geom_label_repel(show.legend = FALSE) +
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






























model_time_result <- list()

model_time_result[[1]] <- 
  as.data.frame(pairs(em_time_richness, adjust = "tukey")
  ) |> mutate(variable = paste0("richness")
  ) |> 
  rename(estimate_ratio = estimate)
model_time_result[[2]] <- 
  as.data.frame(pairs(em_time_abundance, adjust = "tukey")
  ) |> mutate(variable = paste0("abundance"))|> 
  rename(estimate_ratio = estimate)
model_time_result[[3]]  <- 
  as.data.frame(pairs(em_time_evenness, adjust = "tukey")
  ) |> mutate(variable = paste0("evenness"))|> 
  rename(estimate_ratio = estimate)
model_time_result[[4]] <- 
  as.data.frame(pairs(em_time_sla, adjust = "tukey")
  ) |> mutate(variable = paste0("SLA"))|> 
  rename(estimate_ratio = ratio) |> 
  select(-null)
model_time_result[[5]]  <- 
  as.data.frame(pairs(em_time_ldmc, adjust = "tukey")
  ) |> mutate(variable = paste0("LDMC"))|> 
  rename(estimate_ratio = estimate)
model_time_result[[6]]  <- 
  as.data.frame(pairs(em_time_leafN, adjust = "tukey")
  ) |> mutate(variable = paste0("leafN"))|> 
  rename(estimate_ratio = estimate)
model_time_result[[7]]  <- 
  as.data.frame(pairs(em_time_biomass, adjust = "tukey")
  ) |> mutate(variable = paste0("biomass"))|> 
  rename(estimate_ratio = ratio) |> 
  select(-null)

model_time <- do.call(rbind, model_time_result) |> 
  filter(!contrast %in% c("p - w", "p / w", "w - wp", "w /wp")) |> 
  rename(diff_value = estimate_ratio,
         p_value = p.value) |> 
  mutate(
    effect_sign = ifelse(diff_value > 0 , "negative", "positive"), 
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
    ),
    model = paste0("GLMM")
  ) |>  
  select(-contrast)





lrr_table_dyn <- read.csv("results/effect_size_dynamics.csv") |> 
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





models_dynamics <- full_join(model_time, lrr_table_dyn) |> 
  mutate(variable.bis = variable) |> 
  select(variable, sampling, eff_descriptor, model, 
         p_value, effect_significance, diff_value,
         effect_sign, variable.bis
  ) |> 
  mutate(
    variable_model = paste0(variable, "-", model)
  )


# Juntar en una columna variable-model para usar en el eje Y y así verlo todo junto
models_dynamics |> 
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


models_dynamics |> 
  filter(eff_descriptor == "p_vs_c") |> 
  ggplot(aes(x = sampling, y = variable_model, color = effect_significance, shape = effect_sign)) +
  #facet_wrap(~model) + 
  geom_point(size = 8) +
  scale_color_manual(values = palette_sig) +
  scale_shape_manual(values = palette_shape) +
  labs(title = "Recovery(p/c)", 
       x = "Sampling", 
       y = "Variable and model")



models_dynamics |> 
  filter(eff_descriptor == "w_vs_c") |> 
  ggplot(aes(x = sampling, y = variable_model, color = effect_significance, shape = effect_sign)) +
  #facet_wrap(~model) + 
  geom_point(size = 8) +
  scale_color_manual(values = palette_sig) +
  scale_shape_manual(values = palette_shape) +
  labs(title = "Warming effect(w/c)", 
       x = "Sampling", 
       y = "Variable and model")



models_dynamics |> 
  filter(eff_descriptor == "wp_vs_c") |> 
  ggplot(aes(x = sampling, y = variable_model, color = effect_significance, shape = effect_sign)) +
  #facet_wrap(~model) + 
  geom_point(size = 8) +
  scale_color_manual(values = palette_sig) +
  scale_shape_manual(values = palette_shape) +
  labs(title = "Combined effect(wp/c)", 
       x = "Sampling", 
       y = "Variable and model")

