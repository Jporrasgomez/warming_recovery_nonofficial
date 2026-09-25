

rm(list = ls(all.names = TRUE))
pacman::p_load(
  dplyr, reshape2, tidyverse, lubridate, ggplot2, ggpubr, gridExtra, ggrepel,
  car, ggsignif, dunn.test, rstatix, performance, emmeans, DHARMa, mgcv, glmmTMB
)

source("code/palettes_labels.R")

## Opening LOG RESPONSE RATIO results datasets #############################

lrr_reading <- function(data) {
  data |> 
    select(-scale, -X) |> 
    rename(estimate = eff_value) |> 
    mutate(model = paste0("LRR"),
           variable = as.factor(variable)) |> 
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
    select(-null_effect)
}


lrr_table <- read.csv("results/effect_size_aggregated.csv") |> 
  lrr_reading()
lrr_table_dyn <- read.csv("results/effect_size_dynamics.csv") |> 
  filter(sampling != "0") |> 
  mutate(sampling = as.factor(sampling)) |> 
  lrr_reading()


## Opening GLMM results datasets ###############################################

glmm_treatment <-  read.csv("results/GLMM_agg.csv") |>  select(-X)
glmm_dynamics  <-  read.csv("results/GLMM_dyn.csv") |>  select(-X) |>  mutate(sampling = as.factor(sampling))

####### Joining GLMM, GAM and LRR in specific way for comparison visualization ############

model_agg <- full_join(glmm_treatment,lrr_table ) |> 
  mutate(variable.bis = variable) |> 
  select(variable, eff_descriptor, upper_limit, lower_limit, model, estimate, z.ratio, p_value, effect_significance,
         effect_sign, variable.bis)


model_dyn <- full_join(glmm_dynamics, lrr_table_dyn) |>
  mutate(
    variable.bis = variable,
    variable_model = paste0(variable, "-", model),
    sampling = fct_reorder(sampling, as.numeric(sampling))
  ) |>
  select(
    variable, eff_descriptor, sampling, eff_descriptor, upper_limit, lower_limit, model, estimate, z.ratio, p_value, effect_significance,
    effect_sign, variable.bis, variable_model)

################### RESULTS COMPARISON VISUALIZATION ###########################

gg_wpp_time<- 
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
       shape = "Sign of effect", color = "Effect significance")
print(gg_wpp_time)

gg_pc_time <- 
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
       shape = "Sign of effect", color = "Effect significance")
print(gg_pc_time)


gg_wc_time <- 
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
       shape = "Sign of effect", color = "Effect significance")
print(gg_wc_time)

gg_wpc_time <-
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
       shape = "Sign of effect", color = "Effect significance")
print(gg_wpc_time)



#ggsave("results/model_comparison_wp_vs_p.png", plot = gg_wpp_time, dpi = 600)
#ggsave("results/model_comparison_p_vs_c.png",  plot = gg_pc_time, dpi = 600)
#ggsave("results/model_comparison_w_vs_c.png",  plot = gg_wc_time, dpi = 600)
#ggsave("results/model_comparison_wp_vs_c.png", plot = gg_wpc_time, dpi = 600)







#### TIME SERIES COMPARISON ###

gg_wpp_dyn <- 
  model_dyn |> 
  filter(eff_descriptor == "wp_vs_p") |> 
  ggplot(aes(x = sampling, y = variable_model, color = effect_significance, shape = effect_sign)) +
  geom_point(size = 8) +
  geom_hline(yintercept = seq(2.5, 16.5, by = 2), linetype = "solid", linewidth = 0.1) +
  scale_color_manual(values = palette_significance) +
  scale_shape_manual(values = shape_significance) +
  labs(title = "Warming effect on post-disturbance dynamics (WP vs. P)",
       subtitle = "Model: variable ~ treatment * sampling + ar(sampling) + (1 | plot)",
       x = "Sampling", 
       y = "Variable and model",
       color = "Effect significance",
       shape = "Sign of effect")
print(gg_wpp_dyn)




gg_pc_dyn <- 
  model_dyn |> 
  filter(eff_descriptor == "p_vs_c") |> 
  ggplot(aes(x = sampling, y = variable_model, color = effect_significance, shape = effect_sign)) +
  geom_point(size = 8) +
  geom_hline(yintercept = seq(2.5, 16.5, by = 2), linetype = "solid", linewidth = 0.1) +
  scale_color_manual(values = palette_significance) +
  scale_shape_manual(values = shape_significance) + 
  labs(title = "Recovery (P vs. C)",
       subtitle = "Model: variable ~ treatment * sampling + ar(sampling) + (1 | plot)",
       x = "Sampling", 
       y = "Variable and model",
       color = "Effect significance",
       shape = "Sign of effect")
print(gg_pc_dyn)



gg_wc_dyn <- 
  model_dyn |> 
  filter(eff_descriptor == "w_vs_c") |> 
  ggplot(aes(x = sampling, y = variable_model, color = effect_significance, shape = effect_sign)) +
  geom_point(size = 8) +
  geom_hline(yintercept = seq(2.5, 16.5, by = 2), linetype = "solid", linewidth = 0.1) +
  scale_color_manual(values = palette_significance) +
  scale_shape_manual(values = shape_significance) +
  labs(title = "Warming effect(W vs. C)",
       subtitle = "Model: variable ~ treatment * sampling + ar(sampling) + (1 | plot)",
       x = "Sampling", 
       y = "Variable and model",
       color = "Effect significance",
       shape = "Sign of effect")
print(gg_wc_dyn)


gg_wpc_dyn <- 
  model_dyn |> 
  filter(eff_descriptor == "wp_vs_c") |> 
  ggplot(aes(x = sampling, y = variable_model, color = effect_significance, shape = effect_sign)) +
  geom_point(size = 8) +
  geom_hline(yintercept = seq(3.5, 18.5, by = 3), linetype = "solid", linewidth = 0.1) +
  scale_color_manual(values = palette_significance) +
  scale_shape_manual(values = shape_significance) +
  labs(title = "Combined effect(wp/c)", 
       subtitle = "Model: variable ~ treatment * sampling + ar(sampling) + (1 | plot)",
       x = "Sampling", 
       y = "Variable and model",
       color = "Effect significance",
       shape = "Sign of effect")
print(gg_wpc_dyn)



#ggsave("results/model_comparison_dynamics_wp_vs_p.png", plot = gg_wpp_dyn, dpi = 600)
#ggsave("results/model_comparison_dynamics_w_vs_c.png", plot = gg_pc_dyn, dpi = 600)
#ggsave("results/model_comparison_dynamics_p_vs_c.png", plot = gg_wc_dyn, dpi = 600)
#ggsave("results/model_comparison_dynamics_wp_vs_c.png", plot = gg_wpc_dyn, dpi = 600)



