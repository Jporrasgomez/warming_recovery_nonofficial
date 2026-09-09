

rm(list = ls(all.names = TRUE))
pacman::p_load(
  dplyr, reshape2, tidyverse, lubridate, ggplot2, ggpubr, gridExtra, ggrepel, patchwork
)

source("code/palettes_labels.R")


notime         <- read.csv("results/model_comparison_notime.csv")
time_treatment <- read.csv("results/model_comparison_time_treatment.csv")
time_dynamics  <- read.csv("results/model_comparison_time_dynamics.csv")


### TREATMENTS  ###

# Simple model

gg_wpp <- 
  notime |> 
  filter(eff_descriptor == "wp_vs_p") |> 
  ggplot(aes(x = model, y = variable, color = effect_significance, shape = effect_sign,
             label = round(p_value, 2))) +
  geom_point(size = 10) +
  geom_label_repel(show.legend = FALSE) +
  scale_color_manual(values = palette_significance) +
  scale_shape_manual(values = shape_significance) +
  labs(title = "Warming effect on recovery (wp vs p)",
       subtitle = "Model: variable ~ treatment + (1 | plot)",
       x = "Statistical analysis", y = NULL,
       shape = "Sign of effect", color = "Effect significance")
print(gg_wpp)

gg_pc <- 
  notime |> 
  filter(eff_descriptor == "p_vs_c") |> 
  ggplot(aes(x = model, y = variable, color = effect_significance, shape = effect_sign,
             label = round(p_value, 2))) +
  geom_point(size = 10) +
  geom_label_repel(show.legend = FALSE) +
  scale_color_manual(values = palette_significance) +
  scale_shape_manual(values = shape_significance) +
  labs(title = "Perturbation effect (p vs c)",
       subtitle = "Model: variable ~ treatment + (1 | plot)",
       x = "Statistical analysis", y = NULL,
       shape = "Sign of effect", color = "Effect significance")
print(gg_pc)

gg_wpc <- 
  notime |> 
  filter(eff_descriptor == "wp_vs_c") |> 
  ggplot(aes(x = model, y = variable, color = effect_significance, shape = effect_sign,
             label = round(p_value, 2))) +
  geom_point(size = 10) +
  geom_label_repel(show.legend = FALSE) +
  scale_color_manual(values = palette_significance) +
  scale_shape_manual(values = shape_significance) +
  labs(title = "Combined effect (wp vs c)",
       subtitle = "Model: variable ~ treatment + (1 | plot)",
       x = "Statistical analysis", y = NULL,
       shape = "Sign of effect", color = "Effect significance")
print(gg_wpc)


gg_wc <- 
  notime |> 
  filter(eff_descriptor == "w_vs_c") |> 
  ggplot(aes(x = model, y = variable, color = effect_significance, shape = effect_sign,
             label = round(p_value, 2))) +
  geom_point(size = 10) +
  geom_label_repel(show.legend = FALSE) +
  scale_color_manual(values = palette_significance) +
  scale_shape_manual(values = shape_significance) +
  labs(title = "Warming effect on assembly (w vs c)",
       subtitle = "Model: variable ~ treatment + ( 1 | plot)",
       x = "Statistical analysis", y = NULL,
       shape = "Sign of effect", color = "Effect significance")
print(gg_wc)


#ggsave("results/model_comparison_wp_vs_p.png", plot = gg_wpp, dpi = 600)
#ggsave("results/model_comparison_p_vs_c.png", plot  = gg_pc, dpi = 600)
#ggsave("results/model_comparison_w_vs_c.png", plot  = gg_wc, dpi = 600)
#ggsave("results/model_comparison_wp_vs_c.png", plot = gg_wpc, dpi = 600)




# Complex model




gg_wpp_time<- 
  time_treatment |> 
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
  time_treatment |> 
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
  time_treatment |> 
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
  time_treatment |> 
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


gg_wpp_final <- 
gg_wpp + gg_wpp_time +
     plot_layout(guides = "collect",
                 widths = c(1, 1)) +
  plot_annotation(tag_levels = "A",
                  theme = theme(legend.position = "bottom"))



gg_pc_final <- 
gg_pc + gg_pc_time + theme(legend.position = "none") + 
  plot_layout(guides = "collect",
              widths = c(1, 1)) +
  plot_annotation(tag_levels = "A",
                  theme = theme(legend.position = "bottom"))



gg_wc_final <- 
gg_wc + gg_wc_time + theme(legend.position = "none") + 
  plot_layout(guides = "collect",
              widths = c(1, 1)) +
  plot_annotation(tag_levels = "A",
                  theme = theme(legend.position = "bottom"))



gg_wpc_final <- 
gg_wpc + gg_wpc_time + theme(legend.position = "none") + 
  plot_layout(guides = "collect",
              widths = c(1, 1)) +
  plot_annotation(tag_levels = "A",
                  theme = theme(legend.position = "bottom"))




ggsave("results/model_comparison_wp_vs_p.png", plot = gg_wpp_final, dpi = 600)
ggsave("results/model_comparison_p_vs_c.png", plot  = gg_pc_final, dpi = 600)
ggsave("results/model_comparison_w_vs_c.png", plot  = gg_wc_final, dpi = 600)
ggsave("results/model_comparison_wp_vs_c.png", plot = gg_wpc_final, dpi = 600)








#### TIME SERIES COMPARISON ###

gg_wpp_dyn <- 
  time_dynamics |> 
  filter(eff_descriptor == "wp_vs_p") |> 
  ggplot(aes(x = sampling, y = variable_model, color = effect_significance, shape = effect_sign)) +
  geom_point(size = 8) +
  geom_hline(yintercept = seq(3.5, 18.5, by = 3), linetype = "solid", linewidth = 0.1) +
  scale_color_manual(values = palette_significance) +
  scale_shape_manual(values = shape_significance) +
  labs(title = "Warming effect on recovery (wp/p)",
       subtitle = "Model: variable ~ treatment * sampling + ar(sampling) + (1 | plot)",
       x = "Sampling", 
       y = "Variable and model",
       color = "Effect significance",
       shape = "Sign of effect")
print(gg_wpp_dyn)




gg_pc_dyn <- 
  time_dynamics |> 
  filter(eff_descriptor == "p_vs_c") |> 
  ggplot(aes(x = sampling, y = variable_model, color = effect_significance, shape = effect_sign)) +
  geom_point(size = 8) +
  geom_hline(yintercept = seq(3.5, 18.5, by = 3), linetype = "solid", linewidth = 0.1) +
  scale_color_manual(values = palette_significance) +
  scale_shape_manual(values = shape_significance) + 
  labs(title = "Recovery (p/c)",
       subtitle = "Model: variable ~ treatment * sampling + ar(sampling) + (1 | plot)",
       x = "Sampling", 
       y = "Variable and model",
       color = "Effect significance",
       shape = "Sign of effect")
print(gg_pc_dyn)



gg_wc_dyn <- 
  time_dynamics |> 
  filter(eff_descriptor == "w_vs_c") |> 
  ggplot(aes(x = sampling, y = variable_model, color = effect_significance, shape = effect_sign)) +
  geom_point(size = 8) +
  geom_hline(yintercept = seq(3.5, 18.5, by = 3), linetype = "solid", linewidth = 0.5) +
  scale_color_manual(values = palette_significance) +
  scale_shape_manual(values = shape_significance) +
  labs(title = "Warming effect(w/c)",
       subtitle = "Model: variable ~ treatment * sampling + ar(sampling) + (1 | plot)",
       x = "Sampling", 
       y = "Variable and model",
       color = "Effect significance",
       shape = "Sign of effect")
print(gg_wc_dyn)


gg_wpc_dyn <- 
  time_dynamics |> 
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



ggsave("results/model_comparison_dynamics_wp_vs_p.png", plot = gg_wpp_dyn, dpi = 600)
ggsave("results/model_comparison_dynamics_w_vs_c.png", plot = gg_pc_dyn, dpi = 600)
ggsave("results/model_comparison_dynamics_p_vs_c.png", plot = gg_wc_dyn, dpi = 600)
ggsave("results/model_comparison_dynamics_wp_vs_c.png", plot = gg_wpc_dyn, dpi = 600)

