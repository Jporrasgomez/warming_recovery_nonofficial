



rm(list = ls(all.names = TRUE))
pacman::p_load(
  dplyr, reshape2, tidyverse, lubridate, ggplot2, ggpubr, gridExtra, ggrepel,
  car, ggsignif, dunn.test, rstatix, performance, emmeans, DHARMa, mgcv, glmmTMB, patchwork
)

source("code/palettes_labels.R")

{

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
  arrange(plot, sampling)



arkaute_richness  <- arkaute |> filter(!is.na(richness))
arkaute_abundance <- arkaute |> filter(!is.na(abundance))
arkaute_evenness  <- arkaute |> filter(!is.na(Y_zipf))
arkaute_sla       <- arkaute |> filter(!is.na(SLA))
arkaute_ldmc      <- arkaute |> filter(!is.na(LDMC))
arkaute_leafN     <- arkaute |> filter(!is.na(leafN))
arkaute_biomass   <- arkaute |> filter(!is.na(biomass_mice_lm))
arkaute_biomass_mice <- arkaute |> filter(!is.na(biomass_mice))
arkaute_biomass_raw <- arkaute |> filter(!is.na(biomass_raw), biomass_raw != 0 )

data_list <- list(arkaute_richness, arkaute_abundance, arkaute_evenness, arkaute_sla, 
                  arkaute_ldmc, arkaute_leafN, arkaute_biomass, arkaute_biomass_mice, arkaute_biomass_raw)


################# UNDERSTANDING INTERACTION TEST BETWEEN FACTORS #####################


#GLMMs following a factorial design 2 x 2 (OTC / Ambient  and Tilled / Undisturbed) ###############

#     2 x 2 factorial       Treatment
# OTC contains the levels "OTC" and "Ambient"
# TILLING contains the levels "Tilled" and "Undisturbed"

# Ambient, Undisturbed == CONTROL
# Ambient, Tilled      == PERTURBED-ONLY
# OTC, Undisturbed     == WARMED-ONLY
# OTC, Tilled          == COMBINED


#  I use "richness" as an example


# 1. Fit the GLMM with the full three-way interaction: OTC * TILLING * sampling
# This model evaluates whether the effect of warming on disturbed plots differs from undisturbed plots,
# and whether this interaction changes over time.
glmm_richness <- glmmTMB(richness ~ OTC * TILLING * sampling + (1 | plot),
                          data = arkaute_richness,family = genpois())


# 2. Calculate Estimated Marginal Means (EMMEANS) for the global 2x2 factorial design
# NOTE: Results are averaged over the level of 'sampling' (i.e., temporal variability is collapsed)

em_treat_richness <- emmeans(glmm_richness, ~ OTC * TILLING, type = "response")

# 3. Simple pairwise comparisons between treatment combinations (Averaged over time)
# NOTE: These represent the "simple main effects" (e.g., WP vs P, or W vs C)
as.data.frame(pairs(em_treat_richness, adjust = "tukey"))  
# This is the comparison between treatments (same thing we do with variable ~ treatment * sampling + (1|plot))


# 4. The Global Interaction Contrast (Comparison of Comparisons / Difference of Differences)
# NOTE: This provides the "direct statistical assessment of the warming × disturbance interaction",
# evaluating if the effect of warming under disturbance differs from ambient.
# Since it is a genpois family (log link) and we use type = "response", the contrast represents the ratio: (P/WP) / (C/W).
# Note that due to alphabetical factor ordering by default (e.g., 'Ambient' before 'OTC'), 
# the internal calculation evaluates:
# Interaction = [Ln(Ambient, Tilled) - Ln(OTC, Tilled)] - [Ln(Ambient, Undisturbed) - Ln(OTC, Undisturbed)]
# What is the same that <- Interaction = (P-WP) - (C-W)

richness_interaction_contrast <- contrast(em_treat_richness, interaction = "pairwise")
summary(richness_interaction_contrast, type = "response")



# 5. Type III ANOVA Table (Evaluating Temporal Dynamics)
# NOTE: The term 'OTC:TILLING:sampling' provides the formal test for the "temporal variation" 
# of the interaction. It evaluates whether the warming × disturbance 
# interaction changes significantly across the temporal trajectories.
# Under standard dummy coding ('contr.treatment'), the two-way term 'OTC:TILLING' tests the interaction 
# ONLY at the reference sampling date, not the global average calculated in step 4.
Anova(glmm_richness, type = "III")


# 6. Time-Conditioned EMMEANS (Slicing the 3-way interaction by sampling)
# Same result  as a model with treatment as factor
# variable ~ treatment * sampling + (1 | plot)
em_time_richness <- emmeans(glmm_richness, ~ OTC * TILLING | sampling, type = "response")
as.data.frame(pairs(em_time_richness, adjust = "tukey")) 


# 7. Time-Conditioned Interaction Contrasts
# This calculates the interaction contrast (difference of differences) as in step 4, but independently for every sampling event.
# It mathematically demonstrates the "temporal variation" by pinpointing exactly during which specific sampling 
# dates (e.g., the early recolonization phase) the warming treatment exerted a significantly different effect 
# on perturbed plots compared to undisturbed plots.

richness_time_interaction_contrasts <- contrast(em_time_richness, interaction = "pairwise")
summary(richness_time_interaction_contrasts, type = "response")






############ APPLYING FOR EVERY VARIABLE #################################################


variables <- c("richness", "abundance", "evenness", "SLA", "LDMC", "LeafN", 
               "biomass", "biomass_mice", "biomass_raw")


############# GLMM MODELS ##############################################################
## Richness ## 
glmm_richness <- glmmTMB( richness ~ OTC * TILLING * sampling + (1 | plot),
                          data = arkaute_richness,family = genpois(),
                          contrasts = list(OTC = "contr.sum", TILLING = "contr.sum", sampling = "contr.sum"))
#contrasts = list(OTC = "contr.sum", TILLING = "contr.sum", sampling = "contr.sum")
# if we add this line  in the glmm it makes ANOVA Type III OTC:TILLING to be inferred globally, otherwise it uses only the first sampling


## ABUNDANCE ##
# We use gaussian because it tolerates 0 (present in sampling 1 for treatments p and wp)
glmm_abundance <- glmmTMB(abundance ~ OTC * TILLING * sampling + (1 | plot),
                          data = arkaute_abundance, family = gaussian(),
                          contrasts = list(OTC = "contr.sum", TILLING = "contr.sum", sampling = "contr.sum"))
## EVENNESS ##
glmm_evenness <- glmmTMB(Y_zipf ~ OTC * TILLING * sampling + (1 | plot),
                         dispformula = ~treatment + sampling, 
                         data = arkaute_evenness, family = gaussian(),
                         contrasts = list(OTC = "contr.sum", TILLING = "contr.sum", sampling = "contr.sum"))
## SLA ##
glmm_sla <- glmmTMB(SLA ~ OTC * TILLING * sampling + (1 | plot), 
                    dispformula = ~ treatment, 
                    data = arkaute_sla, family = gaussian(),
                    contrasts = list(OTC = "contr.sum", TILLING = "contr.sum", sampling = "contr.sum"))
## LDMC ##
glmm_LDMC <- glmmTMB(LDMC ~ OTC * TILLING * sampling + (1 | plot), 
                     dispformula = ~ treatment, 
                     data = arkaute_ldmc, family = gaussian(),
                     contrasts = list(OTC = "contr.sum", TILLING = "contr.sum", sampling = "contr.sum"))
## Leaf nitrogen ##
glmm_leafN <-  glmmTMB(leafN ~ OTC * TILLING * sampling + (1 | plot),
                       dispformula = ~ treatment,
                       data = arkaute_leafN, family = gaussian(),
                       contrasts = list(OTC = "contr.sum", TILLING = "contr.sum", sampling = "contr.sum"))
### BIOMASS ###
glmm_biomass <- glmmTMB(biomass_mice_lm ~ OTC * TILLING * sampling  + (1 | plot),
                        dispformula = ~treatment + sampling, 
                        data = arkaute_biomass, family = Gamma(link = "log"),
                        contrasts = list(OTC = "contr.sum", TILLING = "contr.sum", sampling = "contr.sum"))
## BIOMASS - just MICE ####
glmm_biomass_mice <- glmmTMB(biomass_mice ~ OTC * TILLING * sampling +  (1 | plot),
                             dispformula = ~treatment + sampling, 
                             data = arkaute_biomass_mice,  family = Gamma(link = "log"),
                             contrasts = list(OTC = "contr.sum", TILLING = "contr.sum", sampling = "contr.sum"))
## BIOMASS raw ##
glmm_biomass_raw <- glmmTMB(biomass_raw ~ OTC * TILLING * sampling + (1 | plot),
                            dispformula = ~treatment + sampling, 
                            data = arkaute_biomass_raw,  family = Gamma(link = "log"),
                            contrasts = list(OTC = "contr.sum", TILLING = "contr.sum", sampling = "contr.sum"))

glmm_list <- list(glmm_richness, glmm_abundance, glmm_evenness, glmm_sla, 
                  glmm_LDMC, glmm_leafN, glmm_biomass, glmm_biomass_mice, 
                  glmm_biomass_raw)




############ LOOP 1:  GLOBAL INTERACTION (PROMEDIADA EN EL TIEMPO) ############

emmeans_agg_list <- list()

for(i in seq_along(glmm_list)) {
  em_agg <- emmeans(glmm_list[[i]], ~ OTC * TILLING, type = "response")
  inter_df <-
    as.data.frame(
      # We add this list to ensure the interaction to be "WP-P" and "W-C", otherwise it would be
      #"P-WP" and "C-W" (due to alphabetical sorting of the models (Ambient before OTC))
      summary(contrast(em_agg, interaction = list(OTC = "revpairwise", TILLING = "pairwise")), 
                                    infer = c(TRUE, TRUE))) |> 
    mutate(variable = variables[i])
  
  if ("ratio" %in% names(inter_df)) {
    emmeans_agg_list[[i]] <- inter_df |> 
      rename(estimate = ratio, lower_CL = asymp.LCL, upper_CL = asymp.UCL) |> 
      mutate(estimate_type = "ratio") |> 
      select(-any_of("null"))
  } else {
    emmeans_agg_list[[i]] <- inter_df |> 
      rename(lower_CL = asymp.LCL, upper_CL = asymp.UCL) |> 
      mutate(estimate_type = "difference") |> 
      select(-any_of("null"))
  }
}

interaction_test_agg <- bind_rows(emmeans_agg_list) |> 
  mutate(
    ecological_effect = case_when(
      estimate_type == "difference" & estimate < 0 ~ "suppression_by_warming",
      estimate_type == "difference" & estimate > 0 ~ "enhancement_by_warming",
      estimate_type == "ratio"      & estimate < 1 ~ "suppression_by_warming",
      estimate_type == "ratio"      & estimate > 1 ~ "enhancement_by_warming"
    ),
    significance = case_when(
      p.value < 0.05 ~ "p < 0.05",
      p.value < 0.10 ~ "0.05 <= p < 0.10",
      TRUE           ~ "Non-significant"
    )
  )

write.csv(interaction_test_agg, "results/Interaction_test_agg.csv", row.names = FALSE)


############ LOOP 2:  ANOVA  ############

anova_list <- list()

for(i in seq_along(glmm_list)) {
  aov <- as.data.frame(Anova(glmm_list[[i]], type = "III"))
  aov$factor <- rownames(aov)
  rownames(aov) <- NULL
  anova_list[[i]] <- aov |> mutate(variable = variables[i])
}

anova_list <- list()

for(i in seq_along(glmm_list)) {
  
  aov_df <- as.data.frame(Anova(glmm_list[[i]], type = "III"))
  
  anova_list[[i]] <- aov_df |> 
    tibble::rownames_to_column(var = "factor") |> 
    mutate(variable = variables[i])
}

interaction_test_anova <- bind_rows(anova_list)

write.csv(interaction_test_anova, "results/Interaction_test_ANOVA.csv", row.names = FALSE)


############ LOOP 3:  TEMPORAL INTERACTION (x SAMPLING) ############

emmeans_time_list <- list()

for(i in seq_along(glmm_list)) {
  em_time <- emmeans(glmm_list[[i]], ~ OTC * TILLING | sampling, type = "response")
  inter_time_df <- 
    as.data.frame(
      # We add this list to ensure the interaction to be "WP-P" and "W-C", otherwise it would be
      #"P-WP" and "C-W" (due to alphabetical sorting of the models (Ambient before OTC))
      summary(contrast(em_time, interaction = list(OTC = "revpairwise", TILLING = "pairwise")), 
                                         infer = c(TRUE, TRUE))) |> 
    mutate(variable = variables[i])
  
  if ("ratio" %in% names(inter_time_df)) {
    emmeans_time_list[[i]] <- inter_time_df |> 
      rename(estimate = ratio, lower_CL = asymp.LCL, upper_CL = asymp.UCL) |> 
      mutate(estimate_type = "ratio") |> 
      select(-any_of("null"))
  } else {
    emmeans_time_list[[i]] <- inter_time_df |> 
      rename(lower_CL = asymp.LCL, upper_CL = asymp.UCL) |> 
      mutate(estimate_type = "difference") |> 
      select(-any_of("null"))
  }
}

interaction_test_time <- bind_rows(emmeans_time_list) |> 
  mutate(
    ecological_effect = case_when(
      estimate_type == "difference" & estimate < 0 ~ "suppression_by_warming",
      estimate_type == "difference" & estimate > 0 ~ "enhancement_by_warming",
      estimate_type == "ratio"      & estimate < 1 ~ "suppression_by_warming",
      estimate_type == "ratio"      & estimate > 1 ~ "enhancement_by_warming"
    ),
    significance = case_when(
      p.value < 0.05 ~ "p < 0.05",
      p.value < 0.10 ~ "0.05 <= p < 0.10",
      TRUE           ~ "Non-significant"
    )
  )


# Integrating sampling dates

sampling_dates <- read.csv("data/sampling_dates.csv") |> 
  mutate(sampling          = as.factor(sampling),
         date              = ymd(date), 
         year              = year(date),
         date_label        = format(date, "%b %d %y"),
         date_label        = factor(date_label, levels = format(sort(unique(date)), "%b %d %y"), ordered = TRUE),
         date_label_noyear = substr(as.character(date_label), 1, nchar(as.character(date_label)) - 3),
         date_label_noyear = factor( date_label_noyear,levels = unique(date_label_noyear[order(date)]),ordered = TRUE)
  ) |> 
  select(sampling, date, year, date_label_noyear) |> 
  mutate(across(where(is.character), as.factor))


interaction_test_time <- interaction_test_time |>  left_join(sampling_dates)

write.csv(interaction_test_time, "results/Interaction_test_sampling.csv", row.names = FALSE)






############ VISUALS ############





labels_main_variables <- c("richness" = "Richness",            # 1
                           "abundance" = "Cover",              # 2     
                           "evenness" = "Evenness",              # 3     
                           "SLA" = "SLA",                      # 4     
                           "LDMC" = "LDMC",                    # 5     
                           "LeafN"= "Leaf - N",                      # 6     
                           "biomass" = "Biomass"       # 7
)    

# When k = 2
labels_biomass_variables <- c("biomass_raw" = "No imputation", "biomass_mice" = "MICE",
                              "biomass" = "MICE + LM")      
sorting_function <- function(data){
  data |> filter(variable %in% lvls) |> 
    mutate(variable = factor(variable, levels = lvls, labels = labs))
}



labels_list <- list(labels_main_variables, labels_biomass_variables)

gg_ratio_list <- list()
gg_diff_list <- list()

for (k in seq_along(labels_list)){
  
lvls <- names(labels_list[[k]])
labs <- unname(labels_list[[k]])

gg_agg_ratio <- interaction_test_agg |> 
  filter(estimate_type == "ratio") |> 
  sorting_function() |>  
  ggplot(aes(x = variable, y = estimate, color = significance, fill = significance)) +
  facet_wrap(~ variable, scales = "free", ncol = 1, strip.position = "left") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_linerange(aes(ymin = lower_CL, ymax = upper_CL), alpha = 1, linewidth = 1) +
  geom_point(aes(shape = ecological_effect), size = 4) +
  scale_y_log10() +
  scale_color_manual(values = palette_significance2) +
  scale_fill_manual(values = palette_significance2) +
  scale_shape_manual(
    values = shape_values,
    labels = shape_labels
  ) +
  labs(
    y = "Interaction Ratio [(WP / P) / (W / C)] (log scale)",
    #y = NULL, 
    x = NULL,
    color = "Significance", fill = "Significance", shape = "Interaction effect"
  ) +
  theme1 +
  theme(
    legend.position = "none", 
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    #strip.text = element_text(color = "white")
    strip.placement    = "outside",
    strip.text         = element_text(face = "bold", size = 12),
    strip.text.y.left  = element_text(face = "bold"),
  )

gg_agg_diff <- interaction_test_agg |> 
  filter(estimate_type == "difference") |> 
  sorting_function() |> 
  ggplot(aes(x = variable, y = estimate, color = significance, fill = significance)) +
  facet_wrap(~ variable, scales = "free", ncol = 1, strip.position = "left") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_linerange(aes(ymin = lower_CL, ymax = upper_CL), alpha = 1, linewidth = 1) +
  geom_point(aes(shape = ecological_effect), size = 4) +
  scale_color_manual(values = palette_significance2) +
  scale_fill_manual(values = palette_significance2) +
  scale_shape_manual(
    values = shape_values,
    labels = shape_labels
  ) +
  labs(
    y = "Interaction Difference [(WP - P) - (W - C)]",
    #y = NULL, 
    x = NULL,
    color = "Significance", fill = "Significance", shape = "Interaction effect"
  ) +
  theme1 +
  theme(
    legend.position = "none", 
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    #strip.text = element_text(color = "white")
    strip.placement    = "outside",
    strip.text         = element_text(face = "bold", size = 12),
    strip.text.y.left  = element_text(face = "bold"),
  )

gg_dyn_ratio <- interaction_test_time |> 
  filter(estimate_type == "ratio") |> 
  sorting_function() |> 
  ggplot(aes(x = date_label_noyear, y = estimate, color = significance, fill = significance)) +
  facet_wrap(~ variable, scales = "free_y", ncol = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 10.5, linetype = "solid", color = "grey50") +
  geom_linerange(aes(ymin = lower_CL, ymax = upper_CL), alpha = 1, linewidth = 0.6) +
  geom_path(group = 1, color = "grey50") + 
  geom_point(aes(shape = ecological_effect), size = 2.5) +
  scale_y_log10() +
  scale_color_manual(values = palette_significance2) +
  scale_fill_manual(values = palette_significance2) +
  scale_shape_manual(
    values = shape_values,
    labels = shape_labels
  ) +
  labs(
    y = NULL, x = NULL,
    color = "Significance", fill = "Significance", shape = "Interaction effect"
  ) +
  theme1 +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_blank(), 
        strip.background = element_blank())


gg_dyn_diff <- interaction_test_time |> 
  sorting_function() |> 
  filter(estimate_type == "difference") |> 
  ggplot(aes(x = date_label_noyear, y = estimate, color = significance, fill = significance)) +
  facet_wrap(~ variable, scales = "free_y", ncol = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 10.5, linetype = "solid", color = "grey50") +
  geom_linerange(aes(ymin = lower_CL, ymax = upper_CL), alpha = 1, linewidth = 0.6) +
  geom_path(group = 1, color = "grey50") + 
  geom_point(aes(shape = ecological_effect), size = 2.5) +
  scale_color_manual(values = palette_significance2) +
  scale_fill_manual(values = palette_significance2) +
  scale_shape_manual(
    values = shape_values,
    labels = shape_labels
  ) +
  labs(
    y = NULL, x = NULL,
    color = "Significance", fill = "Significance", shape = "Interaction effect"
  ) +
  theme1 +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_blank(), 
        strip.background = element_blank())



gg_ratio_list[[k]] <- gg_agg_ratio + 
  gg_dyn_ratio + 
  plot_layout(guides = "collect", widths = c(1, 3)) +
  plot_annotation(theme = theme(legend.position = "bottom"))

gg_diff_list[[k]] <- gg_agg_diff + 
  gg_dyn_diff + 
  plot_layout(guides = "collect", widths = c(1, 3)) +
  plot_annotation(theme = theme(legend.position = "bottom"))

}

}

print(interaction_test_agg)
tail(interaction_test_time)
print(interaction_test_anova)
print(interaction_test_anova |>  filter(factor %in% ("OTC:TILLING:sampling")))
print(interaction_test_time |>  filter(p.value <= 0.1))

# All variables
gg_ratio_list[[1]]
ggsave("results/Interaction_effect_ratio.png", plot = gg_ratio_list[[1]], dpi = 600)
ggsave("results/Interaction_effect_ratio.svg", plot = gg_ratio_list[[1]], dpi = 600)

gg_diff_list[[1]]
ggsave("results/Interaction_effect_difference.png", plot = gg_diff_list[[1]], dpi = 600)
ggsave("results/Interaction_effect_difference.svg", plot = gg_diff_list[[1]], dpi = 600)

# Biomass imputation
gg_ratio_list[[2]]
ggsave("results/Interaction_effect_BIOMASS_ratio.png", plot = gg_ratio_list[[2]], dpi = 600)
ggsave("results/Interaction_effect_BIOMASS_ratio.svg", plot = gg_ratio_list[[2]], dpi = 600)
  