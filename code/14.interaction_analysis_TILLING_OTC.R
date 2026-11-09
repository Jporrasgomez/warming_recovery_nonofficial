



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

# Function to check the need of implementing an ar1 (temporal autocorrelation term)

ar_test <- function(model, data){
  testTemporalAutocorrelation(recalculateResiduals(simulateResiduals(model),
                                                   group = data$sampling_num),
                              time = sort(unique(data$sampling_num)))
  
}# If p-value is higher than 0.05 we do not need autocorrelation term






str(arkaute)

####### GLMMs following a factorial design 2 x 2 (OTC / Ambient  and Tilled / Undisturbed) ###############

#     2 x 2 factorial       Treatment

# Ambient, Undisturbed == CONTROL
# Ambient, Tilled      == PERTURBED-ONLY
# OTC, Undisturbed     == WARMED-ONLY
# OTC, Tilled          == COMBINED


########## RICHNESS ###########


glmm_richness <- glmmTMB(richness ~ OTC * TILLING * sampling + (1 | plot),
                          data = arkaute_richness,family = genpois())

diagnose_glmm(glmm_richness)
ar_test(glmm_richness, arkaute_richness)

# Calculating Estimated Marginal Means (EMMEANS) of factorial design 2 x 2
# Results are averaged over the level of: sampling (i.e. no temporal variability included)
em_treat_richness <- emmeans(glmm_richness, ~ OTC * TILLING, type = "response")
as.data.frame(pairs(em_treat_richness, adjust = "tukey"))  
# This is the comparison between treatments (same thing we do with variable ~ treatment * sampling + (1|plot))

# This is the key difference: COMPARISON OF COMPARISONS (or difference of differences) 
# but WITHOUT SAMPLING inclusion
# Since it is a genpois family, the contrast is calculated as
# Contrast = [Ln(Ambient, Tilled) - Ln(OTC, Tilled)] - [Ln(Ambient, Undisturbed) - Ln(OTC, Undisturbed)]
richness_interaction_contrast <- contrast(em_treat_richness, interaction = "pairwise")
summary(richness_interaction_contrast, type = "response")


# Considering temporal dynamics: 
Anova(glmm_richness, type = "III")

# Same thing as variable ~ treatment * sampling + (1 | plot)
em_time_richness <- emmeans(glmm_richness, ~ OTC * TILLING | sampling, type = "response")
as.data.frame(pairs(em_time_richness, adjust = "tukey")) 


# COMPARISONS OF COMPARISONS ACROSS SAMPLINGS: 
richness_time_interaction_contrasts <- contrast(em_time_richness, interaction = "pairwise")
summary(richness_time_interaction_contrasts, type = "response")




####### ABUNDANCE ############

glmm_abundance <- glmmTMB(abundance ~ OTC * TILLING * sampling + (1 | plot),
                          data = arkaute_abundance, family = gaussian())
diagnose_glmm(glmm_abundance)
ar_test(glmm_abundance, arkaute_abundance)

em_treat_abundance <- emmeans(glmm_abundance, ~ OTC * TILLING, type = "response")
as.data.frame(pairs(em_treat_abundance, adjust = "tukey")) 

abundance_interaction_contrast <- contrast(em_treat_abundance, interaction = "pairwise")
summary(abundance_interaction_contrast, type = "response")


# Considering temporal dynamics: 
Anova(glmm_abundance, type = "III")

em_time_abundance <- emmeans(glmm_abundance, ~ OTC * TILLING | sampling)
abundance_time_interaction_contrasts <- contrast(em_time_abundance, interaction = "pairwise")
as.data.frame(summary(abundance_time_interaction_contrasts, type = "response"))






####### BIOMASS ############

glmm_biomass <- glmmTMB(biomass_mice_lm ~ OTC * TILLING * sampling + mean_vwc + (1 | plot),
                        dispformula = ~treatment + sampling, 
                        data = arkaute_biomass, family = Gamma(link = "log"))

AIC(glmm_biomass)
diagnose_glmm(glmm_biomass)
ar_test(glmm_biomass, arkaute_biomass)

em_treat_biomass <- emmeans(glmm_biomass, ~ OTC * TILLING, type = "response")
as.data.frame(pairs(em_treat_biomass, adjust = "tukey")) 

biomass_interaction_contrast <- contrast(em_treat_biomass, interaction = "pairwise")
summary(biomass_interaction_contrast, type = "response")


# Considering temporal dynamics: 
Anova(glmm_biomass, type = "III")

em_time_biomass <- emmeans(glmm_biomass, ~ OTC * TILLING | sampling)
biomass_time_interaction_contrasts <- contrast(em_time_biomass, interaction = "pairwise")
as.data.frame(summary(biomass_time_interaction_contrasts, type = "response"))


