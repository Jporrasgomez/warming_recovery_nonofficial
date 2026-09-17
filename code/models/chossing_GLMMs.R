





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
    sampling_num = as.numeric(as.character(sampling)), # Necesario para splines en GAM
    sampling = factor(sampling, levels = as.character(sort(unique(as.numeric(as.character(sampling)))))),
    sampling_f = numFactor(sampling_num),
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








####################### CHOOSING GLMM's ########################

### Richness ####

m0_richness <- glmmTMB( richness ~ treatment * sampling + (1 | plot),
                        data = arkaute_richness,family = genpois())
AIC(m0_richness)
diagnose_glmm(m0_richness)
ar_test(m0_richness, arkaute_richness)# If p-value is higher than 0.05 we do not need autocorrelation term

as.data.frame(pairs(emmeans(m0_richness, ~ treatment, type = "response"), adjust = "tukey"))




### Abundance ###

hist(arkaute_abundance$abundance)

m0_abundance <- glmmTMB(abundance ~ treatment * sampling + (1 | plot),
                        data = arkaute_abundance, family = gaussian())

AIC(m0_abundance)
diagnose_glmm(m0_abundance)
ar_test(m0_abundance, arkaute_abundance)
as.data.frame(pairs(emmeans(m0_abundance, ~ treatment, type = "response"), adjust = "tukey"))



#### Evenness ###

hist(arkaute_evenness$Y_zipf)
quantile(arkaute_evenness$Y_zipf)

m0_evenness <- glmmTMB(Y_zipf ~ treatment * sampling + (1 | plot),
                       data = arkaute_evenness, family = gaussian())

AIC(m0_evenness)
diagnose_glmm(m0_evenness)

m1_evenness <- glmmTMB(Y_zipf ~ treatment * sampling + (1 | plot),
                       dispformula = ~treatment + sampling, 
                       data = arkaute_evenness, family = gaussian())
AIC(m1_evenness)
diagnose_glmm(m1_evenness) # Best
ar_test(m1_evenness, arkaute_evenness)
as.data.frame(pairs(emmeans(m1_evenness, ~ treatment, type = "response"), adjust = "tukey"))



### SLA ###

hist(arkaute_sla$SLA, breaks = 50)

m0_sla <- glmmTMB(SLA ~ treatment * sampling + (1 | plot), 
                  data = arkaute_sla, family = gaussian())
AIC(m0_sla)
diagnose_glmm(m0_sla)
ar_test(m0_sla, arkaute_sla)

m1_sla <- glmmTMB(SLA ~ treatment * sampling + (1 | plot), 
                  dispformula = ~ treatment, 
                  data = arkaute_sla, family = gaussian())
AIC(m1_sla) # BEST
diagnose_glmm(m1_sla)
ar_test(m0_sla, arkaute_sla)
as.data.frame(pairs(emmeans(m0_sla, ~ treatment, type = "response"), adjust = "tukey"))



### LDMC ###

hist(arkaute_ldmc$LDMC)
quantile(arkaute_ldmc$LDMC)

m0_ldmc <- glmmTMB(LDMC ~ treatment * sampling + (1 | plot), 
                  data = arkaute_ldmc, family = gaussian())

AIC(m0_ldmc)
diagnose_glmm(m0_ldmc)
ar_test(m0_ldmc, arkaute_ldmc)

m1_ldmc <- glmmTMB(LDMC ~ treatment * sampling + (1 | plot), 
                   dispformula = ~ treatment, 
                   data = arkaute_ldmc, family = gaussian())

AIC(m1_ldmc) # BEST
diagnose_glmm(m1_ldmc)
ar_test(m1_ldmc, arkaute_ldmc)
as.data.frame(pairs(emmeans(m1_ldmc, ~ treatment, type = "response"), adjust = "tukey"))



### Leaf N ###

hist(arkaute_leafN$leafN)
quantile(arkaute_leafN$leafN)


m0_leafN <- glmmTMB(leafN ~ treatment * sampling + (1 | plot),
                   data = arkaute_leafN, family = gaussian())
AIC(m0_leafN)
diagnose_glmm(m0_leafN)
ar_test(m0_leafN, arkaute_leafN)


m1_leafN <- glmmTMB(leafN ~ treatment * sampling + (1 | plot),
                    dispformula = ~ treatment,
                    data = arkaute_leafN, family = gaussian())
AIC(m1_leafN)
diagnose_glmm(m1_leafN)
ar_test(m1_leafN, arkaute_leafN)
as.data.frame(pairs(emmeans(m1_leafN, ~ treatment, type = "response"), adjust = "tukey"))




### Biomass ###
hist(arkaute_biomass$biomass_mice_lm)

m0_biomass <- glmmTMB(biomass_mice_lm ~ treatment * sampling + (1 | plot),
                      data = arkaute_biomass, family = Gamma(link = "log"))
AIC(m0_biomass)
diagnose_glmm(m0_biomass)
ar_test(m0_biomass, arkaute_biomass)

m1_biomass <- glmmTMB(biomass_mice_lm ~ treatment * sampling + (1 | plot),
                      dispformula = ~treatment + sampling, 
                      data = arkaute_biomass, family = Gamma(link = "log"))
AIC(m1_biomass)
diagnose_glmm(m1_biomass)
ar_test(m1_biomass, arkaute_biomass)

as.data.frame(pairs(emmeans(m1_biomass, ~ treatment, type = "response"), adjust = "tukey"))




# Biomass with mice imputation only #

hist(arkaute_biomass_mice$biomass_mice)

m0_biomass_mice <- glmmTMB(biomass_mice ~ treatment * sampling + (1 | plot),
                           data = arkaute_biomass_mice,  family = Gamma(link = "log"))
AIC(m0_biomass_mice)
diagnose_glmm(m0_biomass_mice)
ar_test(m0_biomass_mice, arkaute_biomass_mice)

m1_biomass_mice <- glmmTMB(biomass_mice ~ treatment * sampling + (1 | plot),
                           dispformula = ~treatment + sampling, 
                           data = arkaute_biomass_mice,  family = Gamma(link = "log"))
AIC(m1_biomass_mice) #BEST
diagnose_glmm(m1_biomass_mice)
ar_test(m1_biomass_mice, arkaute_biomass_mice)

as.data.frame(pairs(emmeans(m1_biomass_mice, ~ treatment, type = "response"), adjust = "tukey"))


# Biomass raw #
hist(arkaute_biomass_raw$biomass_raw)
min(arkaute_biomass_raw$biomass_raw)

m0_biomass_raw <- glmmTMB(biomass_raw ~ treatment * sampling + (1 | plot),
                           data = arkaute_biomass_raw,  family = Gamma(link = "log"))
AIC(m0_biomass_raw)
diagnose_glmm(m0_biomass_raw)
ar_test(m0_biomass_raw, arkaute_biomass_raw)

m1_biomass_raw <- glmmTMB(biomass_raw ~ treatment * sampling + (1 | plot),
                           dispformula = ~treatment + sampling, 
                           data = arkaute_biomass_raw,  family = Gamma(link = "log"))
AIC(m1_biomass_raw) #BEST
diagnose_glmm(m1_biomass_raw)
ar_test(m1_biomass_raw, arkaute_biomass_raw)

as.data.frame(pairs(emmeans(m1_biomass_raw, ~ treatment, type = "response"), adjust = "tukey"))


