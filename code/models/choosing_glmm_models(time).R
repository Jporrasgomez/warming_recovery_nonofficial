



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


### Richness ####

m1_richness <- glmmTMB( richness ~ treatment * sampling + ar1(sampling + 0 | plot),
               data = arkaute_richness,family = genpois())


m2_richness <- glmmTMB( richness ~ treatment * poly(sampling_num, 3) + ar1(sampling + 0 | plot),
                          data = arkaute_richness,family = genpois())


m3_richness <- glmmTMB( richness ~ treatment * sampling + ar1(sampling + 0 | plot),
                        dispformula = ~ treatment, 
                        data = arkaute_richness,family = genpois())


as.data.frame(pairs(emmeans(m3_richness, ~ treatment, type = "response"), adjust = "tukey"))
pairs(emmeans(m3_richness, ~ sampling, type = "response"), adjust = "tukey")

diagnose_glmm(m1_richness)
diagnose_glmm(m2_richness)
diagnose_glmm(m3_richness)
diagnose_glmm(m4_richness)
  
AIC(m1_richness)
AIC(m2_richness)
AIC(m3_richness)
AIC(m4_richness)



### Abundance ###


m1_abundance <- glmmTMB(abundance ~ treatment * sampling + ar1(sampling + 0 | plot),
                          data = arkaute_abundance, family = tweedie(link = "log"))


m2_abundance <- glmmTMB(abundance ~ treatment * poly(sampling_num, 3) + ar1(sampling_f + 0 | plot),
                        data = arkaute_abundance, family = tweedie(link= "log"))

m3_abundance <- glmmTMB(
  abundance ~ treatment * poly(sampling_num, 3) + ar1(sampling_f + 0 | plot),
  data = arkaute_abundance, family = Gamma(link = "log")
)



m4_abundance <- glmmTMB(
  abundance ~ treatment * sampling + ar1(sampling + 0 | plot),
  dispformula = ~ treatment, # Allows variance/dispersion of treatments to be estimated independently
  data = arkaute_abundance,
  family = lognormal(link = "log")
)


diagnose_glmm(m1_abundance)
diagnose_glmm(m2_abundance)
diagnose_glmm(m3_abundance)
diagnose_glmm(m4_abundance)

AIC(m1_abundance)
AIC(m2_abundance)
AIC(m3_abundance)
AIC(m4_abundance)




#### Evenness ###

m1_evenness <- glmmTMB(Y_zipf ~ treatment * sampling + ar1(sampling + 0 | plot),
                         data = arkaute_evenness, family = gaussian())


m2_evenness <- glmmTMB(Y_zipf ~ treatment * poly(sampling_num, 3) + ar1(sampling + 0 | plot),
                      data = arkaute_evenness, family = gaussian())


m3_evenness <- glmmTMB(
  Y_zipf ~ treatment * poly(sampling_num, 3) + ar1(sampling_f + 0 | plot),
  data = arkaute_evenness,
  family = t_family(link = "identity")
)


diagnose_glmm(m1_evenness)
diagnose_glmm(m2_evenness)
diagnose_glmm(m3_evenness)

AIC(m1_evenness)
AIC(m2_evenness)
AIC(m3_evenness)





### SLA ###

hist(arkaute_sla$SLA, na.rm = T, breaks = 50)

m1_sla <- glmmTMB(SLA ~ treatment * sampling + ar1(sampling + 0 | plot),
                    data = arkaute_sla, family = gaussian(link = "log"))

m2_sla <- glmmTMB(SLA ~ treatment * sampling + ar1(sampling + 0 | plot),
                  data = arkaute_sla, family = Gamma(link = "log"))

m3_sla <- glmmTMB(
  SLA ~ treatment * sampling + ar1(sampling + 0 | plot),
  dispformula = ~ treatment,
  data = arkaute_sla, family = lognormal(link = "log")
)

m4_sla <- glmmTMB(SLA ~ treatment * poly(sampling_num, 3) + ar1(sampling_f + 0 | plot),
                  dispformula = ~ treatment,
                  data = arkaute_sla, family = Gamma(link = "log"))


m5_sla <- glmmTMB(SLA ~ treatment * poly(sampling_num, 3) + ar1(sampling_f + 0 | plot),
                  dispformula = ~ treatment,
                  data = arkaute_sla, family = lognormal(link = "log"))


diagnose_glmm(m1_sla)
diagnose_glmm(m2_sla)
diagnose_glmm(m3_sla)
diagnose_glmm(m4_sla)
diagnose_glmm(m5_sla)

AIC(m1_sla)
AIC(m2_sla)
AIC(m3_sla)
AIC(m4_sla)
AIC(m5_sla)



### LDMC ###

hist(arkaute_ldmc$LDMC)
quantile(arkaute_ldmc$LDMC)

m1_LDMC <- glmmTMB( LDMC ~ treatment * sampling + ar1(sampling + 0 | plot),
                      data = arkaute_ldmc, family = gaussian())



m2_LDMC <- glmmTMB( LDMC ~ treatment * poly(sampling_num, 3) + ar1(sampling_f + 0 | plot),
                    dispformula = ~ treatment,
                    data = arkaute_ldmc, family = Gamma(link = "log"))

m3_LDMC <- glmmTMB( LDMC ~ treatment * poly(sampling_num, 3) + ar1(sampling_f + 0 | plot),
                    dispformula = ~treatment, 
                    data = arkaute_ldmc, family = Gamma(link = "log"))


m4_LDMC <- glmmTMB( LDMC ~ treatment * poly(sampling_num, 3) + ar1(sampling_f + 0 | plot),
                    dispformula = ~treatment, 
                    data = arkaute_ldmc, family = lognormal(link = "log"))

m5_LDMC <- glmmTMB( LDMC ~ treatment * poly(sampling_num, 3) + ar1(sampling_f + 0 | plot),
                    dispformula = ~treatment + poly(sampling_num, 2), 
                    data = arkaute_ldmc, family = lognormal(link = "log"))


diagnose_glmm(m1_LDMC)
diagnose_glmm(m2_LDMC)
diagnose_glmm(m3_LDMC)
diagnose_glmm(m4_LDMC)
diagnose_glmm(m5_LDMC)

AIC(m1_LDMC)
AIC(m2_LDMC)
AIC(m3_LDMC)
AIC(m4_LDMC)
AIC(m5_LDMC)



### Leaf N ###

hist(arkaute_leafN$leafN)
quantile(arkaute_leafN$leafN)

m1_leafN <- glmmTMB(leafN ~ treatment * sampling + ar1(sampling + 0 | plot),
                      data = arkaute_leafN, family = gaussian(link = "identity"))

m2_leafN <- glmmTMB(leafN ~ treatment * poly(sampling_num, 3) + ar1(sampling_f + 0 | plot),
                    data = arkaute_leafN, family = gaussian(link = "identity"))

m3_leafN <- glmmTMB(leafN ~ treatment * poly(sampling_num, 3) + ar1(sampling_f + 0 | plot),
                    dispformula = ~treatment, 
                    data = arkaute_leafN, family = gaussian(link = "identity"))



diagnose_glmm(m1_leafN)
diagnose_glmm(m2_leafN)
diagnose_glmm(m3_leafN)

AIC(m1_leafN)
AIC(m2_leafN)
AIC(m3_leafN)


### Biomass ###


m1_biomass <- glmmTMB(biomass_mice_lm ~ treatment * sampling + ar1(sampling + 0 | plot),
                        data = arkaute_biomass,  family = Gamma(link = "log"))

m2_biomass <- glmmTMB(biomass_mice_lm ~ treatment * poly(sampling_num, 3)  + ar1(sampling_f + 0 | plot),
                      data = arkaute_biomass,  family = Gamma(link = "log"))

m3_biomass <- glmmTMB(biomass_mice_lm ~ treatment * sampling + ar1(sampling + 0 | plot),
                      dispformula = ~treatment,
                      data = arkaute_biomass,  family = Gamma(link = "log"))

m4_biomass <- glmmTMB(biomass_mice_lm ~ treatment * sampling + ar1(sampling + 0 | plot),
                      dispformula = ~treatment,
                      data = arkaute_biomass,  family = lognormal(link = "log"))



diagnose_glmm(m1_biomass)
diagnose_glmm(m2_biomass)
diagnose_glmm(m3_biomass)
diagnose_glmm(m4_biomass)

diagnose_glmm(m3_biomass, data = arkaute_biomass, group_var = "treatment")


AIC(m1_biomass)
AIC(m2_biomass)
AIC(m3_biomass)
AIC(m4_biomass)



