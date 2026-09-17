



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

# We have to look at Radj (the higher, the better, to k-index (should be >1)
# and to DHARMA dispersion tests


diagnose_gam <- function(model, data = NULL, group_var = NULL) {
  
  cat("========================================\n")
  cat("           GAM MODEL SUMMARY            \n")
  cat("========================================\n")
  print(summary(model))
  
  # GAM-specific diagnostic (checking k dimension and basis fit)
  cat("\n========================================\n")
  cat("      GAM BASIS DIAGNOSTICS (mgcv)      \n")
  cat("========================================\n")
  
  AIC(model)
  cat("-> Evaluates if the basis dimension (k) is adequate (k-index < 1 or p < 0.05 suggests increasing 'k').\n\n")
  mgcv::gam.check(model)
  
  # Simulate scaled residuals with DHARMa
  sim <- DHARMa::simulateResiduals(fittedModel = model)
  plot(sim)
  
  # --- DHARMa DIAGNOSTIC TESTS ---
  cat("\n========================================\n")
  cat("        DHARMa DIAGNOSTIC TESTS          \n")
  cat("========================================\n")
  
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
    cat("\n----------------------------------------\n")
    cat("5. Levene Test on DHARMa residuals by", group_var, ":\n")
    cat("   -> Evaluates variance homogeneity. Tests if residual spread differs significantly across levels of the grouping factor.\n")
    print(car::leveneTest(sim$scaledResiduals ~ grp))
    DHARMa::plotResiduals(sim, form = grp)
  }
  
  invisible(sim)
}





## Richness ##

m1_richness <- gam(
  richness ~ treatment + s(sampling_num, by = treatment, k = 10) + s(plot, bs = "re"),
  data = arkaute_richness, family = gaussian())

m2_richness <- gam(
  richness ~ treatment + s(sampling_num, by = treatment, k = 15) + s(plot, bs = "re"),
  data = arkaute_richness, family = poisson())

m3_richness <- gam(
  richness ~ treatment +  s(sampling_num, by = treatment, k = 15) +  s(plot, bs = "re"),
  data = arkaute_richness,family = gaussian()
)


diagnose_gam(m1_richness)
diagnose_gam(m2_richness)
diagnose_gam(m3_richness)


AIC(m1_richness)
AIC(m2_richness)
AIC(m3_richness)



## Abundance ##

m1_abundance <- gam(abundance ~ treatment + s(sampling_num, by = treatment, k = 10) + s(plot, bs = "re"),
                     data = arkaute_abundance, family = gaussian(link = "identity"))

m2_abundance <- gam(
  abundance ~ treatment + s(sampling_num, by = treatment, k = 18) + 
    s(plot, bs = "re"), data = arkaute_abundance,  family = gaussian(link = "identity"))

m3_abundance <- gam(
  abundance ~ treatment + s(sampling_num, by = treatment, k = 18) + 
    s(plot, bs = "re"), data = arkaute_abundance,  family = gaussian(link = "log"))



diagnose_gam(m1_abundance)
diagnose_gam(m2_abundance)
diagnose_gam(m3_abundance)

AIC(m1_abundance)
AIC(m2_abundance)
AIC(m3_abundance)


## Evenness ## 
quantile(arkaute_evenness$Y_zipf)
m1_evenness <- gam(Y_zipf ~ treatment + s(sampling_num, by = treatment, k = 10) + s(plot, bs = "re"),
                    data = arkaute_evenness, family = gaussian())

m2_evenness <- gam(Y_zipf ~ treatment + s(sampling_num, by = treatment, k = 15) + s(plot, bs = "re"),
                   data = arkaute_evenness, family = gaussian())

m3_evenness <- gam( Y_zipf ~ treatment + s(sampling_num, by = treatment, k = 15) + s(plot, bs = "re"),
  data = arkaute_evenness, family = scat()
)

diagnose_gam(m1_evenness)
diagnose_gam(m2_evenness)
diagnose_gam(m3_evenness)

AIC(m1_evenness)
AIC(m2_evenness)
AIC(m3_evenness)



## SLA ##

m1_sla <- gam(SLA ~ treatment + s(sampling_num, by = treatment, k = 15) + s(plot, bs = "re"),
              data = arkaute_sla, family = gaussian(link = "log"))


m2_sla <- gam(SLA ~ treatment + s(sampling_num, by = treatment, k = 10) + s(plot, bs = "re"),
               data = arkaute_sla, family = gaussian(link = "log"))


diagnose_gam(m1_sla)
diagnose_gam(m2_sla)

AIC(m1_sla)
AIC(m2_sla)


## LDMC ##

m1_ldmc <- gam(LDMC ~ treatment + s(sampling_num, by = treatment, k = 10) + s(plot, bs = "re"),
                data = arkaute_ldmc, family = gaussian())

m2_ldmc <- gam(LDMC ~ treatment + s(sampling_num, by = treatment, k = 18) + s(plot, bs = "re"),
               data = arkaute_ldmc, family = gaussian())

diagnose_gam(m1_ldmc)
diagnose_gam(m2_ldmc)

AIC(m1_ldmc)
AIC(m2_ldmc)


## LeafN ## 

m1_leafN <- gam(leafN ~ treatment + s(sampling_num, by = treatment, k = 10) + s(plot, bs = "re"),
                data = arkaute_leafN, family = gaussian(link = "identity"))

m2_leafN <- gam(leafN ~ treatment + s(sampling_num, by = treatment, k = 18) + s(plot, bs = "re"),
                data = arkaute_leafN, family = gaussian(link = "identity"))

diagnose_gam(m1_leafN)
diagnose_gam(m2_leafN)

AIC(m1_leafN)
AIC(m2_leafN)



## Biomass ##

m1_biomass <- gam(biomass_mice_lm ~ treatment + s(sampling_num, by = treatment, k = 10) + s(plot, bs = "re"),
                   data = arkaute_biomass, family = tw(link = "log"))

m2_biomass <- gam(biomass_mice_lm ~ treatment + s(sampling_num, by = treatment, k = 15) + s(plot, bs = "re"),
                  data = arkaute_biomass, family = tw(link = "log"))


diagnose_gam(m1_biomass)
diagnose_gam(m2_biomass)

AIC(m1_biomass)
AIC(m2_biomass)

