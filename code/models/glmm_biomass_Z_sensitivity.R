




rm(list = ls(all.names = TRUE))
pacman::p_load(
  dplyr, reshape2, tidyverse, lubridate, ggplot2, ggpubr, gridExtra, ggrepel,
  car, ggsignif, dunn.test, rstatix, performance, emmeans, DHARMa, mgcv, glmmTMB
)

source("code/palettes_labels.R")

biomass_data_all <-   read.csv("data/processed_data/biomass_data_Z.csv") |> select(-X) |> 
  mutate(
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



raw_z_1_2 <- biomass_data_all |>  filter(biomass_level == "biomass_raw") |>
  select(-z_5_6, -z_2_3_original) |>  rename(biomass = z_1_2) |> filter(!is.na(biomass), biomass != 0) 

mice_z_1_2 <- biomass_data_all |>  filter(biomass_level == "biomass_mice") |>
  select(-z_5_6, -z_2_3_original)|>  rename(biomass = z_1_2) |> filter(!is.na(biomass), biomass != 0)

final_z_1_2 <- biomass_data_all |>  filter(biomass_level == "biomass_mice_lm") |>
  select(-z_5_6, -z_2_3_original)|>  rename(biomass = z_1_2) |> filter(!is.na(biomass), biomass != 0)



raw_z_5_6 <- biomass_data_all |>  filter(biomass_level == "biomass_raw") |>
  select(-z_1_2, -z_2_3_original)|>  rename(biomass = z_5_6) |> filter(!is.na(biomass), biomass != 0)

mice_z_5_6 <- biomass_data_all |>  filter(biomass_level == "biomass_mice") |>
  select(-z_1_2, -z_2_3_original)|>  rename(biomass = z_5_6) |> filter(!is.na(biomass), biomass != 0)

final_z_5_6 <- biomass_data_all |>  filter(biomass_level == "biomass_mice_lm") |>
  select(-z_1_2, -z_2_3_original)|>  rename(biomass = z_5_6) |> filter(!is.na(biomass), biomass != 0)



raw_z_2_3 <- biomass_data_all |>  filter(biomass_level == "biomass_raw") |>
  select(-z_1_2, -z_5_6)|>  rename(biomass = z_2_3_original) |> filter(!is.na(biomass), biomass != 0)

mice_z_2_3 <- biomass_data_all |>  filter(biomass_level == "biomass_mice") |>
  select(-z_1_2, -z_5_6)|>  rename(biomass = z_2_3_original) |> filter(!is.na(biomass), biomass != 0)

final_z_2_3 <- biomass_data_all |>  filter(biomass_level == "biomass_mice_lm") |>
  select(-z_1_2, -z_5_6)|>  rename(biomass = z_2_3_original) |> filter(!is.na(biomass), biomass != 0)





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
  print(suppressWarnings(DHARMa::testQuantiles(sim)))
  
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


ar_test <- function(model, data){
  testTemporalAutocorrelation(recalculateResiduals(simulateResiduals(model),
                                                   group = data$sampling_num),
                              time = sort(unique(data$sampling_num)))
  
}# If p-value is higher than 0.05 we do not need autocorrelation term


# I will use the glmm biomass model we use for the main analysis
glmm_biomass_str <- function(dataset) {
  model <- glmmTMB(biomass ~ treatment * sampling + (1 | plot),
                   dispformula = ~ treatment + sampling, 
                   data = dataset,  family = Gamma(link = "log"))
  diagnose_glmm(model)
  ar_test(model, dataset)
}


glmm_biomass_str(raw_z_2_3) # Good fit, no need of ar1 term
glmm_biomass_str(mice_z_2_3) # Good fit, no need of ar1 term
glmm_biomass_str(final_z_2_3) # Good fit, no need of ar1 term

glmm_biomass_str(raw_z_5_6) # Good fit, no need of ar1 term
glmm_biomass_str(mice_z_5_6) # Good fit, no need of ar1 term
glmm_biomass_str(final_z_5_6) # Good fit, no need of ar1 term

glmm_biomass_str(raw_z_1_2) # Good fit, no need of ar1 term
glmm_biomass_str(mice_z_1_2) # Good fit, no need of ar1 term
glmm_biomass_str(final_z_1_2) # Could improve fit, but will accept it. No need of ar1 term









