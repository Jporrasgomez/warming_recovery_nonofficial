


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
    sampling = factor(sampling, levels = as.character(sort(unique(as.numeric(as.character(sampling)))))),
    plot      = factor(plot),
    treatment = factor(treatment)
  ) %>%
  filter(sampling != "0") |>
  filter(sampling != "1") 
#filter(!(sampling == "1" & treatment %in% c("p", "wp")))


  # databases for models:
  arkaute_richness <- arkaute |>  filter(!is.na(richness))
  arkaute_abundance <- arkaute |>  filter(!is.na(abundance))
  arkaute_evenness <- arkaute |> filter(!is.na(Y_zipf))
  arkaute_sla <- arkaute |> filter(!is.na(SLA))
  arkaute_ldmc <- arkaute |> filter(!is.na(LDMC))
  arkaute_leafN <- arkaute |> filter(!is.na(leafN))
  arkaute_biomass <- arkaute |> filter(!is.na(biomass_mice_lm))
  
  # Diagnosis functions
  
  

  
  
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
  
  
  ########## Choosing GLMMs ###########
  
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
  
  
  
  ##### RICHNESS #####
  
  m1_glmm_richness <- glmmTMB(richness ~ treatment + (1 | plot), data = arkaute, family = gaussian)
  
  diagnose_glmm(m1_glmm_richness)
 
  
  ##### ABUNDANCE #####
  m1_glmm_abundance <- glmmTMB(abundance ~ treatment + (1 | plot),
                               data = arkaute, family = gaussian(link = "identity"))
  
  m2_glmm_abundance <- glmmTMB(abundance ~ treatment + (1 | plot),
                               dispformula =  ~ treatment, 
                               data = arkaute, family = gaussian(link = "identity"))
  
  m3_glmm_abundance <- glmmTMB(abundance ~ treatment + (1 | plot), 
                       data = arkaute, 
                       family = tweedie())
  
  # 2. Modelo Binomial Negativa con varianza por tratamiento
  m4_glmm_abundance <- glmmTMB(abundance ~ treatment + (1 | plot), 
                     data = arkaute, 
                     family = Gamma(link = "log"))
  
  diagnose_glmm(m1_glmm_abundance)
  diagnose_glmm(m2_glmm_abundance)
  diagnose_glmm(m3_glmm_abundance)
  diagnose_glmm(m4_glmm_abundance)
  
  AIC(m1_glmm_abundance)
  AIC(m2_glmm_abundance)
  AIC(m3_glmm_abundance)

  
  ##### EVENNESS #####
  
  m1_glmm_evenness <- glmmTMB(Y_zipf ~ treatment + (1 | plot), data = arkaute_evenness, family = gaussian())
  m2_glmm_evenness <- glmmTMB(Y_zipf ~ treatment + (1 | plot),
                              dispformula = ~ treatment, 
                              data = arkaute_evenness, family = gaussian())
  diagnose_glmm(m1_glmm_evenness)
  diagnose_glmm(m2_glmm_evenness)

  
  ### SLA ###
  m1_glmm_sla <- glmmTMB(SLA ~ treatment + (1 | plot), data = arkaute_sla, family = gaussian(link = "log"))
  
  m2_glmm_sla <- glmmTMB(SLA ~ treatment + (1 | plot), data = arkaute_sla, family = Gamma(link = "log"))
  
  m3_glmm_sla <- glmmTMB(SLA ~ treatment + (1 | plot), data = arkaute_sla, family = lognormal())
  
  m4_glmm_sla <- glmmTMB(SLA ~ treatment + (1 | plot), 
                         dispformula = ~ treatment, 
                         data = arkaute_sla, family = gaussian(link = "log"))
  

  
  AIC(m1_glmm_sla, m2_glmm_sla, m3_glmm_sla, m4_glmm_sla, m5_glmm_sla)
  diagnose_glmm(null_glmm_sla)
  diagnose_glmm(m4_glmm_sla)
  
  
  ### LDMC ###
  m1_glmm_LDMC <- glmmTMB(LDMC ~ treatment + (1 | plot), data = arkaute_ldmc, family = gaussian())


  m2_glmm_LDMC <- glmmTMB(LDMC ~ treatment + (1 | plot), 
                          dispformula = ~ treatment, 
                          data = arkaute_ldmc, 
                          family = gaussian())

  AIC(m1_glmm_LDMC, m2_glmm_LDMC)
  
  diagnose_glmm(m1_glmm_LDMC)
  diagnose_glmm(m2_glmm_LDMC)
  
  
  ### Leaf nitrogen ###
  m1_glmm_leafN <- glmmTMB(leafN ~ treatment + (1 | plot), data = arkaute_leafN, family = gaussian(link = "identity"))
  m2_glmm_leafN <- glmmTMB(leafN ~ treatment + (1 | plot),
                           dispformula = ~ treatment, 
                           data = arkaute_leafN, family = gaussian(link = "identity"))
  AIC(m1_glmm_leafN, m2_glmm_leafN)
  diagnose_glmm(m1_glmm_leafN)
  diagnose_glmm(m2_glmm_leafN)

  
  ### BIOMASS ###
  #m1_glmm_biomass <- glmmTMB(biomass_mice_lm ~ treatment + (1 | plot), data = arkaute_biomass, family = Gamma(link = "log"))
  m1_glmm_biomass <- glmmTMB(biomass_mice_lm ~ treatment + (1 | plot), data = arkaute_biomass, family = Gamma(link = "log"))
  m2_glmm_biomass <- glmmTMB(biomass_mice_lm ~ treatment + (1 | plot),
                             dispformula = ~ treatment, 
                             data = arkaute_biomass, family = Gamma(link = "log"))
  diagnose_glmm(m1_glmm_biomass)
  diagnose_glmm(m2_glmm_biomass)
  AIC(m1_glmm_biomass, m2_glmm_biomass)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ########## Choosing GAMs ###########
  
  
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
  
  
  
  
  ############# GAM MODELS ##############
  ## Richness ## 
  m1_gam_richness <- gam( richness ~ treatment + s(plot, bs = "re"), data = arkaute_richness, family = gaussian())
  diagnose_gam(m1_gam_richness)
 
  
  ## Abundance ##
  m1_gam_abundance <- gam(abundance ~ treatment + s(plot, bs = "re"), data = arkaute_abundance,
                          family = gaussian(link = "identity"))
  
  m2_gam_abundance <- gam(abundance ~ treatment + s(plot, bs = "re"), data = arkaute_abundance,
                          family = Gamma(link = "log"))
  
  AIC(m1_gam_abundance, m2_gam_abundance)
  diagnose_gam(m1_gam_abundance)
  diagnose_gam(m2_gam_abundance)
  
  
  ## Evenness ##
  m1_gam_evenness <- gam(Y_zipf ~ treatment + s(plot, bs = "re"), data = arkaute_evenness, family = gaussian())
  
  m1_gam_reml <- gam(Y_zipf ~ treatment + s(plot, bs = "re"), 
                     data = arkaute_evenness, 
                     family = gaussian(), 
                     method = "REML")
  AIC(m1_gam_evenness, m1_gam_reml)
  diagnose_gam(m1_gam_evenness)
  diagnose_gam(m1_gam_reml)
  
  
  ## SLA ##
  m1_gam_sla <- gam(SLA ~ treatment + s(plot, bs = "re"), data = arkaute_sla, family = gaussian(link = "log"))
  diagnose_gam(m1_gam_sla)
  
  
  ## LDMC ##
  m1_gam_ldmc <- gam(LDMC ~ treatment + s(plot, bs = "re"), data = arkaute_ldmc, family = gaussian())
  diagnose_gam(m1_gam_ldmc)
  
  
  ## Leaf Nitrogen ##
  m1_gam_leafN <- gam(leafN ~ treatment + s(plot, bs = "re"), data = arkaute_leafN, family = gaussian(link = "identity"))
  diagnose_gam(m1_gam_leafN)
  
  
  ## Biomass ##
  m1_gam_biomass <- gam(biomass_mice_lm ~ treatment + s(plot, bs = "re"), data = arkaute_biomass, family = Gamma(link = "log"))
  diagnose_gam(m1_gam_biomass)
  
  
  