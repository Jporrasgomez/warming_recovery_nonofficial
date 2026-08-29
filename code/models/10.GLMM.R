




# =================================================================================
# 0. Clear environment
# =================================================================================
rm(list = ls(all.names = TRUE))

# =================================================================================
# 1. Load packages
# =================================================================================
pacman::p_load(
  dplyr, reshape2, tidyverse, lubridate,
  ggplot2, ggpubr, gridExtra,
  car, ggsignif, dunn.test, rstatix,
  lme4, nlme, glmmTMB, performance,
  emmeans, DHARMa
)

# Custom palettes and labels
source("code/palettes_labels.R")

# Global ggplot2 theme
theme_set(
  theme_bw() +
    theme(
      legend.position    = "right",
      panel.grid         = element_blank(),
      strip.background   = element_blank(),
      strip.text         = element_text(face = "bold"),
      text               = element_text(size = 11)
    )
)

# =================================================================================
# 2. Data import & preprocessing
# =================================================================================
arkaute <- read.csv("data/processed_data/arkaute.csv") %>%
  mutate(
    year             = factor(year),
    date             = ymd(date),
    #omw_date         = factor(omw_date),
    #one_month_window = factor(one_month_window),
    sampling         = factor(sampling),
    plot             = factor(plot),
    treatment        = factor(treatment)
  ) %>%
  # drop pre-sampling and noisy first sampling for p & wp
  filter(sampling != "0") %>%
  filter(!(sampling == "1" & treatment %in% c("p", "wp")))


# =================================================================================
# 3. Exploratory check: collinearity of mean_vwc and treatment
# =================================================================================
# Boxplots
ggplot(arkaute, aes(treatment, mean_vwc, fill = treatment)) +
  geom_boxplot(alpha = 0.6) +
  geom_jitter(width = 0.1, alpha = 0.3) +
  labs(title = "Soil moisture by treatment",
       x = "Treatment", y = "mean_vwc")

ggplot(arkaute, aes(plot, mean_vwc)) +
  geom_boxplot(fill = "lightblue") +
  labs(title = "Soil moisture by plot",
       x = "Plot", y = "mean_vwc")

# Linear model R²
colmod <- lm(mean_vwc ~ treatment, data = arkaute)
summary(colmod)  # R² ~0.08

# =================================================================================
# 4. Diagnostic helper functions
# =================================================================================

# 4.1 Overdispersion test for nlme::lme
test_overdispersion_lme <- function(model) {
  resid_p <- residuals(model, type = "pearson")
  rdf     <- length(resid_p) - length(fixef(model))
  chi2    <- sum(resid_p^2)
  ratio   <- chi2 / rdf
  p_value <- pchisq(chi2, df = rdf, lower.tail = FALSE)
  cat("Overdispersion test:\n")
  cat("  Chi2 =", round(chi2, 2),
      " on", rdf, "df → dispersion =", round(ratio, 2), "\n")
  cat("  P =", signif(p_value, 3), "\n\n")
}

# 4.2 GLMM diagnostics via DHARMa (including heteroscedasticity)
diagnose_glmm <- function(model, data = NULL, group_var = NULL) {
  # 1) Summary
  print(summary(model))
  # 2) Simulate
  sim <- DHARMa::simulateResiduals(fittedModel = model)
  # 3) DHARMa standard plots
  plot(sim)
  # 4) Dispersion
  print(DHARMa::testDispersion(sim))
  
  # 5) Heteroscedasticity via Levene on scaled residuals
  if (!is.null(data) && !is.null(group_var)) {
    grp  <- data[[group_var]]
    # extract the scaled (uniform) residuals
    resu <- sim$scaledResiduals
    cat("Levene test on DHARMa residuals by", group_var, ":\n")
    print(car::leveneTest(resu ~ grp))
    # also plot residuals vs group
    plotResiduals(sim, form = grp)
  }
  
  invisible(sim)
}



# =================================================================================
# 5. Modeling workflow
#    For each response: GLMM → diagnose_glmm() → emmeans
#                      LME  → diagnose_lme()  → emmeans
# =================================================================================

# 5.1 Richness
glmm_richness <- glmmTMB(
  richness ~ treatment + scale(mean_vwc)+ (1 | plot),
  data   = arkaute,
  family = nbinom2(link = "log")
)

diagnose_glmm(
  glmm_richness,
  data      = arkaute,
  group_var = "treatment"
)

emmeans(glmm_richness, pairwise ~ treatment)



# 5.2 Abundance
# GLMM Gamma
glmm_abundance <- glmmTMB(
  abundance ~ treatment * scale(mean_vwc) + (1 | plot),
  data   = arkaute,
  family = Gamma(link = "log")
)

diagnose_glmm(
  glmm_abundance,
  data      = arkaute,
  group_var = "treatment"
)

# GLMM Tweedie
glmm_abundance_tw <- glmmTMB(
  abundance ~ treatment + scale(mean_vwc) + (1 | plot),
  data   = arkaute,
  family = tweedie(link = "log")
)
diagnose_glmm(
  glmm_abundance_tw,
  data      = arkaute,
  group_var = "treatment"
)

# GLMM log(abundance+1)
arkaute <- arkaute %>% mutate(log_abundance = log(abundance + 1))
glmm_abundance_log <- glmmTMB(
  log_abundance ~ treatment + scale(mean_vwc) + (1 | plot),
  data   = arkaute,
  family = gaussian()
)

diagnose_glmm(
  glmm_abundance_log,
  data      = arkaute,
  group_var = "treatment"
)

emmeans(glmm_abundance_log, pairwise ~ treatment, type = "response")


# 5.3 Evenness (Y_zipf)
arkaute_Yzipf      <- arkaute %>% filter(!is.na(Y_zipf))


glmm_evenness <- glmmTMB(
  Y_zipf ~ treatment + scale(mean_vwc) + (1 | plot),
  data   = arkaute_Yzipf,
  family = gaussian()
)
diagnose_glmm(
  glmm_evenness,
  data      = arkaute_Yzipf,
  group_var = "treatment"
)
summary(glmm_evenness)
emmeans(glmm_evenness, pairwise ~ treatment)




# 5.4 Biomass
arkaute_biomass   <- arkaute %>% filter(!is.na(biomass_mice_lm))

glmm_biomass <- glmmTMB(
  biomass_mice_lm ~ treatment + scale(mean_vwc) + (1 | plot),
  data   = arkaute_biomass,
  family = Gamma(link = "log")
)

diagnose_glmm(
  glmm_biomass,
  data      = arkaute_biomass,
  group_var = "treatment"
)
summary(glmm_biomass)
emmeans(glmm_biomass, pairwise ~ treatment)
emmeans(glmm_biomass, pairwise ~ treatment, type = "response")





# SLA

arkaute_SLA   <- arkaute %>% filter(!is.na(SLA))

glmm_SLA <- glmmTMB(
  SLA ~ treatment + scale(mean_vwc) + (1 | plot),
  data   = arkaute_biomass,
  family = Gamma(link = "log")
)

diagnose_glmm(
  glmm_SLA,
  data      = arkaute_biomass,
  group_var = "treatment"
)
summary(glmm_SLA)
emmeans(glmm_SLA, pairwise ~ treatment, type = "response")




# LDMC

arkaute_LDMC   <- arkaute %>% filter(!is.na(LDMC))

glmm_LDMC <- glmmTMB(
  LDMC ~ treatment + scale(mean_vwc) + (1 | plot),
  data   = arkaute_biomass,
  family = Gamma(link = "log")
)

diagnose_glmm(
  glmm_LDMC,
  data      = arkaute_biomass,
  group_var = "treatment"
)
summary(glmm_LDMC)
emmeans(glmm_LDMC, pairwise ~ treatment, type = "response")



# leafN

arkaute_leafN   <- arkaute %>% filter(!is.na(leafN))

glmm_leafN <- glmmTMB(
  leafN ~ treatment + scale(mean_vwc) + (1 | plot),
  data   = arkaute_biomass,
  family = Gamma(link = "log")
)

diagnose_glmm(
  glmm_leafN,
  data      = arkaute_biomass,
  group_var = "treatment"
)
summary(glmm_leafN)
emmeans(glmm_leafN, pairwise ~ treatment, type = "response")



