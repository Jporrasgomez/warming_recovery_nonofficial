


rm(list = ls(all.names = TRUE))
pacman::p_load(dplyr, reshape2, tidyverse, lubridate, ggplot2, ggpubr, gridExtra,
  car, ggsignif, dunn.test, rstatix, lme4, nlme, glmmTMB, performance,
  emmeans, DHARMa, glmm
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


arkaute <- read.csv("data/processed_data/arkaute.csv") %>%
  mutate(
    year             = factor(year),
    date             = ymd(date),
    sampling         = factor(sampling),
    plot             = factor(plot),
    treatment        = factor(treatment)
  ) %>%
  # drop pre-sampling and noisy first sampling for p & wp
  filter(sampling != "0") %>%
  filter(!(sampling == "1" & treatment %in% c("p", "wp")))



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

#############
# 4. Diagnostic helper functions
#############

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

glmm_richness_time <- glmmTMB(
  richness ~ treatment * sampling + (1 | plot),
  data   = arkaute,
  family = genpois(link = "log")
)

diagnose_glmm(
  glmm_richness_time,
  data      = arkaute,
  group_var = "treatment"
)

emmeans(glmm_richness_time, pairwise ~ treatment)
emmeans(glmm_richness_time, pairwise ~ treatment | sampling, type = "response")




glmm_richness <- glmmTMB(
  richness ~ treatment + (1 | plot),
  data   = arkaute,
  family = genpois(link = "log")
)

diagnose_glmm(
  glmm_richness,
  data      = arkaute,
  group_var = "treatment"
)

emmeans(glmm_richness, pairwise ~ treatment)



# 5.2 Abundance

# No consigo ajustar el modelo a la variable abundance

# Consiste en hallar el primer cuartil (Q₁) y el tercer cuartil (Q₃),
# calcular el IQR = Q₃ - Q₁, y definir límites: cualquier valor menor a
# Q₁ - 1.5 × IQR o mayor a Q₃ + 1.5 × IQR se considera un outlier.

hist(arkaute$abundance, breaks = 50)

abundance_q1 <- quantile(arkaute$abundance)[1]
abundance_q3 <- quantile(arkaute$abundance)[3]
abundance_IQR <- abundance_q3 - abundance_q1


abundance_noOutliers <- arkaute |> 
  select(sampling, treatment, plot, abundance) |> 
  filter(
    abundance < (abundance_q3 + 1.5 * abundance_IQR)
  )

hist(abundance_noOutliers$abundance, breaks = 50)
  
glmm_abundance <- glmmTMB(
  abundance ~ treatment  + (1 | plot),
  data   = arkaute,
  family = tweedie(link = "log")
)

diagnose_glmm(
  glmm_abundance,
  data      = arkaute,
  group_var = "treatment"
)


emmeans(glmm_abundance, pairwise ~ treatment, type = "response")


# 5.3 Evenness (Y_zipf)
arkaute_Yzipf <- arkaute %>% filter(!is.na(Y_zipf)) |>  mutate(Y_zipf = Y_zipf*-1)

hist(arkaute_Yzipf$Y_zipf, breaks = 50)

glmm_evenness <- glmmTMB(
  Y_zipf ~ treatment  + (1 | plot),
  data   = arkaute_Yzipf,
  family = gaussian(link = "log")
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
hist(arkaute_biomass$biomass_mice_lm, breaks = 50)

glmm_biomass <- glmmTMB(
  biomass_mice_lm ~ treatment + (1 | plot),
  data   = arkaute_biomass,
  family = Gamma(link = "log")
)

diagnose_glmm(
  glmm_biomass,
  data      = arkaute_biomass,
  group_var = "treatment"
)
summary(glmm_biomass)
emmeans(glmm_biomass, pairwise ~ treatment, type = "response")
emmeans(glmm_richness, pairwise ~ treatment | sampling, type = "response")




# SLA

arkaute_SLA   <- arkaute %>% filter(!is.na(SLA))

hist(arkaute_SLA$SLA, breaks = 50)

glmm_SLA <- glmmTMB(
  SLA ~ treatment + (1 | plot),
  data   = arkaute_SLA,
  family = gaussian(link = "identity")
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
hist(arkaute_LDMC$LDMC, breaks = 50)

glmm_LDMC <- glmmTMB(
  LDMC ~ treatment + (1 | plot),
  data   = arkaute_LDMC,
  family = gaussian(link = "identity")
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
hist(arkaute_leafN$leafN, breaks = 50)



glmm_leafN <- glmmTMB(
  leafN ~ treatment + (1 | plot),
  data   = arkaute_leafN,
  family = gaussian(link = "identity")
)

diagnose_glmm(
  glmm_leafN,
  data      = arkaute_biomass,
  group_var = "treatment"
)
summary(glmm_leafN)
emmeans(glmm_leafN, pairwise ~ treatment, type = "response")



