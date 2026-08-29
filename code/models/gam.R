

rm(list = ls(all.names = TRUE))

# Cargar paquetes (se añade mgcv)
pacman::p_load(
  dplyr, reshape2, tidyverse, lubridate, ggplot2, ggpubr, gridExtra,
  car, ggsignif, dunn.test, rstatix, performance, emmeans, DHARMa, mgcv
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

arkaute <- read.csv("data/processed_data/arkaute.csv") %>%
  mutate(
    year      = factor(year),
    date      = ymd(date),
    sampling  = as.numeric(as.character(sampling)), # Para suavizados s(sampling) conviene que sea numérico
    plot      = factor(plot),
    treatment = factor(treatment)
  ) %>%
  # drop pre-sampling and noisy first sampling for p & wp
  filter(sampling != 0) %>%
  filter(!(sampling == 1 & treatment %in% c("p", "wp")))

# Visualizaciones iniciales
ggplot(arkaute, aes(treatment, mean_vwc, fill = treatment)) +
  geom_boxplot(alpha = 0.6) +
  geom_jitter(width = 0.1, alpha = 0.3) +
  labs(title = "Soil moisture by treatment", x = "Treatment", y = "mean_vwc")

ggplot(arkaute, aes(plot, mean_vwc)) +
  geom_boxplot(fill = "lightblue") +
  labs(title = "Soil moisture by plot", x = "Plot", y = "mean_vwc")

# Linear model R²
colmod <- lm(mean_vwc ~ treatment, data = arkaute)
summary(colmod)  # R² ~0.08

#############
# Diagnósticos adaptados a GAM (mgcv)
#############

diagnose_gam <- function(model, data = NULL, group_var = NULL) {
  # 1) Resumen del GAM
  print(summary(model))
  
  # 2) Simulación de residuos con DHARMa
  sim <- DHARMa::simulateResiduals(fittedModel = model)
  
  # 3) Diagnósticos gráficos DHARMa
  plot(sim)
  
  # 4) Test de dispersión
  print(DHARMa::testDispersion(sim))
  
  # 5) Diagnósticos nativos de mgcv
  gam.check(model)
  
  # 6) Heterocedasticidad vía Levene
  if (!is.null(data) && !is.null(group_var)) {
    grp  <- data[[group_var]]
    resu <- sim$scaledResiduals
    cat("Levene test on DHARMa residuals by", group_var, ":\n")
    print(car::leveneTest(resu ~ grp))
    plotResiduals(sim, form = grp)
  }
  
  invisible(sim)
}


# =================================================================================
# Workflow con GAMs (mgcv)
# Nota: (1|plot) pasa a ser s(plot, bs = "re")
# =================================================================================

# 5.1 Richness

# GAM temporal (curva suave por tratamiento + efecto aleatorio de parcela)
gam_richness_time <- gam(
  richness ~ treatment + s(sampling, by = treatment, k = 4) + s(plot, bs = "re"),
  data   = arkaute,
  family = poisson(link = "log") # Alternativa no-entera/sobredispersa: nb() o tw()
)

diagnose_gam(gam_richness_time, data = arkaute, group_var = "treatment")
emmeans(gam_richness_time, pairwise ~ treatment)


# GAM sin componente temporal
gam_richness <- gam(
  richness ~ treatment + s(plot, bs = "re"),
  data   = arkaute,
  family = nb(link = "log")
)

diagnose_gam(gam_richness, data = arkaute, group_var = "treatment")
emmeans(gam_richness, pairwise ~ treatment)


# 5.2 Abundance

hist(arkaute$abundance, breaks = 50)

abundance_q1  <- quantile(arkaute$abundance, 0.25, na.rm = TRUE)
abundance_q3  <- quantile(arkaute$abundance, 0.75, na.rm = TRUE)
abundance_IQR <- abundance_q3 - abundance_q1

abundance_noOutliers <- arkaute %>%
  select(sampling, treatment, plot, abundance) %>%
  filter(abundance < (abundance_q3 + 1.5 * abundance_IQR))

hist(abundance_noOutliers$abundance, breaks = 50)

# Tweedie (tw) en mgcv maneja perfectamente ceros y asimetría en abundancias
gam_abundance <- gam(
  abundance ~ treatment 
  #+ s(sampling, by = treatment, k = 4) 
  + s(plot, bs = "re"),
  data   = abundance_noOutliers,
  family = Tweedie(p = 1.5, link = "log")
)


diagnose_gam(gam_abundance, data = abundance_noOutliers, group_var = "treatment")
emmeans(gam_abundance, pairwise ~ treatment, type = "response")


# 5.3 Evenness (Y_zipf)

arkaute_Yzipf <- arkaute %>% 
  filter(!is.na(Y_zipf)) %>% 
  mutate(Y_zipf = Y_zipf * -1)

hist(arkaute_Yzipf$Y_zipf, breaks = 50)

gam_evenness <- gam(
  Y_zipf ~ treatment + s(plot, bs = "re"),
  data   = arkaute_Yzipf,
  family = gaussian(link = "log")
)

diagnose_gam(gam_evenness, data = arkaute_Yzipf, group_var = "treatment")
emmeans(gam_evenness, pairwise ~ treatment, type = "response")


# 5.4 Biomass

arkaute_biomass <- arkaute %>% filter(!is.na(biomass_mice_lm))
hist(arkaute_biomass$biomass_mice_lm, breaks = 50)

gam_biomass <- gam(
  biomass_mice_lm ~ treatment + s(plot, bs = "re"),
  data   = arkaute_biomass,
  family = Gamma(link = "log")
)

diagnose_gam(gam_biomass, data = arkaute_biomass, group_var = "treatment")
emmeans(gam_biomass, pairwise ~ treatment, type = "response")


# 5.5 SLA

arkaute_SLA <- arkaute %>% filter(!is.na(SLA))
hist(arkaute_SLA$SLA, breaks = 50)

gam_SLA <- gam(
  SLA ~ treatment + s(plot, bs = "re"),
  data   = arkaute_SLA,
  family = Gamma(link = "log")
)

# Corrección: pasamos la tabla correcta 'arkaute_SLA' a la función de diagnóstico
diagnose_gam(gam_SLA, data = arkaute_SLA, group_var = "treatment")
emmeans(gam_SLA, pairwise ~ treatment, type = "response")


# 5.6 LDMC

arkaute_LDMC <- arkaute %>% filter(!is.na(LDMC))
hist(arkaute_LDMC$LDMC, breaks = 50)

gam_LDMC <- gam(
  LDMC ~ treatment + s(plot, bs = "re"),
  data   = arkaute_LDMC,
  family = gaussian(link = "identity")
)

diagnose_gam(gam_LDMC, data = arkaute_LDMC, group_var = "treatment")
emmeans(gam_LDMC, pairwise ~ treatment, type = "response")


# 5.7 leafN

arkaute_leafN <- arkaute %>% filter(!is.na(leafN))
hist(arkaute_leafN$leafN, breaks = 50)

gam_leafN <- gam(
  leafN ~ treatment + s(plot, bs = "re"),
  data   = arkaute_leafN,
  family = gaussian(link = "identity")
)

diagnose_gam(gam_leafN, data = arkaute_leafN, group_var = "treatment")
emmeans(gam_leafN, pairwise ~ treatment, type = "response")
