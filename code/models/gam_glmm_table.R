

rm(list = ls(all.names = TRUE))

# Cargar paquetes
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

# Carga y limpieza de datos
arkaute <- read.csv("data/processed_data/arkaute.csv") %>%
  mutate(
    year      = factor(year),
    date      = ymd(date),
    sampling  = factor(sampling),
    plot      = factor(plot),
    treatment = factor(treatment)
  ) %>%
  filter(sampling != "0") %>%
  filter(!(sampling == "1" & treatment %in% c("p", "wp")))


#############
# Diagnóstico Helper para GAM
#############

diagnose_gam <- function(model, data = NULL, group_var = NULL) {
  print(summary(model))
  sim <- DHARMa::simulateResiduals(fittedModel = model)
  plot(sim)
  print(DHARMa::testDispersion(sim))
  gam.check(model)
  
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
# Lista para almacenar los resultados de emmeans de cada variable
# =================================================================================
results_gam_list <- list()

# Función para extraer los contrastes de emmeans en formato largo
extract_emmeans_contrasts <- function(model, var_name) {
  # Extraer emmeans y contrastes pareados
  emm <- emmeans(model, pairwise ~ treatment)
  contrasts_df <- as.data.frame(emm$contrasts)
  
  # Identificar si la columna de test ratio es t.ratio o z.ratio según la familia
  ratio_col <- if ("t.ratio" %in% colnames(contrasts_df)) "t.ratio" else "z.ratio"
  
  # Contrastes objetivo
  target_contrasts <- c("c - p", "c - w", "c - wp", "p - wp")
  
  # Extraer metadatos del modelo
  fam_name <- model$family$family
  link_name <- model$family$link
  
  # Filtrar y formatear la tabla larga
  res_df <- contrasts_df %>%
    filter(contrast %in% target_contrasts) %>%
    mutate(
      variable = var_name,
      model    = "GAM",
      family   = fam_name,
      link     = link_name,
      contrast = gsub(" ", "", contrast), # Limpia los espacios: "c - p" -> "c-p"
      ratio    = .data[[ratio_col]],
      p_value  = p.value
    ) %>%
    select(variable, model, family, link, contrast, ratio, p_value)
  
  return(res_df)
}


# =================================================================================
# Modelado GAM (sin componente temporal)
# =================================================================================

# 1. Richness
gam_richness <- gam(
  richness ~ treatment + s(plot, bs = "re"),
  data   = arkaute,
  family = poisson(link = "log")
)
diagnose_gam(gam_richness, data = arkaute, group_var = "treatment")
results_gam_list[["richness"]] <- extract_emmeans_contrasts(gam_richness, "richness")


# 2. Abundance
abundance_q1  <- quantile(arkaute$abundance, 0.25, na.rm = TRUE)
abundance_q3  <- quantile(arkaute$abundance, 0.75, na.rm = TRUE)
abundance_IQR <- abundance_q3 - abundance_q1

abundance_noOutliers <- arkaute %>%
  select(sampling, treatment, plot, abundance) %>%
  filter(abundance < (abundance_q3 + 1.5 * abundance_IQR))

gam_abundance <- gam(
  abundance ~ treatment + s(plot, bs = "re"),
  data   = abundance_noOutliers,
  family = tw(link = "log")
)
diagnose_gam(gam_abundance, data = abundance_noOutliers, group_var = "treatment")
results_gam_list[["abundance"]] <- extract_emmeans_contrasts(gam_abundance, "abundance")


# 3. Evenness (Y_zipf)
arkaute_Yzipf <- arkaute %>% 
  filter(!is.na(Y_zipf)) %>% 
  mutate(Y_zipf = Y_zipf * -1)

gam_evenness <- gam(
  Y_zipf ~ treatment + s(plot, bs = "re"),
  data   = arkaute_Yzipf,
  family = gaussian(link = "log")
)
diagnose_gam(gam_evenness, data = arkaute_Yzipf, group_var = "treatment")
results_gam_list[["evenness"]] <- extract_emmeans_contrasts(gam_evenness, "evenness")


# 4. Biomass
arkaute_biomass <- arkaute %>% filter(!is.na(biomass_mice_lm))

gam_biomass <- gam(
  biomass_mice_lm ~ treatment + s(plot, bs = "re"),
  data   = arkaute_biomass,
  family = Gamma(link = "log")
)
diagnose_gam(gam_biomass, data = arkaute_biomass, group_var = "treatment")
results_gam_list[["biomass"]] <- extract_emmeans_contrasts(gam_biomass, "biomass")


# 5. SLA
arkaute_SLA <- arkaute %>% filter(!is.na(SLA))

gam_SLA <- gam(
  SLA ~ treatment + s(plot, bs = "re"),
  data   = arkaute_SLA,
  family = gaussian(link = "identity")
)
diagnose_gam(gam_SLA, data = arkaute_SLA, group_var = "treatment")
results_gam_list[["SLA"]] <- extract_emmeans_contrasts(gam_SLA, "SLA")


# 6. LDMC
arkaute_LDMC <- arkaute %>% filter(!is.na(LDMC))

gam_LDMC <- gam(
  LDMC ~ treatment + s(plot, bs = "re"),
  data   = arkaute_LDMC,
  family = gaussian(link = "identity")
)
diagnose_gam(gam_LDMC, data = arkaute_LDMC, group_var = "treatment")
results_gam_list[["LDMC"]] <- extract_emmeans_contrasts(gam_LDMC, "LDMC")


# 7. leafN
arkaute_leafN <- arkaute %>% filter(!is.na(leafN))

gam_leafN <- gam(
  leafN ~ treatment + s(plot, bs = "re"),
  data   = arkaute_leafN,
  family = gaussian(link = "identity")
)
diagnose_gam(gam_leafN, data = arkaute_leafN, group_var = "treatment")
results_gam_list[["leafN"]] <- extract_emmeans_contrasts(gam_leafN, "leafN")


# =================================================================================
# Tabla final consolidada
# =================================================================================

emmeans_gam_summary_table <- bind_rows(results_gam_list)

# Mostrar la tabla en consola
print(emmeans_gam_summary_table)

# Opcional: exportar la tabla a un archivo CSV
# write.csv(emmeans_summary_table, "results/emmeans_gam_summary.csv", row.names = FALSE)

















##### GLMM ###########



#############
# Helper para diagnósticos de GLMM
#############

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


# =================================================================================
# Lista y función para extraer contrastes de emmeans en GLMM (formato largo)
# =================================================================================
results_glmm_list <- list()

extract_emmeans_glmm <- function(model, var_name) {
  # Obtener contrastes de emmeans
  emm <- emmeans(model, pairwise ~ treatment)
  contrasts_df <- as.data.frame(emm$contrasts)
  
  # Identificar dinámicamente si el test estadístico es t.ratio o z.ratio
  ratio_col <- if ("t.ratio" %in% colnames(contrasts_df)) "t.ratio" else "z.ratio"
  
  # Contrastes objetivo
  target_contrasts <- c("c - p", "c - w", "c - wp", "p - wp")
  
  # Extraer metadatos de la familia y el enlace desde glmmTMB
  fam_info  <- family(model)
  fam_name  <- fam_info$family
  link_name <- fam_info$link
  
  # Filtrar y formatear la tabla larga
  res_df <- contrasts_df %>%
    filter(contrast %in% target_contrasts) %>%
    mutate(
      variable = var_name,
      model    = "GLMM",
      family   = fam_name,
      link     = link_name,
      contrast = gsub(" ", "", contrast), # Limpia los espacios: "c - p" -> "c-p"
      ratio    = .data[[ratio_col]],
      p_value  = p.value
    ) %>%
    select(variable, model, family, link, contrast, ratio, p_value)
  
  return(res_df)
}


# =================================================================================
# Modelado GLMM (glmmTMB)
# =================================================================================

# 1. Richness
glmm_richness <- glmmTMB(
  richness ~ treatment + (1 | plot),
  data   = arkaute,
  family = genpois(link = "log")
)
diagnose_glmm(glmm_richness, data = arkaute, group_var = "treatment")
results_glmm_list[["richness"]] <- extract_emmeans_glmm(glmm_richness, "richness")


# 2. Abundance
abundance_q1  <- quantile(arkaute$abundance, 0.25, na.rm = TRUE)
abundance_q3  <- quantile(arkaute$abundance, 0.75, na.rm = TRUE)
abundance_IQR <- abundance_q3 - abundance_q1

abundance_noOutliers <- arkaute %>%
  select(sampling, treatment, plot, abundance) %>%
  filter(abundance < (abundance_q3 + 1.5 * abundance_IQR))

glmm_abundance <- glmmTMB(
  abundance ~ treatment + (1 | plot),
  data   = abundance_noOutliers,
  family = tweedie(link = "log")
)
diagnose_glmm(glmm_abundance, data = abundance_noOutliers, group_var = "treatment")
results_glmm_list[["abundance"]] <- extract_emmeans_glmm(glmm_abundance, "abundance")


# 3. Evenness (Y_zipf)
arkaute_Yzipf <- arkaute %>% 
  filter(!is.na(Y_zipf)) %>% 
  mutate(Y_zipf = Y_zipf * -1)

glmm_evenness <- glmmTMB(
  Y_zipf ~ treatment + (1 | plot),
  data   = arkaute_Yzipf,
  family = gaussian(link = "log")
)
diagnose_glmm(glmm_evenness, data = arkaute_Yzipf, group_var = "treatment")
results_glmm_list[["evenness"]] <- extract_emmeans_glmm(glmm_evenness, "evenness")


# 4. Biomass
arkaute_biomass <- arkaute %>% filter(!is.na(biomass_mice_lm))

glmm_biomass <- glmmTMB(
  biomass_mice_lm ~ treatment + (1 | plot),
  data   = arkaute_biomass,
  family = Gamma(link = "log")
)
diagnose_glmm(glmm_biomass, data = arkaute_biomass, group_var = "treatment")
results_glmm_list[["biomass"]] <- extract_emmeans_glmm(glmm_biomass, "biomass")


# 5. SLA
arkaute_SLA <- arkaute %>% filter(!is.na(SLA))

glmm_SLA <- glmmTMB(
  SLA ~ treatment + (1 | plot),
  data   = arkaute_SLA,
  family = gaussian(link = "identity")
)
diagnose_glmm(glmm_SLA, data = arkaute_SLA, group_var = "treatment")
results_glmm_list[["SLA"]] <- extract_emmeans_glmm(glmm_SLA, "SLA")


# 6. LDMC
arkaute_LDMC <- arkaute %>% filter(!is.na(LDMC))

glmm_LDMC <- glmmTMB(
  LDMC ~ treatment + (1 | plot),
  data   = arkaute_LDMC,
  family = gaussian(link = "identity")
)
diagnose_glmm(glmm_LDMC, data = arkaute_LDMC, group_var = "treatment")
results_glmm_list[["LDMC"]] <- extract_emmeans_glmm(glmm_LDMC, "LDMC")


# 7. leafN
arkaute_leafN <- arkaute %>% filter(!is.na(leafN))

glmm_leafN <- glmmTMB(
  leafN ~ treatment + (1 | plot),
  data   = arkaute_leafN,
  family = gaussian(link = "identity")
)
diagnose_glmm(glmm_leafN, data = arkaute_leafN, group_var = "treatment")
results_glmm_list[["leafN"]] <- extract_emmeans_glmm(glmm_leafN, "leafN")


# =================================================================================
# Tabla final consolidada
# =================================================================================

emmeans_glmm_summary_table <- bind_rows(results_glmm_list)

# Imprimir la tabla en consola
print(emmeans_glmm_summary_table)

# Opcional: exportar la tabla a un CSV
# write.csv(emmeans_glmm_summary_table, "results/emmeans_glmm_summary.csv", row.names = FALSE)




emmeans_table <- bind_rows(emmeans_gam_summary_table, emmeans_glmm_summary_table)
