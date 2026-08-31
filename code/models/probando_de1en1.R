


rm(list = ls(all.names = TRUE))

# Cargar paquetes
pacman::p_load(
  dplyr, reshape2, tidyverse, lubridate, ggplot2, ggpubr, gridExtra,
  car, ggsignif, dunn.test, rstatix, performance, emmeans, DHARMa, mgcv, glmmTMB
)



arkaute <- read.csv("data/processed_data/arkaute.csv") %>%
  mutate(
    year      = factor(year),
    date      = ymd(date),
    sampling  = factor(sampling),
    plot      = factor(plot),
    treatment = factor(treatment)
  ) %>%
  filter(sampling != "0")


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


arkaute_Yzipf   <-
  arkaute |> 
  filter(!is.na(Y_zipf),
         Y_zipf < -1e-5
         ) |> 
  mutate(Y_zipf = Y_zipf*-1)

glmm_evenness <- glmmTMB(Y_zipf ~ treatment + (1 | sampling) + (1 |plot),
                         data = arkaute_Yzipf,
                         family = Gamma(link = "log"))

emmeans(glmm_evenness, pairwise ~ treatment, type = "response")

# Comparar todos contra la referencia "p" únicamente (Ajuste de Dunnett)
emmeans(glmm_evenness, trt.vs.ctrl ~ treatment, ref = "p", type = "response")

# O sin ajuste si es una prueba de hipótesis a priori única:
contrast(emmeans(glmm_evenness, ~ treatment), method = "pairwise", adjust = "none")


diagnose_glmm(glmm_evenness)





arkaute_SLA   <-
  arkaute |> 
  filter(!is.na(SLA))

glmm_SLA <- glmmTMB(SLA ~ treatment + (1 | sampling),
                         data = arkaute_SLA,
                         family = Gamma(link = "log"))

emmeans(glmm_SLA, pairwise ~ treatment, type = "response")
emmeans(glmm_SLA, trt.vs.ctrl ~ treatment, ref = "p", type = "response")
emmeans(glmm_SLA, trt.vs.ctrl ~ treatment, ref = "p", type = "response", adjust = "none")
contrast(emmeans(glmm_SLA, ~ treatment), method = "pairwise", adjust = "none")


diagnose_glmm(glmm_evenness)





# 2. Fit GLMM with AR(1) temporal autocorrelation

arkaute_sla <- arkaute |> filter(!is.na(SLA))
arkaute_sla$treatment  <- factor(arkaute_sla$treatment, levels = c("c", "w", "p", "wp"))
str(arkaute_sla)

glmm_sla <- glmmTMB(
  SLA ~ treatment * sampling + (1 | plot) + ar1(sampling + 0 | plot),
  data = arkaute_sla,
  family = gaussian(link = "identity")
)

# 3. Contrastes globales excluyendo el muestreo 1 (para evitar el problema de nonEst)
samplings_validos <- setdiff(levels(arkaute_sla$sampling), "1")
samplings_validos

emm_trt <- emmeans(glmm_sla, ~ treatment, at = list(sampling = samplings_validos))
pairs(emm_trt, adjust = "tukey")

# 4. Contrastes entre tratamientos en cada fecha de muestreo individual
emm_time <- emmeans(glmm_sla, ~ treatment | sampling_f)

# Contraste específico entre p y wp a lo largo del tiempo
p_vs_wp_time <- contrast(emm_time, method = list("p_vs_wp" = c(c = 0, w = 0, p = 1, wp = -1)))
summary(p_vs_wp_time, infer = c(TRUE, TRUE))
















# 3. Extract key post-hoc contrasts (wp vs. p over time)
emm_time <- emmeans(glmm_sla, ~ treatment | sampling)

# Contrast matrix: c=0, w=0, p=-1, wp=1
wp_vs_p <- contrast(emm_time, method = list("wp_vs_p" = c(0, 0, -1, 1)))
summary(wp_vs_p, adjust = "fdr")



# 1. Calcular las medias marginales solo para la variable tratamiento
emm_general <- emmeans(glmm_sla, ~ treatment)

# 2. Ver TODAS las comparaciones posibles entre los 4 tratamientos (c, w, p, wp)
# Por defecto aplicará una corrección de Tukey para comparaciones múltiples
comparaciones_todas <- pairs(emm_general)
summary(comparaciones_todas)

# 3. Si quieres aplicar SOLO tu contraste específico (wp vs p) a nivel global:
wp_vs_p_general <- contrast(emm_general, method = list("wp_vs_p" = c(0, 0, -1, 1)))
summary(wp_vs_p_general)





