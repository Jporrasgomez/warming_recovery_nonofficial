


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
  filter(sampling != "0") %>%
  filter(!(sampling == "1" & treatment %in% c("p", "wp")))


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



