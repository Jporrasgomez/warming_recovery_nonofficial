



rm(list = ls(all.names = TRUE))  #Se limpia el environment
pacman::p_unload(pacman::p_loaded(), character.only = TRUE) #se quitan todos los paquetes (limpiamos R)


pacman::p_load(dplyr,reshape2,tidyverse, lubridate, ggplot2, ggpubr, gridExtra,
               car, ggsignif, dunn.test, rstatix, ggbreak, effsize) #Cargamos los paquetes que necesitamos

source("code/palettes_labels.R")



arkaute <- read.csv("data/processed_data/arkaute.csv") %>% 
  mutate(
    year = as.factor(year),
    date = ymd(date),
    omw_date = as.factor(omw_date),
    one_month_window = as.factor(one_month_window),
    sampling = as.factor(sampling),
    plot = as.factor(plot),
    treatment = as.factor(treatment))  %>%
  mutate(
    date_label = format(date, "%b %d %y"),  # p.ej. "04-May-2023"
    date_label = factor(
      date_label,
      levels = format(sort(unique(date)), "%b %d %y"),
      ordered = TRUE
    )
  )


arkaute_no0 <- arkaute %>% 
  filter(sampling != "0")



variables <- c("richness",                         # 1     
               "abundance",                        # 2     
               "Y_zipf",                           # 3     
               "biomass_raw",                      # 4     
               "biomass_mice",                     # 5     
               "biomass_mice_lm",                  # 6     
               "SLA",                              # 7     
               "LDMC",                             # 8     
               "leafN"                             # 9     
)



library(metafor)

for(i in seq_along(variables)){

data <- arkaute_no0 %>%
  distinct(sampling, date, plot, treatment, .data[[variables[i]]]) %>%
  group_by(treatment) %>%
  mutate(
    mean_variable = mean(.data[[paste0(variables[i])]], na.rm = T),   # Usamos .data para referirnos a la columna
    sd_variable = sd(.data[[paste0(variables[i])]], na.rm = T),
    n = n()
  ) %>%
  ungroup() %>%
  select(treatment, n, mean_variable, sd_variable) %>%
  distinct()
  

rr_data <- data %>%
  filter(treatment != "c") %>%
  mutate(
    mean_c = data$mean_variable[data$treatment == "c"],
    sd_c   = data$sd_variable[data$treatment == "c"],
    n_c    = data$n[data$treatment == "c"]
  ) %>%
  rename(mean_t = mean_variable, sd_t = sd_variable, n_t = n)



rr_es <- escalc(measure = "ROM",
                m1i = mean_t, sd1i = sd_t, n1i = n_t,
                m2i = mean_c, sd2i = sd_c, n2i = n_c,
                data = rr_data)

rr_analysis <- rma(yi = yi, vi = vi, data = rr_es)
#summary(rr_analysis)

rr_es_df <- as.data.frame(rr_es)




  forest(x = rr_es_df$yi,   # Los valores de los tamaños de efecto
       sei = sqrt(rr_es_df$vi),  # Error estándar (raíz cuadrada de la varianza)
       slab = rr_es_df$treatment,   # Etiquetas de los tratamientos
       xlab = paste0("Effect size (ROM)", " - ",variables[i]))   # Etiqueta del eje X


}
 


