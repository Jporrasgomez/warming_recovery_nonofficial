


rm(list = ls(all.names = TRUE))  #Se limpia el environment
pacman::p_unload(pacman::p_loaded(), character.only = TRUE) #se quitan todos los paquetes (limpiamos R)

pacman::p_load(dplyr, reshape2,tidyverse, lubridate, ggplot2,
               ggpubr, gridExtra, stringr, readr, nortest, patchwork)

source("code/palettes_labels.R")


### 1. OPENING DATA #####


data <-  read.csv("data/data_sensors.csv") %>%  select(-X) %>% 
  mutate(
    plot       = factor(plot),
    OTC_label  = factor(OTC_label),
    datetime   = ymd_hm(paste(date, time)),
    date       = as_date(datetime),
    month      = month(datetime, label = TRUE, abbr = TRUE),
    hour       = hour(datetime)
  ) %>% select(-datetime)


data_long <-  read.csv("data/data_sensors_long.csv") %>%  select(-X) %>% 
  mutate(
    plot       = factor(plot),
    OTC_label  = factor(OTC_label),
    datetime   = ymd_hm(paste(date, time)),
    date       = as_date(datetime),
    month      = month(datetime, label = TRUE, abbr = TRUE),
    hour       = hour(datetime)
  )%>% select(-datetime)



## SENSOR DATA FOR DAYS OF SAMPLINGS - for arkaute.csv ##

sampling_dates <- read.csv("data/sampling_dates.csv") %>% 
  mutate(sampling = as.factor(sampling),
         date = ymd(date), 
         month = month(date, label = TRUE),
         day = day(date), 
         year = year(date)) %>% 
  select(sampling, date, day, month, year, one_month_window, omw_date) %>% 
  mutate(across(where(is.character), as.factor))

sampling_dates_vector <- sampling_dates$date

sampling_days_data <- data %>% 
  filter(date %in% sampling_dates_vector) %>% 
  select(t_top, vwc, date, time, plot)  

sampling_days_data$treatment <- gsub("[0-9]", "", as.character(sampling_days_data$plot))
sampling_days_data$plot <- gsub("[^0-9]", "", as.character(sampling_days_data$plot))

sampling_days_data <- sampling_days_data %>% 
  mutate(treatment = as.factor(treatment), 
         plot = as.factor(plot)) %>% 
  group_by(date, plot, treatment) %>% 
  summarize(mean_temperature = mean(t_top), 
            mean_vwc = mean(vwc))

sampling_days_data %>%  write.csv("data/processed_data/temp_vwc_data.csv",  row.names = F)


### 2. OTC EFFECT  ############


daylight_hours <- c(7:19)
growth_season <- c("Apr", "May", "Jun", "Jul", "Aug", "Sep")
summer <- c("Jul", "Aug", "Sep")

data <- data %>% 
  mutate(
    data = "All year - 24h"
  )

data_daylight  <- data %>% 
  filter(hour %in% daylight_hours) %>%   
  mutate(
    data = "All year - daylight"
  )

data_growth <- data %>% 
  filter(month %in% growth_season) %>% 
  mutate(
    data = "Growth season - 24h"
  )

data_growth_daylight <- data %>% 
  filter(month %in% growth_season) %>% 
  filter(hour %in% daylight_hours) %>% 
  mutate(
    data = "Growth season - daylight"
  )


data_summer_24h<- data %>% 
  filter(month %in% summer) %>% 
  mutate(
    data = "Summer - 24h"
  )


data_summer_daylight <- data %>% 
  filter(month %in% summer) %>% 
  filter(hour %in% daylight_hours) %>% 
  mutate(
    data = "Summer - daylight"
  )


data_august_24h<- data %>% 
  filter(month == "Aug") %>% 
  mutate(
    data = "August - 24h"
  )


data_august_daylight <- data %>% 
  filter(month == "Aug") %>% 
  filter(hour %in% daylight_hours) %>% 
  mutate(
    data = "August - daylight"
  )


{data_list <- list()
data_list[["data"]] <- data
data_list[["data_daylight"]] <- data_daylight
data_list[["data_growth"]] <- data_growth
data_list[["data_growth_daylight"]] <- data_growth_daylight
data_list[["data_summer"]] <- data_summer_24h
data_list[["data_summer_daylight"]] <- data_summer_daylight
data_list[["data_august"]] <- data_august_24h
data_list[["data_august_daylight"]] <- data_august_daylight


df_mean_difference <- data.frame(
  data       = rep(NA, 8),
  mean_diff  = rep(NA, 8),
  IC_95      = rep(NA, 8),
  mean_value_OTC = rep(NA, 8),
  sd_value_OTC   = rep(NA, 8),
  mean_value_control   = rep(NA, 8),
  sd_value_control   = rep(NA, 8),
  stringsAsFactors = FALSE
)

counter = 0 

for ( i in 1:8) {
  

  anova_result <- aov(t_top ~ OTC_label, data = data_list[[i]])
  summary(anova_result)
  a <- TukeyHSD(anova_result)
  
  counter = counter + 1
  
  df_mean_difference$data[counter] <- unique(data_list[[i]]$data)
  
  df_mean_difference$mean_diff[counter] <- round(a$OTC_label[1], 2)
  df_mean_difference$IC_95[counter] <- paste0(round(a$OTC_label[2], 2), "-", round(a$OTC_label[3],2))
  
  df_mean_difference$mean_value_OTC[counter] <- mean(data_list[[i]]$t_top[which(data_list[[i]]$OTC_label == "otc")])
  df_mean_difference$sd_value_OTC[counter] <- sd(data_list[[i]]$t_top[which(data_list[[i]]$OTC_label == "otc")])
  
  df_mean_difference$mean_value_control[counter] <-mean(data_list[[i]]$t_top[which(data_list[[i]]$OTC_label == "control")])
  df_mean_difference$sd_value_control[counter] <-sd(data_list[[i]]$t_top[which(data_list[[i]]$OTC_label == "control")])
  
  
}
}

# Seeing anova test for all year, 24h : 
anova_result_temp <- aov(t_top ~ OTC_label, data = data_list[[1]])
summary(anova_result_temp)
TukeyHSD(anova_result_temp)

# Statistical differences in Supplementary Table 4. or: 

#View(df_mean_difference)





### 3. EVOLUTION OF EFFECT THROUGH TIME #########

metadata <- do.call(rbind, data_list)


# Sensors measured 4 variables (See Wild, J. et al. Climate at ecologically 
# relevant scales: A new temperature and soil moisture logger for long-term
# microclimate measurement. Agric. For. Meteorol. 268, 40–47 (2019))

variables = c("t_top", "t_ground", "t_bottom", "vwc")

                # 1        # 2        # 3      # 4

ylabels <- c(
  t_top     = "Temperature at 40 cm (ºC)",
  t_ground  = "Temperature at -6 cm (ºC)",
  t_bottom  = "Temperature at 2 cm (ºC)",
  vwc       = "VWC (%)"
)



{
  

# We use t_top 
i = 1



ytitle      <- ylabels[variables[i]]



gg_boxplot_alldata <- 
  data_long %>% 
  filter(variable == variables[i]) %>% 

  ggboxplot(
    x       = "OTC_label", 
    y       = "value",            
    fill    = "OTC_label",
    width   = 0.6,
    ggtheme = theme2
  ) +
  stat_compare_means(
    method      = "t.test",
    label       = "p.signif",
    comparisons = list(c("control", "otc"))
    )+
  scale_x_discrete(
    name   = NULL,
    labels = labels_OTC
  ) +
  
  scale_fill_manual(
    name   = NULL,
    values = palette_OTC,
    labels = labels_OTC,
    guide  = FALSE
  ) +
  theme1 + 
  theme(legend.position = "none") +
  labs(y = ytitle, x = NULL)



  gg_boxplot_daily_average <- 
    data_long %>% 
    group_by(date, OTC_label, variable) %>% 
    summarise(
      mean_value = round(mean(value, na.rm = TRUE), 4),
      sd_value   = round(sd(value,   na.rm = TRUE), 4)
    ) %>% 
    filter(variable == variables[i]) %>% 
    
    ggboxplot(
      x       = "OTC_label", 
      y       = "mean_value",            
      fill    = "OTC_label",
      width   = 0.6,
      ggtheme = theme2
    ) +
    stat_compare_means(
      method      = "t.test",
      label       = "p.signif",
      comparisons = list(c("control", "otc"))
    ) +
    scale_x_discrete(
      name   = NULL      # <— aquí aplicas tus etiquetas “Without OTC” / “With OTC”
    ) +
    scale_fill_manual(          # Esto sigue afectando sólo a la leyenda (si la tuvieras)
      name   = NULL,
      values = palette_OTC,
      labels = labels_OTC
    ) +
    theme1 + 
    theme(legend.position = "none") +
    labs(y = ytitle, x = NULL)
  
  

# EVOLUTION THROUGH TIME: ALL 652 DAYS MEASURED


gg_allyear <- 
  data_long %>% 
  group_by(date, OTC_label, variable) %>% 
  summarise(mean_value = round(mean(value, na.rm = T), 4),
            sd_value = round(sd(value, na.rm = T), 4)
  )%>% 
  filter(variable == variables[i]) %>% 
  
  ggplot(aes(x = date, y = mean_value, color = OTC_label, fill = OTC_label)) +
  
  geom_line() +
  
  geom_smooth(stat = "smooth", alpha = 0.2) +
  
  scale_x_date(
    date_breaks  = "4 month",            # intervalos de 1 mes
    date_labels  = "%Y-%m-%d",              # ej. "Jan 2024"
    expand       = expansion(add = c(0, 0))  # ajusta márgenes si hace falta
  ) +
  
  scale_color_manual( 
    name = NULL, values = palette_OTC, labels = labels_OTC) +
  
  scale_fill_manual(
    name = NULL, values = palette_OTC, labels = labels_OTC) +
  
  labs ( x = NULL,
         #y = ytitle,
         y = "Temperature (ºC)") +
  
  theme1 




# DAILY AVERAGE DIFFERENCE THROUGH TIME: 24 H

n <- as.numeric(as.Date(max(data$date)) - as.Date(min(data$date)))
# 662 days we sampled

gg_24h_diff <- 
  data %>% 
    group_by(time, OTC_label) %>% 
    summarise(t_top_mean = round(mean(t_top, na.rm = T), 2),
              t_top_sd = round(sd(t_top, na.rm = T), 2), 
              t_bottom_mean = round(mean (t_bottom, na.rm = T), 2),
              t_bottom_sd = round(sd(t_bottom, na.rm = T), 2), 
              t_ground_mean = round(mean(t_ground, na.rm = T), 2),
              t_ground_sd = round(sd(t_ground, na.rm = T), 2),
              vwc_mean = round(mean(vwc, na.rm = T), 6),
              vwc_sd = round(sd(vwc, na.rm = T), 6)
    ) %>%
    pivot_wider(names_from = OTC_label,
                values_from = c(t_top_mean, t_top_sd, t_bottom_mean, t_bottom_sd,
                                t_ground_mean, t_ground_sd, vwc_mean, vwc_sd),
                names_prefix = "") %>%
    select(time, starts_with("t_"), starts_with("vwc_")) %>% 
    mutate(
      t_top_mean_diff = t_top_mean_otc - t_top_mean_control,
      t_top_sd_diff = sqrt((t_top_sd_otc^2 / n) + (t_top_sd_control^2 / n)),  ## equation to calculate SD differences 
      t_bottom_mean_diff = t_bottom_mean_otc - t_bottom_mean_control,
      t_bottom_sd_diff = sqrt((t_bottom_sd_otc^2 / n) + (t_bottom_sd_control^2 / n)),
      t_ground_mean_diff = t_ground_mean_otc - t_ground_mean_control,
      t_ground_sd_diff = sqrt((t_ground_sd_otc^2 / n) + (t_ground_sd_control^2 / n)),
      vwc_mean_diff =  vwc_mean_otc - vwc_mean_control,
      vwc_sd_diff =  sqrt((vwc_sd_otc^2 / n) + (vwc_sd_control^2 / n))
    ) %>% 
    mutate(
      data = names(data_list)[i]
    ) %>% 
    select(ends_with("_diff"), time) %>% 
    pivot_longer(
      cols      = -time,
      names_to  = c("variable", ".value"),
      names_pattern = "(.*)_(mean|sd)_diff"
    ) %>% 
    # renombramos las dos columnas que crea automáticamente: mean y sd
    rename(
      mean_diff_value = mean,
      sd_diff_value   = sd
    ) %>% 
filter(variable == variables[i]) %>%                                       #### Modify in this line the variable or variables we want to see
  ggplot(aes(x = time, y = mean_diff_value, color = variable)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#1FBDC7", linewidth = 1) +
  geom_line(aes(group = variable), color = "#EA6E13", linewidth = 2) +
  geom_line(aes(x = time,
                y = mean_diff_value + sd_diff_value,
                group = variable) , linetype = "dashed", color = "#EA6E13", linewidth = 1) +
  geom_line(aes(x = time,
                y = mean_diff_value - sd_diff_value,
                group = variable) , linetype = "dashed", color = "#EA6E13", linewidth = 1) + 
  scale_x_discrete(breaks = sprintf("%02d:00", c(5, 13, 21))) +
    
  #scale_color_manual(values = palette_sensor) +
    
    
  labs( x = NULL,
        #y = ytitle
        y = "Temperature gain inside OTC (ºC)"
        ) +
  
  theme1 +
  theme(legend.position = "none")



# DAILY AVERAGE DIFFERENCE THROUGH TIME: ALL 652 DAYS MEASURED

gg_year_diff <- 
  data %>% 
  group_by(date, OTC_label) %>% 
  summarise(
    t_top_mean     = round(mean(t_top,     na.rm = TRUE), 2),
    t_bottom_mean  = round(mean(t_bottom,  na.rm = TRUE), 2),
    t_ground_mean  = round(mean(t_ground,  na.rm = TRUE), 2),
    vwc_mean       = round(mean(vwc,       na.rm = TRUE), 6),
    .groups = "drop"
  ) %>%
  
  pivot_wider(
    names_from   = OTC_label,
    values_from  = c(t_top_mean, t_bottom_mean, t_ground_mean, vwc_mean)
  ) %>%
  
  mutate(
    t_top_diff     = t_top_mean_otc    - t_top_mean_control,
    t_bottom_diff  = t_bottom_mean_otc - t_bottom_mean_control,
    t_ground_diff  = t_ground_mean_otc - t_ground_mean_control,
    vwc_diff       = vwc_mean_otc      - vwc_mean_control
  ) %>%
  
  select(date, ends_with("_diff")) %>%
  
  pivot_longer(
    cols       = -date,
    names_to   = "variable",
    values_to  = "mean_diff_value",
    names_pattern = "(.*)_diff"
  ) %>% 
  filter(variable == variables[i]) %>%  
  
  # Plot
  ggplot(aes(x = date, y = mean_diff_value)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#1FBDC7") +
  geom_line(color = "#EA6E13", linewidth = 0.5) +
  geom_smooth(stat = "smooth", alpha = 0.2, color = "#EA6E13", fill = "#EA6E13") +
  labs(
    x = NULL,
    y = "Temperature gain inside OTC (ºC)"
  ) +
  
  scale_x_date(
    date_breaks  = "4 month",            # intervalos de 1 mes
    date_labels  = "%Y-%m-%d",              # ej. "Jan 2024"
    expand       = expansion(add = c(0, 0))  # ajusta márgenes si hace falta
  ) +

  theme1 +
  theme(legend.position = "bottom")

   
}



  print(gg_boxplot_alldata)

  
  print(gg_boxplot_daily_average)
  
  
  print(gg_24h_diff) # Supplementary Fig. 4
  
  
  gg_year_temp <-
    (gg_allyear + gg_year_diff) +
    plot_layout(ncol = 1) + 
    plot_annotation(tag_levels = "a")
  print(gg_year_temp) # Supplementary Fig. 5
  
  
  
  
  
  
## 4. EFFECT OF OTC ON VOLUMETRIC WATER CONTENT


vwc_data <- 
  
  data %>% 
  group_by(date, OTC_label) %>% 
  summarise(t_top_mean = round(mean(t_top, na.rm = T), 2),
            vwc_mean = round(mean(vwc, na.rm = T), 6),
            soil_moisture_mean = round(mean(soil_moisture, na.rm = T), 6))


# Collinearity between transformed VWC and original soil moisture raw data

vwc_data %>% 
  ggplot(aes(y = vwc_mean, x = soil_moisture_mean)) + 
  geom_point(size = 0.5) + 
  labs ( x = "Soil moisture TMS-4 raw signal", y = "VWC (%)") +
  geom_smooth(method = "lm", se = FALSE) +
  
  # primero la ecuación
  stat_regline_equation(
    mapping     = aes(label = after_stat(eq.label)),
    formula     = y ~ x,
    label.x.npc = 0.2,
    label.y.npc = 0.95,
    size        = 5,
    show.legend = FALSE
  ) +
  # luego el R²
  stat_regline_equation(
    mapping     = aes(label = after_stat(rr.label)),
    formula     = y ~ x,
    label.x.npc = 0.2,
    label.y.npc = 0.80,
    size        = 5,
    show.legend = FALSE
  ) +
  
  stat_cor(
    mapping     = aes(label = after_stat(p.label)),
    method      = "pearson",      # test de correlación Pearson
    label.x.npc = 0.2,           # misma X para alinear
    label.y.npc = 0.60,           # un poco más abajo
    size        = 5,
    show.legend = FALSE
  ) +
  
  theme1



gg_vwc_vs_t <- 
  vwc_data %>% 
  ggplot(aes(y = vwc_mean, x = t_top_mean, color = OTC_label, fill = OTC_label)) + 
  geom_point(size = 3, alpha = 0.5) + 
  scale_color_manual( 
    name = NULL,
    values = palette_OTC,
    labels = labels_OTC
  ) +
  scale_fill_manual(
    name = NULL, 
    values = palette_OTC,
    labels = labels_OTC
  ) +
  labs ( x = "Temperature (ºC)", y = "VWC (%)") +
  geom_smooth(method = "lm", se = FALSE, linewidth = 3) +
  
  stat_regline_equation(
    mapping     = aes(label = after_stat(eq.label)),
    formula     = y ~ x,
    label.x.npc = 0.8,
    label.y.npc = 0.95,
    size        = 5,
    show.legend = FALSE
  ) +

  stat_regline_equation(
    mapping     = aes(label = after_stat(rr.label)),
    formula     = y ~ x,
    label.x.npc = 0.8,
    label.y.npc = 0.80,
    size        = 5,
    show.legend = FALSE
  ) +
  
  stat_cor(
    mapping     = aes(label = after_stat(p.label)),
    method      = "pearson",      # test de correlación Pearson
    label.x.npc = 0.8,           # misma X para alinear
    label.y.npc = 0.60,           # un poco más abajo
    size        = 5,
    show.legend = FALSE
  ) +
  
  theme1

  print(gg_vwc_vs_t)  # Supplementary Fig. 6



# Statistics on VWC
  hist(data$vwc)
  anova_result <- aov(vwc ~ OTC_label, data)
  summary(anova_result)
  TukeyHSD(anova_result)











































