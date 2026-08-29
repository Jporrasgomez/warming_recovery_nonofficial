




rm(list = ls(all.names = TRUE))
pacman::p_load(dplyr, reshape2, tidyverse, lubridate, ggplot2, ggpubr, rpivotTable, ggrepel, here)

source("code/palettes_labels.R")


#In this script, all "1.0.processing_data.R" is repeated, integrating
# 1.1.0.biomass_imputation_level_MICE.R", "5.1.biomass_imputation_level_LM.R"
# and "6.EFFECT_SIZE.R" in order to explore robustness of biomass values to the
# variability of parameter Z (power coefficient of biomass equation) and the
# different imputation levels (MICE and LM) via  sensitivity analyses
# 

k = 2 # Z variability of 25% around 2/3

{
  
z_vector_75 <- c(1/6, 1/3, 1/2, 5/6, 3/3, 7/6)
z_vector_25 <- c(1/2, 5/6)

z_vector_list <- list(z_vector_75, z_vector_25)
z_vector <- z_vector_list[[k]]



order_z_75 <-  c(
  "biomass_z_1_6",
  "biomass_z_1_3",
  "biomass_z_1_2",
  "biomass_z_5_6",
  "biomass_z_1",
  "biomass_z_7_6"
)


order_z_25 <- c(
  "biomass_z_1_2",
  "biomass_z_5_6"
)

order_z_list <- list(order_z_75, order_z_25)
order_z <- order_z_list[[k]]


labels_z = c(
  "z_1_6" = "1/6",
  "z_1_3" = "1/3",
  "z_1_2" = "1/2",
  "z_2_3_original" = "2/3",
  "z_5_6" = "5/6",
  "z_1" = "1",
  "z_7_6" = "7/6"
)




#

 
 flora_raw <- read.csv("data/flora_db_raw.csv")
  
  flora_raw <- flora_raw %>%
    mutate(across(where(is.character), as.factor),
           plot = factor(plot),
           sampling = factor(sampling, levels = sort(unique(flora_raw$sampling))),
           treatment =  factor(flora_raw$treatment, levels =  c("c", "w", "p", "wp")),
           plot = factor(plot, levels = sort(unique(flora_raw$plot)))) %>% 
    select(-date, -category, -OBS, -ni.agrupado, -familia, -species_old_1, -species_old_2)
  
  
  # Adding dates
  sampling_dates <- read.csv("data/sampling_dates.csv") %>% 
    mutate(sampling = as.factor(sampling),
           date = ymd(date), 
           month = month(date, label = TRUE),
           day = day(date), 
           year = year(date)) %>% 
    select(sampling, date, day, month, year, one_month_window, omw_date) %>% 
    mutate(across(where(is.character), as.factor))
  
  
  flora_rare <- flora_raw %>% 
    right_join(sampling_dates, by = join_by(sampling)) %>% 
    select(sampling, one_month_window, omw_date, plot, treatment,
           code, abundance, height, Cb, Db, Cm, Dm, date, month)
  
  
  #Adding species information
  species_code <- read.csv("data/species_code.csv") %>% 
    select(species, code, family, genus_level, species_level, growing_type) %>%
    mutate(across(where(is.character), as.factor))
  
  
  # WE WORK WITH IDENTIFIED SPECIES!!
  
  flora_rare <- merge(flora_rare, species_code, by = "code") 
  
  
  # Modification of 0.1 cm in order to avoid caliber error. 
  
  flora_rare <- flora_rare %>%
    mutate(Dm = coalesce(ifelse(Dm < 0.1, 0.1, Dm), Dm))
  
  flora_rare <- flora_rare %>%
    mutate(Db = coalesce(ifelse(Db < 0.1, 0.1, Db), Db))
  
  
  ## BIOMASS AT INDIVIDUAL LEVEL ####
  
  # Ecuación biomasa####                 
  #Transforming diameters into circumferences
  flora_rare$cm <- round(ifelse(!is.na(flora_rare$Dm), flora_rare$Dm * pi, flora_rare$Cm), 2)
  flora_rare$cb <- round(ifelse(!is.na(flora_rare$Db), flora_rare$Db * pi, flora_rare$Cb), 2)
  
  #Calculating area of circle with circumference values
  flora_rare$Ah <- ((flora_rare$cm)^2)/(4*pi)
  flora_rare$Ab <- ((flora_rare$cb)^2)/(4*pi)
  
  
  taxongroups <- flora_rare %>%
    filter(code %in% c("poaceae", "asteraceae", "tosp", "orchidaceae"))%>%
    group_by(code, sampling, one_month_window, omw_date, plot, treatment, date,
             month, species, family, genus_level, species_level) %>%
    summarize(abundance = sum(unique(abundance), na.rm = T), #here we sum abundances 
              height = mean(height, na.rm = T),
              Ah = mean(Ah, na.rm = T),
              Ab = mean(Ab, na.rm = T))
  
  species <- anti_join(flora_rare, taxongroups, by = "code") %>%
    group_by(code, sampling, one_month_window, omw_date, plot, treatment, date,
             month, species, family, genus_level, species_level) %>%
    summarize(abundance = mean(abundance, na.rm = T), #here we mean abundances (fake mean, species in the same plot and sampling have the same abundance info)
              height = mean(height, na.rm = T),
              Ah = mean(Ah, na.rm = T),
              Ab = mean(Ab, na.rm = T))
  
  flora_medium <- bind_rows(taxongroups, species)
  
  
  d <- 1.96
  z <- 2/3
  flora_medium$x <- (flora_medium$height/2)*(flora_medium$Ab + flora_medium$Ah)
  
  flora_medium$biomass_i <- d*(flora_medium$x^z)

  
  
  # Senstitivity analysis for z 
  

  
  for (i in seq_along(z_vector)) {
    
    # etiqueta estable para el nombre de columna
    z_lbl <- gsub("/", "_", as.character(MASS::fractions(z_vector[i])))  # "2/3" -> "2_3"
    colname <- paste0("biomass_z_", z_lbl)
    
    flora_medium[[colname]] <-  d * (((flora_medium$height/2) * (flora_medium$Ab + flora_medium$Ah))^z_vector[i])
  }
  
  
  
  flora_medium <- flora_medium  %>% 
    rename(original_biomass_i = biomass_i) %>% 
    pivot_longer(
      cols = starts_with("biomass"), 
      names_to = "z", 
      values_to = "biomass_i_z"
    ) 
  
  
  quantile(flora_medium$x, na.rm = T)
  
  library(ggdist)
  raincloud_plot <- ggplot(flora_medium, aes(x = 1, y = log(x))) +
    
    # Half-violin (the "cloud")
    stat_halfeye(
      adjust = 0.5,
      width = 0.6,
      justification = -0.3,
      .width = 0,
      point_colour = NA,
      fill = "gray40"
    ) +
    
    
    # Raw data (the "rain")
    geom_jitter(
      width = 0.08,
      alpha = 0.4,
      size = 2
    ) +
    
    # Boxplot (the "box")
    geom_boxplot(
      width = 0.2,
      outlier.shape = NA,
      alpha = 0.4,
      linewidth = 0.7
    ) +
    
    coord_flip() +
    scale_x_continuous(breaks = NULL) +
    labs(
      x = NULL,
      y = "log(x)"
    ) +
    
    theme1
  
  print(raincloud_plot)
  
  #ggsave("results/Plots/protofinal/x_distribution.png", plot = raincloud_plot, dpi = 300)
  
  
  

  
  flora_medium %>% 
    mutate(z = factor(z, levels = order_z)) %>% 
    ggplot(aes(x = original_biomass_i, y = biomass_i_z)) + 
    facet_wrap(~ z, ncol = 3, nrow = 2, scales = "free") + 
    geom_point(alpha = 0.5) + 
    geom_smooth(method = "lm") +
    theme_bw()
  
  
  # RICHNESSS ########
  
  flora_medium <- flora_medium %>% 
    group_by(plot, sampling) %>% 
    mutate(richness = n_distinct(code, na.rm = T)) %>% 
    ungroup()
  
  # ABUNDANCE ########
  
  flora_medium <- flora_medium %>% 
    group_by(plot, sampling) %>% 
    mutate(abundance_community = sum(abundance, na.rm = T)) %>% 
    ungroup()

  
  ## BIOMASS AT SPECIES LEVEL ###########
  
  
  flora_biomass_raw <- flora_medium %>% 
    filter(!sampling %in% c("0", "1", "2", "12")) %>%  # Remove samplings for which we have no morphological measurements
    filter(!is.na(x)) %>%                              # Remove rows where x is NA since we do not have morph. meas. for those datapoints
    filter(height > 5 & height < 200)               # Remove individuals with height < 5 cm or > 200 cm, because the equation does not properly work with them 
  

## NUMBER OF INDIVIDUALS ####


plots <- read.csv("data/plots.csv") %>% 
  select(plot, treatment_code) %>% 
  rename(treatment = treatment_code)

nind1 <- read.csv("data/n_individuals.csv") %>% 
  merge(plots) %>% 
  select(sampling, treatment,  plot, code, abundance, nind_m2) %>%
  group_by(sampling, treatment, plot, code) %>%
  summarize(nind_m2 = sum(nind_m2), abundance = sum(abundance)) %>% 
  mutate(sampling = as.factor(sampling),
         treatment = as.factor(treatment), 
         plot = as.factor(plot), 
         code = as.factor(code), 
         nind_m2 = as.numeric(nind_m2), 
         abundance = as.numeric(abundance))

nind2 <- flora_raw %>%
  filter(!sampling %in% c( "12", "13", "14", "15", "16", "17", "18", "19", "20"))  %>%
  group_by(sampling, treatment, plot, code) %>%
  mutate(n_individuals = n()) %>%
  group_by(sampling, treatment, plot, code, abundance) %>%
  summarize(nind_m2 = mean(n_individuals, na.rm = T)) %>%
  filter(nind_m2 < 5) 

nind <- bind_rows(nind1, nind2)


# Integrating number of individuals data within biomass information

biomass_nind <- flora_biomass_raw %>% 
  merge(species_code) %>% 
  left_join(nind) %>% 
  mutate(
    year = year(date))


####### CALCULATION OF BIOMASS AT SPECIES LEVEL #############

## CALCULATION OF BIOMASS AT SPECIES LEVEL WITHOUT IMPUTATIONS

biomass_just <- biomass_nind %>% 
  mutate(biomass_s_original = nind_m2 * original_biomass_i,
         biomass_s_z = nind_m2 * biomass_i_z) %>% 
  select(sampling, treatment, plot, year, date,  one_month_window, omw_date, code, species_level, genus_level, family,
         abundance, richness, abundance_community, nind_m2, Ah, Ab, height, biomass_s_original, biomass_s_z, z)




# MICE IMPUTATION

#library(mice)
#
#z_unique <- unique(biomass_nind$z)
#mice_list <- list()
#
#for (i in seq_along(z_unique)){
#  
#data <- biomass_nind %>% 
#  filter(z == z_unique[i])
#  
#biomass_mice0 <- data %>% 
#  select(year, sampling, plot, treatment, code,
#         family, richness, abundance, abundance_community,
#         height, Ah, Ab, biomass_i_z, nind_m2) %>% 
#  mutate(across(where(is.character), as.factor))
#
#biomass_mice_imputed <- mice(biomass_mice0, method = "rf", m = 10, maxit = 300)
#
#
#imputed_db <- complete(biomass_mice_imputed, action = "long") %>% 
#  group_by(plot, treatment, sampling, code) %>% 
#  mutate(nind_m2_imputed = round(mean(nind_m2, 0)),
#         sd_imputation = sd(nind_m2)) %>% 
#  ungroup() %>% 
#  select(plot, treatment, sampling, code, nind_m2_imputed, sd_imputation) %>% 
#  distinct()
#
#imputed_db$label_imputation <- ifelse(is.na(data$nind_m2), 1, 0) # If the value nind is imputed, 1, if not, 0. 
#
#biomass_mice_results <- data %>% 
#  left_join(imputed_db) %>% 
#  mutate(biomass_s_imp = nind_m2_imputed * biomass_i_z) %>% 
#  select(sampling, treatment, plot, year, date,  one_month_window, omw_date, code, species_level, genus_level, family,
#         abundance, richness, abundance_community, biomass_s_imp, nind_m2, nind_m2_imputed, label_imputation)
#
#biomass_mice_results$z <- z_unique[i]
#
#mice_list[[i]] <- biomass_mice_results
#
#}
#
#biomass_mice <- do.call(rbind, mice_list)
#
#biomass_mice <- merge(biomass_mice, biomass_just) %>% select(-Ah, -Ab, -height)
#
#biomass_mice %>% write.csv("data/z_mice_biomass_sensitivity.csv", row.names = F)
#
#
#

biomass <- read.csv("data/processed_data/z_mice_biomass_sensitivity.csv") %>% 
  select(-biomass_s_original) %>% 
  rename(biomass_s_z_imp = biomass_s_imp) %>% 
  filter( z %in% order_z)


# removing outliers of imputation data


result_list <- list()
for(i in seq_along(order_z)){
  
data <- biomass %>% 
  filter(z == order_z[i])

Q1 <- quantile(log(data$biomass_s_z_imp), 0.25, na.rm = TRUE)
Q3 <- quantile(log(data$biomass_s_z_imp), 0.75, na.rm = TRUE)
IQR <- Q3 - Q1

result_list[[i]] <- data %>%
  filter(log(biomass_s_z_imp) >= Q1 - 1.5 * IQR & log(biomass_s_z_imp) <= Q3 + 1.5 * IQR)

}

biomass <- do.call(rbind, result_list)



## BIOMASS AT COMMUNITY LEVEL ##

biomass_community <- biomass %>%
  group_by(plot, sampling, treatment, z) %>%
  mutate(biomass_z_mice = sum(biomass_s_z_imp, na.rm = TRUE),
         biomass_z_raw = sum(biomass_s_z, na.rm = TRUE)) %>%
  ungroup() %>% 
  select(year, date, sampling, one_month_window, omw_date, treatment, plot, code, 
         abundance, richness, abundance_community, z, biomass_z_raw,
         biomass_z_mice) %>% 
  mutate(sampling_date = as.factor(format(ymd(date), "%Y-%m-%d"))) %>% 
  distinct(treatment, plot, sampling, date, omw_date, one_month_window, z, biomass_z_raw,
           biomass_z_mice) %>% 
  mutate(z = as.factor(z))


# FIRST DATABASE for raw data (no imputation).  
biomass_community_raw_wide <- biomass_community %>% 
  select(-biomass_z_mice) %>% 
  mutate(z = fct_relabel(z, ~ sub("^biomass_", "", .x))) %>% 
  pivot_wider(
    names_from = z,
    values_from = biomass_z_raw
  ) %>% 
  mutate(biomass_level = paste0("biomass_raw")) %>% 
  select(-one_month_window, -omw_date)




# SECOND DATABASE. 

biomass_community_mice <- biomass_community %>% 
  select(-biomass_z_raw)

biomass_community_mice_wide <- biomass_community_mice %>% 
  mutate(z = fct_relabel(z, ~ sub("^biomass_", "", .x))) %>% 
  pivot_wider(
    names_from = z,
    values_from = biomass_z_mice
  ) %>% 
  mutate(biomass_level = paste0("biomass_mice")) %>% 
  select(-one_month_window, -omw_date)



## Including samplings 0, 1, 2 and 12 for each z subset


abundance_richness <- read.csv("data/processed_data/abrich_db_plot.csv")


data_lm <- abundance_richness %>%
  full_join(biomass_community_mice,  by = c("sampling", "plot", "treatment", "date", "omw_date", "one_month_window")) %>%
  select(sampling, treatment, plot, z, biomass_z_mice, abundance) %>% 
  rename(biomass = biomass_z_mice) %>% 
  select(sampling, treatment, plot, z, biomass, abundance)




lm_nona <- data_lm %>% 
  filter(!is.na(biomass)) %>%
  mutate(z = factor(z, levels = order_z))

lm_nona %>% 
  ggplot(aes(x = abundance, y = biomass)) + 
  facet_wrap(~z, scales = "free") +
  geom_point(aes(color = treatment)) + 
  geom_smooth(method = "lm", se = TRUE, alpha = 0.05, aes(color = treatment, fill = treatment)) + 
  geom_smooth(method = "lm", se = TRUE, alpha = 0.1, color = "black") + 
  geom_smooth(method = "lm", se = TRUE, alpha = 0.1, color = "black") + 
  scale_color_manual(values = palette_CB, labels = labels) +
  theme_bw()

biomass_na <- data_lm %>% 
  filter(is.na(biomass)) %>% 
  select(-z)



library(broom)

# Function to get estimates

fit_one <- function(df) {
  m <- lm(biomass ~ abundance, data = df)
  tibble(
    R2        = summary(m)$r.squared,
    `p-value` = glance(m)$p.value,                 # p-valor global (F-test del modelo)
    slope     = coef(m)[["abundance"]],
    intercept = coef(m)[["(Intercept)"]]
  )
}

treat_levels <- unique(biomass_community_mice$treatment)
z_levels <- order_z

result_list_regression <- list() 
result_list_biomass_lm <- list() 


for(i in seq_along(order_z)){
  
  data_nona <- lm_nona %>% 
    filter(z == order_z[i])

# LM per treatment
lm_by_trt <- data_nona %>%
  filter(treatment %in% treat_levels) %>%
  group_by(treatment) %>%
  group_modify(~ fit_one(.x)) %>%
  ungroup()

# LM for all treatment together
lm_all_row <- data_nona %>%
  fit_one() %>%
  mutate(treatment = "all") %>%
  select(treatment, everything())

# Results
result <- bind_rows(
  lm_by_trt %>% mutate(treatment = factor(treatment, levels = treat_levels)) %>% arrange(treatment),
  lm_all_row
) %>%
  select(treatment, R2, `p-value`, slope, intercept) %>% 
  mutate(z = z_levels[i])

result_list_regression[[i]] <- result

result_treatments <- result %>% 
  filter(treatment != "all")


lm_fill <- biomass_na %>% 
  left_join(result_treatments) %>% 
  mutate(
    biomass_lm = abundance * slope + intercept, 
    biomass_lm_all = abundance* result$slope[5] + result$intercept[5]
    
  ) %>% 
  select(sampling, treatment, plot, biomass_lm, biomass_lm_all, abundance, z)


# we use the biomass extrapolated by using the linear model for all treatments at once 
lm_fill <- lm_fill %>% 
  select(-biomass_lm) %>% 
  rename(biomass = biomass_lm_all)

biomass_lm <- rbind(data_nona, lm_fill) %>% 
  select(sampling, treatment, plot, biomass, z) %>% 
  mutate(
    biomass = ifelse(biomass < 0, 1, biomass)  # Adding 1 unit of biomass to those estimations under 0 
  ) %>% 
  mutate(
    biomass = ifelse(sampling == "1" & treatment %in% c("p", "wp"),
                             0,
                     biomass)
  ) 

result_list_biomass_lm[[i]] <- biomass_lm

}


results_regressions <- do.call(rbind, result_list_regression)
biomass_community_mice_lm <- do.call(rbind, result_list_biomass_lm)



# THIRD DATABASE

dates_samplings_info <- abundance_richness %>% 
  select(plot, sampling, treatment, date)

biomass_community_mice_lm_wide <- biomass_community_mice_lm %>% 
  mutate(z = fct_relabel(z, ~ sub("^biomass_", "", .x))) %>% 
  mutate(z = fct_relabel(z, ~ sub("_mice$", "", .x))) %>% 
  pivot_wider(
    names_from = z,
    values_from = biomass
  ) %>% 
  mutate(biomass_level = paste0("biomass_mice_lm")) %>%
  left_join(dates_samplings_info)



biomass_sensitivity <- rbind(biomass_community_raw_wide, biomass_community_mice_wide) %>% 
  rbind(biomass_community_mice_lm_wide)


#sensitivity_biomass %>%  write.csv("results/sensitivity_biomass.csv")




### Adding original biomass info


biomass_original_raw <- read.csv("data/processed_data/biomass_no_imputation.csv") %>% 
  rename(z_2_3_original = biomass) %>% 
  mutate(biomass_level = paste0("biomass_raw")) %>%
  select(-omw_date, -one_month_window)

biomass_original_mice <- read.csv("data/processed_data/biomass_mice_imputation.csv") %>% 
  rename(z_2_3_original = biomass) %>% 
  mutate(biomass_level = paste0("biomass_mice")) %>%
  select(-omw_date, -one_month_window)

biomass_original_mice_lm <- read.csv("data/processed_data/biomass_mice_lm.csv") %>%
  rename(z_2_3_original = biomass_mice_lm) %>% 
  mutate(biomass_level = paste0("biomass_mice_lm")) %>% 
  select(-X) %>% 
  left_join(dates_samplings_info) 

biomass_original <- rbind(biomass_original_raw, biomass_original_mice) %>% 
  rbind(biomass_original_mice_lm)






biomass_data_all <- merge(biomass_sensitivity, biomass_original) %>% 
  mutate(treatment = as.factor(treatment))



biomass_levels <- unique(biomass_data_all$biomass_level)
z_levels <- colnames(biomass_data_all)[6:ncol(biomass_data_all)]

biomass_no0 <- biomass_data_all %>% 
  filter(sampling != "0")

source("code/functions/eff_size_LRR_function.R")

list_eff <- list()
counter = 0

for (i in seq_along(biomass_levels)){

  data <- biomass_no0 %>% 
    filter(biomass_level == biomass_levels[i])
  
  for(j in seq_along(z_levels)){
    
    counter = counter + 1
    
    LRR_agg(data, z_levels[j])
    
    list_eff[[counter]] <- effsize_data %>% 
      mutate(biomass_level = paste0(biomass_levels[i]))
    
    rm(effsize_data) 
  }
}

eff_size_agg <- do.call(rbind, list_eff)  %>% 
  mutate(
    eff_value = round(eff_value, 2),
    lower_limit = round(lower_limit, 2),
    upper_limit = round(upper_limit, 2)
  ) %>% 
  rename(z_value = variable) %>% 
  select(eff_descriptor, biomass_level, z_value, eff_value, lower_limit, upper_limit, null_effect)





data <- eff_size_agg %>% 
  filter(eff_descriptor == "wp_vs_p") %>% 
  mutate(
    eff_descriptor = as.factor(eff_descriptor),
    biomass_level = factor(
      biomass_level,
      levels = c("biomass_raw", "biomass_mice", "biomass_mice_lm")
    ),
    z_value = factor(
      z_value,
      levels = c("z_1_6", "z_1_3", "z_1_2", "z_2_3_original", "z_5_6", "z_1", "z_7_6")
    )
  )
 

gg_sensitivity_z_wp <-    
ggplot(data, aes(
  x = z_value,                 # centrado en 0 + pequeño desplazamiento
  y = eff_value,
  color = eff_descriptor
)) +
  facet_wrap( ~biomass_level, scales = "free_y", nrow = 3, ncol = 1,
              strip.position = "left",
              labeller = as_labeller(c(
                biomass_raw     = "No imputation",
                biomass_mice    = "MICE",
                biomass_mice_lm = "MICE + LM"
              ))
              ) +
  
  geom_hline(yintercept = 0, linetype = "dashed",
             color = p_CB, linewidth = 0.5) +
  
  
  geom_linerange(aes(ymin = lower_limit, ymax = upper_limit),
                 linewidth = 1, alpha = 1) +
  
  geom_point(size = 2.5) +
  
  geom_text(aes(
    y = ifelse(eff_value < 0, lower_limit - 0.1, upper_limit + 0.1),
    label = ifelse(null_effect == "NO", "*", NA_character_)
  ),
  show.legend = FALSE,
  size = 8) +
  
  scale_color_manual(values = palette_RR_wp, labels = labels_RR_wp2) +
  
  scale_x_discrete(labels = labels_z) +
  
  labs(x = "z value", y = NULL, color = NULL) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = "black"),
    text = element_text(size = 14),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 14),
    strip.placement = "outside",                
    strip.text.y.left = element_text(            
      angle = 90, face = "bold", size = 14
    ),
    axis.text.y = element_text(hjust = 0.5, face = "plain", size = 12),
    axis.text.x = element_text(face = "plain", size = 12),
    legend.position = "bottom",
    legend.text = element_text(size = 14, face = "plain"),
    legend.key = element_rect(fill = "white", colour = NA)
  )



 
  data_c <- eff_size_agg %>% 
    filter(eff_descriptor %in% c("p_vs_c", "w_vs_c")) %>% 
    mutate(
      eff_descriptor = as.factor(eff_descriptor),
      eff_descriptor = factor(eff_descriptor, levels = c("p_vs_c","w_vs_c")),
      biomass_level = factor(
        biomass_level,
        levels = c("biomass_raw", "biomass_mice", "biomass_mice_lm")
      ),
      z_value = factor(
        z_value,
        levels = c("z_1_6", "z_1_3", "z_1_2", "z_2_3_original", "z_5_6", "z_1", "z_7_6")
      ), 
    ) 
  
  pos_d <- position_dodge(width = 0.5)
  
  
 gg_sensitivity_z_c <-   
  ggplot(data_c, aes(
    x = z_value,
    y = eff_value,
    color = eff_descriptor,
    group = eff_descriptor   # importante para que el dodge respete el orden
  )) +
    facet_wrap(~ biomass_level, scales = "free_y", nrow = 3, ncol = 1,
               strip.position = "left", 
               labeller = as_labeller(c(
                 biomass_raw     = "No imputation",
                 biomass_mice    = "MICE",
                 biomass_mice_lm = "MICE + LM"
               ))
               ) +
    
    geom_hline(yintercept = 0, linetype = "dashed",
               color = "black", linewidth = 0.5) +
    
    geom_linerange(
      aes(ymin = lower_limit, ymax = upper_limit),
      linewidth = 1, alpha = 1,
      position = pos_d
    ) +
    
    geom_point(
      size = 2.5,
      position = pos_d
    ) +
    
    geom_text(
      aes(
        y = ifelse(eff_value < 0, lower_limit - 0.1, upper_limit + 0.1),
        label = ifelse(null_effect == "NO", "*", NA_character_)
      ),
      position = pos_d,
      show.legend = FALSE,
      size = 8
    ) +
   
    
    scale_color_manual(values = palette_RR_CB, labels = labels_RR2) +
    
    scale_x_discrete(labels = labels_z) +
    
    labs(x = "z value", y = NULL, color = NULL) +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", color = "black"),
      text = element_text(size = 14),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 14),
      strip.placement = "outside",                
      strip.text.y.left = element_text(            
        angle = 90, face = "bold", size = 14
      ),
      axis.text.y = element_text(hjust = 0.5, face = "plain", size = 12),
      axis.text.x = element_text(face = "plain", size = 12),
      legend.position = "bottom",
      legend.text = element_text(size = 14, face = "plain"),
      legend.key = element_rect(fill = "white", colour = NA)
    )


}


print(gg_sensitivity_z_wp)


print(gg_sensitivity_z_c)



