



rm(list = ls(all.names = TRUE))  #Se limpia el environment
pacman::p_unload(pacman::p_loaded(), character.only = TRUE) #


pacman::p_load(dplyr,reshape2,tidyverse, lubridate) 

source("code/palettes_labels.R")


# Opening database: 
arkaute <- read.csv("data/processed_data/arkaute.csv") |> 
  mutate(
    year = as.factor(year),
    date = ymd(date),
    sampling = as.factor(sampling),
    plot = as.factor(plot),
    treatment = as.factor(treatment))  

# Database without sampling 0 for aggregated analysis
arkaute_no0 <- arkaute |> 
  filter(sampling != "0")



   #  Log Response Ratio Analysis

# Developping functions for Log Response Ratio analysis #############################################################################


# Aggregated level analysis function

LRR_agg <- function(data, variable){
  
  # Step 1: Average by plot first 
  plot_level <- data |> 
    filter(!is.na(.data[[variable]])) |> 
    group_by(plot, treatment) |> 
    summarise(plot_mean = mean(.data[[variable]]), .groups = "drop")
  
  # Step 2: Compute true treatment-level mean, sd, and plot count (N = 4)
  effect <- plot_level |> 
    group_by(treatment) |> 
    summarise(
      mean = mean(plot_mean), sd = sd(plot_mean), n = n(),
      .groups = "drop"
    )
  
  effect_c <- effect |> 
    ungroup() |>                     
    filter(treatment == "c") |> 
    select(n, mean, sd) |>          
    rename(mean_c = mean, sd_c = sd, n_c = n)
  
  
  effect_p <- effect |> 
    filter(treatment == "p") |> 
    ungroup() |> 
    select(n, treatment, mean, sd) |> 
    rename(mean_p = mean,  sd_p = sd, n_p = n)  |>  
    select(-treatment)
  
  
  effect_w <- effect |> 
    filter(treatment == "w") |> 
    ungroup() |> 
    select(n, treatment, mean, sd) |> 
    rename(mean_w = mean,
           sd_w = sd,
           n_w = n) |> 
    select(-treatment)
  
  effect_wp <- effect |> 
    filter(treatment == "wp") |> 
    ungroup() |> 
    select(treatment, n, mean, sd) |> 
    rename(mean_wp = mean,
           sd_wp = sd, 
           n_wp = n) |> 
    select(-treatment)
  
  
  
  ############ LRR = treatment / Control ##########
  
  RR_treat_vs_c <- effect |> 
    filter(treatment != "c") |> 
    cbind(effect_c) |>
    mutate(
      RR = log(mean / mean_c),
      var_RR = (sd^2) / (n * mean^2) + 
        (sd_c^2) / (n * mean_c^2),
      se_RR = sqrt(var_RR)  
    ) |> 
    mutate(
      delta_RR = RR + 0.5 * (
        (sd^2) / (n * mean^2) - 
          (sd_c^2) / (n * mean_c^2)
      ),
      
      var_delta_RR = var_RR + 0.5 * (
        (sd^4) / (n^2 * mean^4) + 
          (sd_c^4) / (n^2 * mean_c^4)
      ),
      se_delta_RR = sqrt(var_delta_RR)
    ) |> 
    filter(! RR == "-Inf") |> 
    rename(eff_descriptor = treatment) |> 
    mutate(
      eff_descriptor = fct_recode(eff_descriptor,
                                  "w_vs_c" = "w",
                                  "p_vs_c" = "p", 
                                  "wp_vs_c" = "wp")) |> 
    select(
      eff_descriptor, delta_RR, se_delta_RR)
  
  
  ############ LRR = Combined / Perturbation  ##########
  
  RR_wp_vs_p <- effect_wp |> 
    cbind(effect_p) |> 
    mutate(
      RR = log(mean_wp / mean_p),
      var_RR = (sd_wp^2) / (n_wp * mean_wp^2) + 
        (sd_p^2) / (n_p * mean_p^2),
      se_RR = sqrt(var_RR)  
      
    ) |> 
    mutate(
      delta_RR = RR + 0.5 * (
        (sd_wp ^2) / (n_wp * mean_wp^2) - 
          (sd_p^2) / (n_p * mean_p^2)
      ),
      var_delta_RR = var_RR + 0.5 * (
        (sd_wp^4) / (n_wp^2 * mean_wp^4) + 
          (sd_p^4) / (n_p^2 * mean_p^4)
      ),
      se_delta_RR = sqrt(var_delta_RR)
    ) |> 
    filter(!RR == "Inf")|> 
    filter(!RR == "-Inf") |> 
    filter(!RR == "NaN") |> 
    mutate(eff_descriptor = paste0("wp_vs_p")) |> 
    select(
      eff_descriptor, delta_RR, se_delta_RR)
  
  
  
  ###### Binding data #####
  effsize_data <- rbind(RR_treat_vs_c, RR_wp_vs_p) |> 
    rename(
      RR = delta_RR, 
      se_RR = se_delta_RR
    ) |> 
    mutate(variable = variable, 
           analysis = paste0("delta_LRR")) |> 
    mutate(
      RR = ifelse(variable == "Y_zipf", RR * -1, RR)    ## If we do not use this, we would show "Unevenness"
    ) |> 
    mutate(
      upper_limit = RR + 1.96 * se_RR,
      lower_limit = RR - 1.96 * se_RR
    ) |> 
    mutate(
      null_effect = ifelse(lower_limit <= 0 & upper_limit >= 0, "YES","NO"),
      scale = (upper_limit - lower_limit) * 0.3
    ) |> 
    rename(
      eff_value = RR
    ) |> 
    select(
      eff_descriptor, eff_value, lower_limit, upper_limit, null_effect, scale,
      variable, analysis
    )
  
  
  effsize_data <<- effsize_data
  
}


# Temporal dynamics function

LRR_dynamics <- function(data, variable){
  
  data <- data |> 
    select(treatment, date, date_label_noyear, sampling, plot, all_of(variable)) |> 
    filter(!is.na(.data[[variable]])) |> 
    rename(
      value = all_of(variable)
    ) |> 
    group_by(treatment, sampling) |> 
    mutate(
      n = n(),
      mean := mean(value),
      sd := sd(value)
    ) |> 
    mutate(
      variable = variable
    ) |> 
    as.data.frame() |> 
    mutate(
      date = ymd(date)
    ) 
  
  effect <- data |> 
    select(date,  date_label_noyear,  sampling, treatment, n, mean, sd) |> 
    distinct()
  
  effect_c <- effect |> 
    filter(treatment == "c") |> 
    rename(mean_c = mean, sd_c = sd, n_c = n,
    ) |> 
    select(-treatment)
  
  effect_wp <- effect |> 
    filter(treatment == "wp") |> 
    rename(mean_wp = mean, sd_wp = sd, n_wp = n
    )
  
  effect_w <- effect |> 
    filter(treatment == "w") |> 
    rename(mean_w = mean, sd_w = sd, n_w = n
    )
  
  effect_p <- effect |> 
    filter(treatment == "p") |> 
    rename(mean_p = mean, sd_p = sd, n_p = n)

  
  ############ LRR = treatment / Control ##########
  
  RR_treat_vs_c <- effect |> 
    filter(treatment != "c") |> 
    left_join(effect_c, by = c("date", "sampling", "date_label_noyear")) |> 
    mutate(
      
      RR = log(mean / mean_c),
      
      var_RR = (sd^2) / (n * mean^2) + 
        (sd_c^2) / (n_c * mean_c^2),
      
      se_RR = sqrt(var_RR)  
      
    ) |> 
    mutate(
      delta_RR = RR + 0.5 * (
        (sd^2) / (n * mean^2) - 
          (sd_c^2) / (n_c * mean_c^2)
      ),
      var_delta_RR = var_RR + 0.5 * (
        (sd^4) / (n^2 * mean^4) + 
          (sd_c^4) / (n_c^2 * mean_c^4)
      ),
      se_delta_RR = sqrt(var_delta_RR)
    ) |> 
    filter(! RR == "-Inf") |> 
    rename(eff_descriptor = treatment) |> 
    mutate(
      eff_descriptor = fct_recode(eff_descriptor,
                                  "w_vs_c" = "w",
                                  "p_vs_c" = "p", 
                                  "wp_vs_c" = "wp")) |> 
    select(
      eff_descriptor, date, sampling, date_label_noyear, delta_RR, se_delta_RR)
  
  
  ############ LRR = Combined / Perturbation  ##########
  
  RR_wp_vs_p <- effect_wp |> 
    left_join(effect_p, by = c("date", "sampling", "date_label_noyear")) |> 
    mutate(
      RR = log(mean_wp / mean_p),
      var_RR = (sd_wp^2) / (n_wp * mean_wp^2) + 
        (sd_p^2) / (n_p * mean_p^2),
      se_RR = sqrt(var_RR)  
    ) |> 
    mutate(
      delta_RR = RR + 0.5 * (
        (sd_wp ^2) / (n_wp * mean_wp^2) - 
          (sd_p^2) / (n_p * mean_p^2)
      ),
      var_delta_RR = var_RR + 0.5 * (
        (sd_wp^4) / (n_wp^2 * mean_wp^4) + 
          (sd_p^4) / (n_p^2 * mean_p^4)
      ),
      se_delta_RR = sqrt(var_delta_RR)
    ) |> 
    filter(!RR == "Inf")|> 
    filter(!RR == "-Inf") |> 
    filter(!RR == "NaN") |> 
    mutate(eff_descriptor = paste0("wp_vs_p")) |> 
    select(
      eff_descriptor, date, sampling, date_label_noyear,  delta_RR, se_delta_RR)
  
  
  ###### Binding data #####
  effsize_dynamics_data <- rbind(RR_treat_vs_c, RR_wp_vs_p) |> 
    mutate(
      variable = variable, 
      analysis = paste0("delta_LRR")
    ) |> 
    mutate(
      delta_RR = ifelse(variable == "Y_zipf", delta_RR * -1, delta_RR)   ## If we do not use this, we would show "Unevenness"
    ) |> 
    filter(delta_RR != "NaN") |> 
    mutate(variable = variable, 
           analysis = paste0("delta_LRR"), 
           upper_limit = delta_RR + 1.96 * se_delta_RR,
           lower_limit = delta_RR - 1.96 * se_delta_RR) |> 
    mutate(
      null_effect = ifelse(lower_limit <= 0 & upper_limit >= 0, "YES","NO"),
      scale = (max(abs(upper_limit)) + max(abs(lower_limit)))/100, 
      year = year(date)
    ) |> 
    rename(
      eff_value = delta_RR
    ) |> 
    select(
      eff_descriptor, sampling, date, date_label_noyear, year, eff_value, lower_limit, upper_limit, null_effect, 
      scale, variable, analysis
    )
  
  
  effsize_dynamics_data <<- effsize_dynamics_data
  
  
}



# Performing LOG RESPONSE RATIO ANALYSIS ###############################################################################################

  variables <- 
    (arkaute |>  
       select(-date, -year, - date_label, -date_label_noyear, -sampling, -plot, -treatment,
              -OTC, -perturbation, -mean_temperature, -mean_vwc
       ) |> 
       colnames()
    ) 
  
  list_arkaute <- list(arkaute_no0, arkaute)
  list_agg <- list()
  list_dyn <- list()
  
  for(i in seq_along(variables)){ 
    
    
    LRR_agg(list_arkaute[[1]], variables[i])
    list_agg[[i]] <- effsize_data 
    
    LRR_dynamics(list_arkaute[[2]], variables[i])
    list_dyn[[i]] <- effsize_dynamics_data
    
  }
  
  # Aggregated level results
  
  agg <- do.call(rbind, list_agg) |> 
    select(eff_descriptor, variable, scale, eff_value, lower_limit, upper_limit, null_effect)
  
  # Temporal dynamics results
  
  dyn <- do.call(rbind, list_dyn) |> 
    mutate(
      date_label_noyear = factor(
        date_label_noyear,
        levels = unique(date_label_noyear[order(date)]),
        ordered = TRUE
      )
    )
  
  
 
  
# Keeping LRR results
agg |> write.csv("results/LRR_results_agg.csv")
dyn |> write.csv("results/LRR_results_dyn.csv")

