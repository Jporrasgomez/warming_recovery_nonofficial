

LRR_dynamics <- function(data, variable){
  
  
  data <- data %>% 
    select(treatment, date, sampling, date_label, plot, all_of(variable)) %>% 
    filter(!is.na(.data[[variable]])) %>% 
    rename(
      value = all_of(variable)
    ) %>% 
    group_by(treatment, sampling) %>% 
    mutate(
      n = n(),
      mean := mean(value),
      sd := sd(value)
    ) %>% 
    mutate(
      variable = variable
    ) %>% 
    as.data.frame() %>% 
    mutate(
      date = ymd(date)
    ) %>% 
    mutate(
      value = ifelse(variable == "Y_zipf", value * -1, value)
    )
  
  
  
  effect <- data %>% 
    select(date,  date_label,  sampling, treatment, mean, sd) %>% 
    distinct()
  
  effect_c <- effect %>% 
    filter(treatment == "c") %>% 
    select(date, sampling, date_label, mean, sd) %>% 
    rename(mean_c = mean,
           sd_c = sd)
  
  effect_wp <- data %>% 
    filter(treatment == "wp") %>% 
    select(date, sampling, date_label, treatment, mean, sd) %>% 
    rename(mean_wp = mean,
           sd_wp = sd) %>% 
    distinct()
  
  effect_w <- data %>% 
    filter(treatment == "w") %>% 
    select(date, sampling, date_label, mean, sd) %>% 
    rename(mean_w = mean,
           sd_w = sd) %>% 
    distinct()
  
  effect_p <- data %>% 
    filter(treatment == "p") %>% 
    select(date, sampling, date_label, mean, sd) %>% 
    rename(mean_p = mean,
           sd_p = sd) %>% 
    distinct()
  
  
  
  
  n = 4
  
  
  ############ LRR = treatment / Control ##########
  
  RR_treat_vs_c <- effect %>% 
    filter(treatment != "c") %>% 
    left_join(effect_c, by = c("date", "sampling", "date_label")) %>% 
    mutate(
      RR = log(mean / mean_c),
      var_RR = (sd^2) / (n * mean^2) + 
        (sd_c^2) / (n * mean_c^2),
      se_RR = sqrt(var_RR)  
    ) %>% 
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
    ) %>% 
    filter(! RR == "-Inf") %>% 
    rename(eff_descriptor = treatment) %>% 
    mutate(
      eff_descriptor = fct_recode(eff_descriptor,
                                 "w_vs_c" = "w",
                                 "p_vs_c" = "p", 
                                 "wp_vs_c" = "wp")) %>% 
    select(
      eff_descriptor, date, sampling, date_label, delta_RR, se_delta_RR)
  
  
  ############ LRR = Combined / Perturbation  ##########
  
  RR_wp_vs_p <- effect_wp %>% 
    left_join(effect_p, by = c("date", "sampling", "date_label")) %>% 
    
    mutate(
      RR = log(mean_wp / mean_p),
      var_RR = (sd_wp^2) / (n * mean_wp^2) + 
        (sd_p^2) / (n * mean_p^2),
      se_RR = sqrt(var_RR)  
      
    ) %>% 
    mutate(
      delta_RR = RR + 0.5 * (
        (sd_wp ^2) / (n * mean_wp^2) - 
          (sd_p^2) / (n * mean_p^2)
      ),
      var_delta_RR = var_RR + 0.5 * (
        (sd_wp^4) / (n^2 * mean_wp^4) + 
          (sd_p^4) / (n^2 * mean_p^4)
      ),
      se_delta_RR = sqrt(var_delta_RR)
    ) %>% 
    filter(!RR == "Inf")%>% 
    filter(!RR == "-Inf") %>% 
    filter(!RR == "NaN") %>% 
    mutate(eff_descriptor = paste0("wp_vs_p")) %>% 
    select(
      eff_descriptor, date, sampling, date_label,  delta_RR, se_delta_RR)
  
  
  
  ############ LRR = Combined / Warming  ##########
  
  RR_wp_vs_w <- effect_wp %>% 
    left_join(effect_w, by = c("date", "sampling", "date_label")) %>% 
    
    mutate(
      RR = log(mean_wp / mean_w),
      var_RR = (sd_wp^2) / (n * mean_wp^2) + 
        (sd_w^2) / (n * mean_w^2),
      se_RR = sqrt(var_RR)  
      
    ) %>% 
    mutate(
      delta_RR = RR + 0.5 * (
        (sd_wp ^2) / (n * mean_wp^2) - 
          (sd_w^2) / (n * mean_w^2)
      ),
      var_delta_RR = var_RR + 0.5 * (
        (sd_wp^4) / (n^2 * mean_wp^4) + 
          (sd_w^4) / (n^2 * mean_w^4)
      ),
      se_delta_RR = sqrt(var_delta_RR)
    ) %>% 
    filter(!RR == "Inf")%>% 
    filter(!RR == "-Inf") %>% 
    filter(!RR == "NaN") %>% 
    mutate(eff_descriptor = paste0("wp_vs_w")) %>% 
    select(
      eff_descriptor, date, sampling, date_label,  delta_RR, se_delta_RR)
  
  
  ###### Binding data #####
  
  
  effsize_dynamics_data <- rbind(RR_treat_vs_c, RR_wp_vs_p) %>% 
    rbind(RR_wp_vs_w) %>% 
    mutate(
      variable = variable, 
      analysis = paste0("delta_LRR")
      ) %>% 
    mutate(
      delta_RR = ifelse(variable == "Y_zipf", delta_RR * -1, delta_RR)   ## If we do not use this, we would show "Unevenness"
    ) %>% 
    filter(delta_RR != "NaN") %>% 
    mutate(variable = variable, 
           analysis = paste0("delta_LRR"), 
           upper_limit = delta_RR + 1.96 * se_delta_RR,
           lower_limit = delta_RR - 1.96 * se_delta_RR) %>% 
    mutate(
      null_effect = ifelse(lower_limit <= 0 & upper_limit >= 0, "YES","NO"),
      scale = (max(abs(upper_limit)) + max(abs(lower_limit)))/100, 
      year = year(date)
    ) %>% 
    rename(
      eff_value = delta_RR
    ) %>% 
    select(
      eff_descriptor, sampling, date, date_label, year, eff_value, lower_limit, upper_limit, null_effect, 
      scale, variable, analysis
    )
  
  
  effsize_dynamics_data <<- effsize_dynamics_data
  
  
  
  
}
