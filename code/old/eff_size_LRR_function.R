


LRR_agg <- function(data, variable){
  
  
  data <- data %>% 
    select(treatment, all_of(variable)) %>% 
    filter(!is.na(.data[[variable]])) %>% 
    rename(
      value = all_of(variable)
    ) %>% 
    group_by(treatment) %>% 
    mutate(
      n = n(), 
      mean := mean(value),
      sd := sd(value)
    ) %>% 
    mutate(variable = variable) 
  
  
  effect <- data %>% 
    select(treatment, n, mean, sd) %>% 
    distinct()
  
  
  effect_c <- effect %>% 
    filter(treatment == "c") %>% 
    select(treatment, n, mean, sd) %>%  # Mantener treatment
    rename(mean_c = mean,
           sd_c = sd,
           n_c = n)
  
  
  effect_p <- effect %>% 
    filter(treatment == "p") %>% 
    ungroup() %>% 
    select(n, treatment, mean, sd) %>% 
    rename(mean_p = mean,
           sd_p = sd,
           n_p = n) %>% 
    #mutate(eff_descriptor = "wp_vs_p") %>% 
    select(-treatment)
  
  
  effect_w <- effect %>% 
    filter(treatment == "w") %>% 
    ungroup() %>% 
    select(n, treatment, mean, sd) %>% 
    rename(mean_w = mean,
           sd_w = sd,
           n_w = n) %>% 
    #mutate(eff_descriptor = "wp_vs_p") %>% 
    select(-treatment)
  
  effect_wp <- effect %>% 
    filter(treatment == "wp") %>% 
    ungroup() %>% 
    select(treatment, n, mean, sd) %>% 
    rename(mean_wp = mean,
           sd_wp = sd, 
           n_wp = n) %>% 
  select(-treatment)
  
  ######## LRR = treatment / CONTROL ###########
  
  RR_treat_vs_c <- effect %>% 
    filter(treatment != "c") %>% 
    mutate(
      mean_c = effect_c$mean_c,
      sd_c = effect_c$sd_c,
      n_c = effect_c$n_c
    ) %>% 
    mutate(
      RR =  log(mean / mean_c),
      se_RR =  sqrt((sd^2) / (n * mean^2) + 
                      (sd_c^2) / (n_c * mean_c^2))
    ) %>% 
    rename(eff_descriptor = treatment) %>% 
    mutate(
      eff_descriptor = fct_recode(eff_descriptor,
                                 "w_vs_c" = "w",
                                 "p_vs_c" = "p", 
                                 "wp_vs_c" = "wp")) %>% 
    select(eff_descriptor, RR, se_RR)
 
    
  
  ############ LRR = Combined / Perturbation ##########
  
  
  RR_wp_vs_p <- effect_wp %>% 
    cbind(effect_p) %>% 
    mutate(
      RR = log(mean_wp / mean_p),
      se_RR = sqrt((sd_wp^2) / (n_wp * mean_wp^2) + 
                     (sd_p^2) / (n_p * mean_p^2))
    ) %>% 
    mutate(
      eff_descriptor = paste0("wp_vs_p")
      ) %>% 
    filter(!RR == "Inf") %>% 
    select(eff_descriptor, RR, se_RR) 
  
  
  
  ############ LRR = Combined / Warming ##########
  
  
  RR_wp_vs_w <- effect_wp %>% 
    cbind(effect_w) %>% 
    mutate(
      RR = log(mean_wp / mean_w),
      se_RR = sqrt((sd_wp^2) / (n_wp * mean_wp^2) + 
                     (sd_w^2) / (n_w * mean_w^2))
    ) %>% 
    mutate(
      eff_descriptor = paste0("wp_vs_w")
    ) %>% 
    filter(!RR == "Inf") %>% 
    select(eff_descriptor, RR, se_RR) 
  
  
  
  
  ###### Binding data #####
  
  
  effsize_data <- rbind(RR_treat_vs_c, RR_wp_vs_p) %>% 
    rbind(RR_wp_vs_w) %>% 
    mutate(variable = variable, 
           analysis = paste0("LRR")) %>% 
    mutate(
      RR = ifelse(variable == "Y_zipf", RR * -1, RR)    ## If we do not use this, we would show "Unevenness"
    ) %>% 
    mutate(
      upper_limit = RR + 1.96 * se_RR,
      lower_limit = RR - 1.96 * se_RR
      ) %>% 
    mutate(
      null_effect = ifelse(lower_limit <= 0 & upper_limit >= 0, "YES","NO"),
      scale = (upper_limit - lower_limit) * 0.3
    ) %>% 
    rename(
      eff_value = RR
    ) %>% 
    select(
      eff_descriptor, eff_value, lower_limit, upper_limit, null_effect, scale,
      variable, analysis
    )
  

  effsize_data <<- effsize_data

  
  
  
  
  
  
}
