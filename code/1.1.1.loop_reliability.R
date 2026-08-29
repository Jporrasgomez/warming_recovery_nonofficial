
# In this script we check how reliable mice imputations are

# Loop takes a substantial amount of time (~ hours)




stability_list <- list()
reliability_db_list <- list()
reliability_plot_list <- list()

stripplotlist <- list()
densityplotlist <- list()
plotlist <- list()


counter = 0

for (i in c(1:100)){
  
  counter = counter + 1
  
  # Sampling random rows of nind_nona each iteration of the loop
  rows_to_NA <- sample(nrow(nind_nona), perc_NA * nrow(biomass))
  
  # Adding NA to the sampled rows in nind_nona
  mice_check <- nind_nona %>%
    mutate(nind_m2 = replace(nind_m2, rows_to_NA, NA)) %>% 
    #Choosing same variables for modeling 
    select(year, sampling, plot, treatment, code, family, richness, abundance, abundance_community,
           height, Ah, Ab, biomass_i, nind_m2) %>% 
    mutate(across(where(is.character), as.factor))
  
  sum(is.na(mice_check$nind_m2)) == perc_NA*nrow(biomass) 
  
  subset_check <- mice_check %>% 
    filter(is.na(nind_m2)) %>% 
    select(sampling, plot, treatment, code, abundance, abundance_community) %>%  # Keeping only necessary identifiers
    left_join(nind_nona %>% select(sampling, plot, treatment, code, abundance, abundance_community, nind_m2), 
              by = c("sampling", "plot", "treatment", "code", "abundance", "abundance_community")) %>% 
    rename(nind_m2_original = nind_m2)
  
  mice_imputed_check <- mice(mice_check, method = "rf", m = 10, maxit = 300) 
  
  plotlist[[counter]] <- plot(mice_imputed_check) # Check if lines stabilize
  stripplotlist[[counter]] <- stripplot(mice_imputed_check, pch = 20, cex = 1.2) #If imputed values (blue dots) align well with observed values, MICE is working well.
  densityplotlist[[counter]] <- densityplot(mice_imputed_check)
  
  stability_check <- complete(mice_imputed_check, action = "long") %>% 
    as.data.frame() %>% 
    group_by(plot, treatment, sampling, code) %>%
    mutate(nind_m2_imputed = round(mean(nind_m2, 0)),
           sd_imputation = sd(nind_m2)) %>% 
    select(plot, treatment, sampling, code, nind_m2_imputed, sd_imputation) %>%
    distinct()
  
  stability_check$label_imputation <- ifelse(is.na(mice_check$nind_m2), 1, 0)
  
  stability_check <- stability_check %>% 
    filter(label_imputation == 1) %>% 
    mutate(CV = sd_imputation / nind_m2_imputed) %>% 
    as.data.frame()
  
  
  reliability_check <- complete(mice_imputed_check, action = "long") %>%  
    mutate(.imp = as.factor(.imp)) %>%      # Convert to factor
    left_join(subset_check) %>%             # Adding original values of nind_m2
    filter(!is.na(nind_m2_original))   # Keeping only original values
  
  reliability_db <- reliability_check %>% 
    group_by(.imp) %>%   # Group by imputation number
    do({
      model <- lm(nind_m2 ~ nind_m2_original, data = .)  # Fit model
      stats <- glance(model) %>% select(r.squared, p.value)  # Extract R² and p-value
      stats$.imp <- unique(.$.imp)  # Add imputation number
      stats
    }) %>% 
    ungroup() %>% 
    rename(R2 = r.squared, p_value = p.value)
  
  
  stability_list[[counter]] <- stability_check
  reliability_db_list[[counter]] <- reliability_db
  reliability_plot_list[[counter]] <- reliability_check
  
  
}
