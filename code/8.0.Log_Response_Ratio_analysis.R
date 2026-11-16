



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

source("code/functions/LRR_agg.R")
source("code/functions/LRR_dyn.R")

# Performing LOG RESPONSE RATIO ANALYSIS ###############################################################################################

  variables <- 
    (arkaute |>  
       select(-date, -year, - date_label, -date_label_noyear, -sampling, -plot, -treatment,
              -OTC, -TILLING, -mean_temperature, -mean_vwc
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

