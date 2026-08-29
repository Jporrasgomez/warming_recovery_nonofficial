


rm(list = ls(all.names = TRUE))
pacman::p_load(dplyr, tidyverse, lubridate, ggplot2)


# Merging databases for main results

abundance_richness <- read.csv("data/processed_data/abrich_db_plot.csv")

evenness           <- read.csv("data/processed_data/radcoeff_db_plot.csv") %>%  select(-one_month_window, -omw_date, -year)

cwm_LES            <-  read.csv("data/processed_data/cwm_plot_db.csv") %>%  select(-X, -year)

biomass_just       <- read.csv("data/processed_data/biomass_no_imputation.csv") %>%  rename(biomass_raw = biomass)

biomass_mice       <- read.csv("data/processed_data/biomass_mice_imputation.csv") %>%  rename(biomass_mice = biomass)

temp_vwc_data      <- read.csv("data/processed_data/temp_vwc_data.csv")

biomass_mice_lm    <- read.csv("data/processed_data/biomass_mice_lm.csv")



arkaute <- abundance_richness %>%
  full_join(biomass_just,     by = c("sampling", "plot", "treatment", "date", "omw_date", "one_month_window")) %>%
  full_join(biomass_mice,  by = c("sampling", "plot", "treatment", "date", "omw_date", "one_month_window")) %>%
  full_join(evenness,    by = c("sampling", "plot", "treatment", "date")) %>% 
  full_join(cwm_LES,        by = c("sampling", "plot", "treatment", "date")) %>% 
  full_join(temp_vwc_data,       by = c("plot", "treatment", "date")) %>%
  full_join(biomass_mice_lm,     by = c("sampling", "plot", "treatment")) %>% 
  mutate(
    OTC               = ifelse(treatment %in% c("w", "wp"), paste0("YES"), paste0("NO")),
    perturbation      = ifelse(treatment %in% c("p", "wp"), paste0("YES"), paste0("NO")),
    date              = as.Date(date),
    year              = year(date),
    date_label        = format(date, "%b %d %y"),
    date_label        = factor(date_label, levels = format(sort(unique(date)), "%b %d %y"), ordered = TRUE),
    date_label_noyear = substr(as.character(date_label), 1, nchar(as.character(date_label)) - 3),
    date_label_noyear = factor( date_label_noyear,levels = unique(date_label_noyear[order(date)]),ordered = TRUE)
    ) %>% 
  select(date, omw_date, one_month_window, date_label, date_label_noyear, year, sampling, plot,
       treatment, OTC, perturbation, richness, abundance, biomass_raw, biomass_mice, biomass_mice_lm, Y_zipf, 
       SLA, LDMC, leafN, mean_temperature, mean_vwc) %>% 
  select( 
    -omw_date, -one_month_window) # Removing unnecessary variables

arkaute %>%
  filter(is.na(biomass_raw)) %>% 
  nrow() %>% 
  print()



arkaute %>%  write.csv("data/processed_data/arkaute.csv", row.names = F)







