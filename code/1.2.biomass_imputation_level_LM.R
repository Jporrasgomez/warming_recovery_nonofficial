



rm(list = ls(all.names = TRUE))
pacman::p_load(dplyr, tidyverse, lubridate, ggplot2, broom)

biomass_mice <- read.csv("data/processed_data/biomass_mice_imputation.csv") %>%  rename(biomass_mice = biomass)
abundance_richness <- read.csv("data/processed_data/abrich_db_plot.csv")

lm0 <- abundance_richness %>%
  full_join(biomass_mice,  by = c("sampling", "plot", "treatment", "date", "omw_date", "one_month_window")) %>%
  select(sampling, treatment, plot, biomass_mice, abundance) %>% 
  rename(biomass = biomass_mice)


lm0 <- lm0 %>% 
  select(sampling, treatment, plot, biomass, abundance)

lm_nona <- lm0 %>% 
  filter(!is.na(biomass))


lm_nona %>% 
  ggplot(aes(x = abundance, y = biomass)) + 
  geom_point(aes(color = treatment)) + 
  geom_smooth(method = "lm", se = TRUE, alpha = 0.05, aes(color = treatment, fill = treatment)) + 
  geom_smooth(method = "lm", se = TRUE, alpha = 0.1, color = "black") + 
  geom_smooth(method = "lm", se = TRUE, alpha = 0.1, color = "black") %>% 
  print()



treat_levels <- c("c", "w", "p", "wp")
dat <- lm_nona

# Function to get estimates

fit_one <- function(df) {
  m <- lm(biomass ~ abundance, data = df)
  tibble(
    R2        = summary(m)$r.squared,
    `p-value` = glance(m)$p.value,                 # global p-value  
    slope     = coef(m)[["abundance"]],
    intercept = coef(m)[["(Intercept)"]]
  )
}


# LM per treatment
by_trt <- dat %>%
  filter(treatment %in% treat_levels) %>%
  group_by(treatment) %>%
  group_modify(~ fit_one(.x)) %>%
  ungroup()

# LM for all treatment together
all_row <- dat %>%
  fit_one() %>%
  mutate(treatment = "all") %>%
  select(treatment, everything())

# Results
result <- bind_rows(
  by_trt %>% mutate(treatment = factor(treatment, levels = treat_levels)) %>% arrange(treatment),
  all_row
) %>%
  select(treatment, R2, `p-value`, slope, intercept)

result_treatments <- result %>% 
  filter(treatment != "all")


lm_fill <- lm0 %>% 
  filter(is.na(biomass)) %>%
  left_join(result_treatments) %>% 
  mutate(
    biomass_lm = abundance * slope + intercept, 
    biomass_lm_all = abundance* result$slope[5] + result$intercept[5]
    
  ) %>% 
  select(sampling, treatment, plot, biomass_lm, biomass_lm_all, abundance)




# we use the biomass extrapolated by using the linear model for all treatments at once 
lm_fill <- lm_fill %>% 
  select(-biomass_lm) %>% 
  rename(biomass = biomass_lm_all)



biomass_lm <- rbind(lm_nona, lm_fill) %>% 
  select(sampling, treatment, plot, biomass) %>% 
  rename(biomass_lm_plot = biomass) %>% 
  mutate(
    biomass_lm_plot = ifelse(biomass_lm_plot < 0, 1, biomass_lm_plot)  # Adding 1 unit of biomass to those estimations under 0 
  ) %>% 
  mutate(
    biomass_lm_plot = ifelse(sampling == "1" & treatment %in% c("p", "wp"),
                              0,
                              biomass_lm_plot)
  ) %>% 
  rename(biomass_mice_lm = biomass_lm_plot)


biomass_lm %>%   write.csv("data/processed_data/biomass_mice_lm.csv")





