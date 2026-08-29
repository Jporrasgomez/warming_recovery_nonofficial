


rm(list = ls(all.names = TRUE))

# Cargar paquetes
pacman::p_load(
  dplyr, reshape2, tidyverse, lubridate, ggplot2, ggpubr, gridExtra,
  car, ggsignif, dunn.test, rstatix, performance, emmeans, DHARMa, mgcv, glmmTMB
)

# Custom palettes and labels
source("code/palettes_labels.R")

# Global ggplot2 theme
theme_set(
  theme_bw() +
    theme(
      legend.position   = "right",
      panel.grid        = element_blank(),
      strip.background  = element_blank(),
      strip.text        = element_text(face = "bold"),
      text              = element_text(size = 11)
    )
)

# Carga y limpieza de datos
arkaute <- read.csv("data/processed_data/arkaute.csv") %>%
  mutate(
    year      = factor(year),
    date      = ymd(date),
    sampling  = factor(sampling),
    plot      = factor(plot),
    treatment = factor(treatment)
  ) %>%
  filter(sampling != "0") %>%
  filter(!(sampling == "1" & treatment %in% c("p", "wp"))) |> 
  filter(year == 2023)


# =================================================================================
# 1. Preparación de Datasets
# =================================================================================

abundance_q1  <- quantile(arkaute$abundance, 0.25, na.rm = TRUE)
abundance_q3  <- quantile(arkaute$abundance, 0.75, na.rm = TRUE)
abundance_IQR <- abundance_q3 - abundance_q1

datasets <- list(
  richness  = arkaute %>% mutate(richness_int = round(richness)),
  abundance = arkaute %>% filter(abundance < (abundance_q3 + 1.5 * abundance_IQR)) %>% mutate(abundance_int = round(abundance)),
  evenness  = arkaute %>% filter(!is.na(Y_zipf)) %>% mutate(Y_zipf = Y_zipf * -1),
  biomass   = arkaute %>% filter(!is.na(biomass_mice_lm)),
  SLA       = arkaute %>% filter(!is.na(SLA)),
  LDMC      = arkaute %>% filter(!is.na(LDMC)),
  leafN     = arkaute %>% filter(!is.na(leafN))
)

get_var_col <- function(var_key, is_count_family = FALSE) {
  if (var_key == "richness")  return(if (is_count_family) "richness_int" else "richness")
  if (var_key == "abundance") return(if (is_count_family) "abundance_int" else "abundance")
  if (var_key == "evenness")  return("Y_zipf")
  if (var_key == "biomass")   return("biomass_mice_lm")
  return(var_key)
}


# =================================================================================
# 2. Catálogo de Candidatos
# =================================================================================

get_glmm_candidates <- function(var_name) {
  if (var_name %in% c("richness", "abundance")) {
    list(
      list(family = genpois(link = "log"), disp = ~1,         is_count = TRUE),
      list(family = nbinom2(link = "log"), disp = ~1,         is_count = TRUE),
      list(family = tweedie(link = "log"), disp = ~1,         is_count = FALSE),
      list(family = tweedie(link = "log"), disp = ~treatment, is_count = FALSE),
      list(family = nbinom2(link = "log"), disp = ~treatment, is_count = TRUE)
    )
  } else {
    list(
      list(family = gaussian(link = "identity"), disp = ~1,         is_count = FALSE),
      list(family = gaussian(link = "identity"), disp = ~treatment, is_count = FALSE),
      list(family = Gamma(link = "log"),          disp = ~1,         is_count = FALSE),
      list(family = Gamma(link = "log"),          disp = ~treatment, is_count = FALSE),
      list(family = gaussian(link = "log"),       disp = ~treatment, is_count = FALSE)
    )
  }
}

get_gam_candidates <- function(var_name) {
  if (var_name %in% c("richness", "abundance")) {
    list(
      list(type = "standard", family = nb(link = "log"),               is_count = TRUE),
      list(type = "standard", family = Tweedie(p = 1.5, link = "log"), is_count = FALSE),
      list(type = "standard", family = quasipoisson(link = "log"),     is_count = FALSE)
    )
  } else {
    list(
      list(type = "standard", family = gaussian(link = "identity"), is_count = FALSE),
      list(type = "gaulss",   family = gaulss(),                    is_count = FALSE),
      list(type = "standard", family = Gamma(link = "log"),          is_count = FALSE)
    )
  }
}


# =================================================================================
# 3. Función Evaluadora Protegida (Sin errores de car::leveneTest)
# =================================================================================

evaluate_and_extract <- function(model, var_name, model_type, fam_name, link_name, data) {
  
  # 1. Simulación DHARMa protegida
  sim_res <- tryCatch({
    DHARMa::simulateResiduals(fittedModel = model, plot = FALSE)
  }, error = function(e) NULL)
  
  if (is.null(sim_res)) return(NULL)
  
  # Test de Dispersión
  disp_p <- tryCatch(DHARMa::testDispersion(sim_res, plot = FALSE)$p.value, error = function(e) NA)
  
  # Test de Uniformidad / KS
  ks_p   <- tryCatch(DHARMa::testUniformity(sim_res, plot = FALSE)$p.value, error = function(e) NA)
  
  # Test de Levene (Alineado usando el model.frame real extraído del modelo ajustado)
  lev_p <- tryCatch({
    mf <- model.frame(model)
    # Extraemos el vector de treatment limpio que coincidió con el modelo
    grp <- mf$treatment
    resu <- sim_res$scaledResiduals
    
    if (length(resu) == length(grp)) {
      lev_out <- car::leveneTest(resu ~ grp)
      lev_out$`Pr(>F)`[1]
    } else {
      NA
    }
  }, error = function(e) NA)
  
  disp_status <- if (!is.na(disp_p) && disp_p >= 0.05) "NS" else "problem"
  ks_status   <- if (!is.na(ks_p)   && ks_p >= 0.05)   "NS" else "problem"
  lev_status  <- if (!is.na(lev_p)  && lev_p >= 0.05)  "NS" else "problem"
  
  # 2. Emmeans y Contrastes protegidos
  emm <- tryCatch({
    emmeans(model, pairwise ~ treatment)
  }, error = function(e) NULL)
  
  if (is.null(emm)) return(NULL)
  
  contrasts_df <- tryCatch(as.data.frame(emm$contrasts), error = function(e) NULL)
  if (is.null(contrasts_df) || nrow(contrasts_df) == 0) return(NULL)
  
  ratio_col <- if ("t.ratio" %in% colnames(contrasts_df)) {
    "t.ratio"
  } else if ("z.ratio" %in% colnames(contrasts_df)) {
    "z.ratio"
  } else {
    NULL
  }
  
  if (is.null(ratio_col)) return(NULL)
  
  target_cts <- c("c - p", "c - w", "c - wp", "p - wp")
  
  res_df <- contrasts_df %>%
    filter(contrast %in% target_cts) %>%
    mutate(
      variable    = var_name,
      model       = model_type,
      family      = fam_name,
      link        = link_name,
      dispersion  = disp_status,
      dharma_ks   = ks_status,
      levene_test = lev_status,
      contrast    = gsub(" ", "", contrast),
      ratio       = .data[[ratio_col]],
      p_value     = p.value
    )
  
  req_cols <- c("variable", "model", "family", "link", "dispersion", "dharma_ks", "levene_test", "contrast", "ratio", "p_value")
  if (!all(req_cols %in% colnames(res_df))) return(NULL)
  
  res_df <- res_df %>% dplyr::select(dplyr::all_of(req_cols))
  
  num_problems <- sum(c(disp_status, ks_status, lev_status) == "problem")
  
  return(list(res_df = res_df, problems = num_problems))
}
# =================================================================================
# 4. Bucle Automático de Selección
# =================================================================================

all_model_results <- list()

for (var_key in names(datasets)) {
  dat <- datasets[[var_key]]
  
  # ----------------------------------------------------
  # A) GLMM
  # ----------------------------------------------------
  glmm_candidates <- get_glmm_candidates(var_key)
  best_glmm_res   <- NULL
  min_glmm_probs  <- Inf
  
  for (cand in glmm_candidates) {
    col_name  <- get_var_col(var_key, cand$is_count)
    f_formula <- as.formula(paste(col_name, "~ treatment + (1 | plot)"))
    
    fit <- tryCatch({
      suppressWarnings(
        glmmTMB(f_formula, dispformula = cand$disp, data = dat, family = cand$family)
      )
    }, error = function(e) NULL)
    
    if (!is.null(fit) && isTRUE(fit$sdr$pdHess)) {
      fam_info <- family(fit)
      f_name   <- if(cand$disp != ~1) paste0(fam_info$family, " (disp~treat)") else fam_info$family
      
      eval_out <- evaluate_and_extract(fit, var_key, "GLMM", f_name, fam_info$link, dat)
      
      if (!is.null(eval_out) && eval_out$problems < min_glmm_probs) {
        min_glmm_probs <- eval_out$problems
        best_glmm_res  <- eval_out$res_df
      }
      if (min_glmm_probs == 0) break
    }
  }
  all_model_results[[paste0(var_key, "_GLMM")]] <- best_glmm_res
  
  # ----------------------------------------------------
  # B) GAM
  # ----------------------------------------------------
  gam_candidates <- get_gam_candidates(var_key)
  best_gam_res   <- NULL
  min_gam_probs  <- Inf
  
  for (cand in gam_candidates) {
    col_name <- get_var_col(var_key, cand$is_count)
    
    fit <- tryCatch({
      suppressWarnings({
        if (cand$type == "standard") {
          f_formula <- as.formula(paste(col_name, "~ treatment + s(plot, bs = 're')"))
          gam(f_formula, data = dat, family = cand$family)
        } else if (cand$type == "gaulss") {
          f1 <- as.formula(paste(col_name, "~ treatment + s(plot, bs = 're')"))
          gam(list(f1, ~ treatment), data = dat, family = gaulss())
        }
      })
    }, error = function(e) NULL)
    
    if (!is.null(fit)) {
      f_name   <- fit$family$family
      l_name   <- fit$family$link
      eval_out <- evaluate_and_extract(fit, var_key, "GAM", f_name, l_name, dat)
      
      if (!is.null(eval_out) && eval_out$problems < min_gam_probs) {
        min_gam_probs <- eval_out$problems
        best_gam_res  <- eval_out$res_df
      }
      if (min_gam_probs == 0) break
    }
  }
  all_model_results[[paste0(var_key, "_GAM")]] <- best_gam_res
}


# =================================================================================
# 5. Consolidación de Resultados y Exportación
# =================================================================================

emmeans_table <- bind_rows(all_model_results) %>%
  arrange(variable, contrast, model) %>%
  mutate(
    effect_sign = ifelse(ratio > 0, "negative", "positive"), 
    effect_significance = case_when(
      p_value < 0.05                     ~ "significant",
      p_value >= 0.05  & p_value < 0.10  ~ "marginal",
      TRUE                               ~ "non-significant"
    ),
    eff_descriptor = case_when(
      contrast == "c-p"  ~ "p_vs_c", 
      contrast == "c-w"  ~ "w_vs_c", 
      contrast == "c-wp" ~ "wp_vs_c", 
      contrast == "p-wp" ~ "wp_vs_p"
    ), 
    variable = fct_recode(variable, "dominance" = "evenness")
  ) %>% 
  rename(diff_value = ratio) %>% 
  select(-contrast)


# Integración con LRR
lrr_table <- read.csv("results/effect_size_aggregated.csv") %>% 
  select(-scale, -X) %>% 
  rename(diff_value = eff_value) %>% 
  mutate(model = "LRR", variable = as.factor(variable)) %>% 
  filter(!variable %in% c("biomass_raw", "biomass_mice")) %>%
  mutate(
    variable = fct_recode(variable,
                          "evenness" = "Y_zipf",
                          "biomass"  = "biomass_mice_lm"
    ),
    variable = droplevels(variable),
    effect_sign = ifelse(diff_value > 0, "positive", "negative"),
    effect_significance = case_when(
      null_effect == "YES" ~ "non-significant",
      TRUE                 ~ "significant"
    )
  ) %>% 
  select(-upper_limit, -lower_limit, -null_effect)


models <- full_join(emmeans_table, lrr_table) %>% 
  mutate(variable.bis = variable) %>% 
  select(variable, eff_descriptor, model, family, link, dispersion, 
         dharma_ks, levene_test, p_value, effect_significance, diff_value,
         effect_sign, variable.bis)

#models %>% write.csv("results/model_comparison.csv", row.names = FALSE)




models_p_vs_c  <- models |> filter (eff_descriptor == "p_vs_c")
models_w_vs_c  <- models |> filter (eff_descriptor == "w_vs_c")
models_wp_vs_c <- models |> filter (eff_descriptor == "wp_vs_c")
models_wp_vs_p <- models |> filter (eff_descriptor == "wp_vs_p")

palette_sig <- 
  c("significant" = "blue", "marginal" = "orange", "non-significant" = "grey")
palette_shape <- c(
  "positive" = "+",  # o pch 43 / 3
  "negative" = "-"   # o pch 45
)

models |> 
  filter(eff_descriptor == "wp_vs_p") |> 
  ggplot(aes(x = model, y = variable, color = effect_significance, shape = effect_sign)) +
  geom_point(size = 10) +
  scale_color_manual(values = palette_sig) +
  scale_shape_manual(values = palette_shape) +
  labs(title = "2023: Warming effect on recovery (wp vs p)", x = "Statistical analysis", y = NULL,
       shape = "Sign of effect", color = "Effect stat. significance")

models |> 
  filter(eff_descriptor == "p_vs_c") |> 
  ggplot(aes(x = model, y = variable, color = effect_significance, shape = effect_sign)) +
  geom_point(size = 10) +
  scale_color_manual(values = palette_sig) +
  scale_shape_manual(values = palette_shape) +
  labs(title = "2023: Perturbation effect (p vs c)", x = "Statistical analysis", y = NULL,
       shape = "Sign of effect", color = "Effect stat. significance")


models |> 
  filter(eff_descriptor == "wp_vs_c") |> 
  ggplot(aes(x = model, y = variable, color = effect_significance, shape = effect_sign)) +
  geom_point(size = 10) +
  scale_color_manual(values = palette_sig) +
  scale_shape_manual(values = palette_shape) +
  labs(title = "2023: Combined effect (wp vs c)", x = "Statistical analysis", y = NULL,
       shape = "Sign of effect", color = "Effect stat. significance")


models |> 
  filter(eff_descriptor == "w_vs_c") |> 
  ggplot(aes(x = model, y = variable, color = effect_significance, shape = effect_sign)) +
  geom_point(size = 10) +
  scale_color_manual(values = palette_sig) +
  scale_shape_manual(values = palette_shape) +
  labs(title = "2023: Warming effect on assembly (w vs c)", x = "Statistical analysis", y = NULL,
       shape = "Sign of effect", color = "Effect stat. significance")


