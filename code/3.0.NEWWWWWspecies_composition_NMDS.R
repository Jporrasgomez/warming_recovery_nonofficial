





rm(list = ls(all.names = TRUE))
pacman::p_load(dplyr, tidyverse, DT, viridis, ggrepel, codyn, vegan, eulerr, ggplot2, ggthemes, ggpubr, ggforce )#


flora_abrich <- read.csv("data/processed_data/flora_abrich.csv")



source("code/palettes_labels.R")
palette <- palette_CB
labels <- labels1



####### DIFFERENCES AT SAMPLING x TREATMENT LEVEL #############

{
species_ab_sampling <- flora_abrich |> 
  group_by(code, date, sampling, treatment) |> 
  summarise(abundance = mean(abundance_s, na.rm = T),
            .groups = "drop_last") |> 
  ungroup() |> 
  mutate(id = paste0(as.character(treatment), "/" , as.character(sampling))) |> 
  filter(!(sampling == "1" & treatment %in% c("p", "wp")))



treats <- c("c", "w", "p", "wp")
palette_treats <- c("#12D08C", "#E4345A", "#1CC4FD", "#842F8A")
labels_treats <- c("Control", "Warmed-only", "Perturbed-only", "Combined")
  
list1 <- list()
gglist1 <- list()
gglist2 <- list()
count = 0

for(i in 1:length(treats)){
  
  count = count + 1
  
  list1[[count]] <- subset(species_ab_sampling, treatment == treats[i])
  
  sp_wide <- list1[[count]] |>
    pivot_wider(id_cols = sampling,
                names_from = code,
                values_from = abundance,
                values_fill = list(abundance = 0)) |> 
    arrange(as.numeric(sampling)) |> 
    column_to_rownames(var = "sampling")
  
  
  #Relative abundance
  
  # Perform NMDS using Bray-Curtis distance
  nmds_res <- metaMDS(sp_wide, distance = "bray", k = 2, trymax = 250, maxit = 999) 
  
  # Extract NMDS sample scores
  nmds_samples <- as.data.frame(scores(nmds_res, display = "sites"))
  
  # Extract NMDS species scores (optional)
  nmds_species <- as.data.frame(scores(nmds_res, display = "species"))
  
  gglist1[[count]] <- ggplot() +
    geom_text_repel(data = nmds_species |> 
                      rownames_to_column(var = "sp"),
                    aes(x = NMDS1, y = NMDS2, label = sp),
                    color = "grey",
                    max.overlaps = 30) +
    geom_point(data = nmds_samples |> 
                 rownames_to_column(var = "sampling"),
               aes(x = NMDS1, y = NMDS2),
               size = 3,
               color = palette_treats[i]) +
    geom_text_repel(data = nmds_samples |> 
                      rownames_to_column(var = "sampling"),
                    aes(x = NMDS1, y = NMDS2, label = sampling),
                    max.overlaps = 9,
                    color = palette_treats[i]) +
    geom_path(data = nmds_samples |> 
                rownames_to_column(var = "sampling"),
              aes(x = NMDS1, y = NMDS2),
              color = palette_treats[i]) +
    geom_hline(aes(yintercept = 0), color = "gray52", linetype = "dashed") +
    geom_vline(aes(xintercept = 0), color = "gray52", linetype = "dashed") +
    labs(title = paste(labels_treats[i]),
         subtitle = paste0("Stress = ", round(nmds_res$stress, 3)),
         x = "NMDS1",
         y = "NMDS2") + 
    theme1
  
  
}


gg_sp_comp_dynamics <- 
ggarrange(
  gglist1[[1]],
  gglist1[[2]],
  gglist1[[3]],
  gglist1[[4]],
  labels = c("A", "B", "C", "D"),
  ncol = 2, nrow = 2)


#As a rule of thumb literature has identified the following cut-off values for stress-level:
#  
#  Higher than 0.2 is poor (risks for false interpretation).
#0.1 - 0.2 is fair (some distances can be misleading for interpretation).
#0.05 - 0.1 is good (can be confident in inferences from plot).
#Less than 0.05 is excellent (this can be rare).


## SORENSEN 

list_sorensen <- list()
list_sorensen_df <- list()
samps <- sort(unique(flora_abrich$sampling))

for(i in 1:length(samps)){
  
  sp_wide <- species_ab_sampling |> 
    filter(sampling == samps[i]) |> 
    pivot_wider(id_cols = c(sampling, date, treatment, id),
                names_from = code,
                values_from = abundance,
                values_fill = list(abundance = 0)) |> 
    ungroup() |> 
    as_tibble()
  
  
  abundance_data <- sp_wide |>
    select(-treatment, -sampling, -date, -id)
  
  sorensen <- vegdist(abundance_data, method = "bray", binary = TRUE)
  sorensen <- as.matrix(sorensen)  
  rownames(sorensen) <- sp_wide$id
  colnames(sorensen) <- sp_wide$id
  sorensen[upper.tri(sorensen)] <- NA
  
  sorensen_df <- sorensen |> 
    as.data.frame() |>
    rownames_to_column(var = "row_name") |>
    pivot_longer(-row_name, names_to = "col_name", values_to = "value") |> 
    filter(!is.na(value))
  
  list_sorensen[[i]] <- sorensen
  list_sorensen_df[[i]] <- sorensen_df
  
  
}



sorensen_df <- bind_rows(list_sorensen_df) |> 
  filter(!value == "0") |>
  separate(col_name, into = c("treatment_x", "sampling_x"), sep = "/") |>
  separate(row_name, into = c("treatment_y", "sampling_y"), sep = "/") |> 
  select(-sampling_x) |> 
  rename(sampling = sampling_y) |> 
  mutate(sampling = factor(as.numeric(sampling), levels = sort(unique(as.numeric(sampling))))) |> 
  arrange(sampling) |> 
  mutate(comparison = paste0(treatment_x, "-", treatment_y)) |> 
  select(-treatment_x, -treatment_y) |> 
  mutate(comparison = ifelse(comparison == "p-c", "c-p", comparison),
         comparison = ifelse(comparison == "w-c", "c-w", comparison),
         comparison = ifelse(comparison == "wp-c", "c-wp", comparison))





##############################################################################
#### 1. SPECIES COMPOSITION DIFFERENCES AT SAMPLING LEVEL #######################
##############################################################################

# Abundance matrix at treatment level (4 dots : 4 replicates per treatment)

species_ab_treatment <- flora_abrich |> 
  group_by(code, plot, treatment) |> 
  summarise(abundance = mean(abundance_s, na.rm = T),
            .groups = "drop_last") |> 
  ungroup() |> 
  mutate(id_treat_plot = paste0(as.character(treatment), "/" , as.character(plot)))


sp_wide_treatment <- species_ab_treatment |>
  pivot_wider(id_cols = c(treatment, id_treat_plot),
              names_from = code,
              values_from = abundance,
              values_fill = list(abundance = 0)) |> 
  ungroup() |> 
  select(-`NA`) |> 
  as_tibble()

abudance_data_treatment <- sp_wide_treatment |> 
  select(-treatment, -id_treat_plot)

abundance_data_treatment_hellinger <- vegan::decostand(abudance_data_treatment, method = "hellinger")


abundance_treatment_list <- list(abudance_data_treatment, abundance_data_treatment_hellinger)
# Compute Bray-Curtis distance matrix
 i = 1
distance_matrix_treatment_bc <- vegan::vegdist(as.data.frame(abundance_treatment_list[[i]]), method = "bray")

# Run NMDS (2 dimensions, 100 tries)
nmds_bc_treatment <- metaMDS(distance_matrix_treatment_bc, k = 2, trymax = 250, maxit = 999)
nmds_df_treatment <- data.frame(
  NMDS1 = nmds_bc_treatment$points[, 1],
  NMDS2 = nmds_bc_treatment$points[, 2],
  treatment = sp_wide_treatment$treatment
)

gg_sp_comp_treatment <- 
  
  ggplot(nmds_df_treatment, aes(x = NMDS1, y = NMDS2)) +  # <- IMPORTANTE: Quita el color general de aquí
  
  # Stat ellipse
  stat_ellipse(geom = "polygon", aes(fill = treatment, color = treatment),
               alpha = 0.20, show.legend = FALSE, level = 0.95) + 
  

  geom_point(aes(fill = treatment, shape = treatment, color = treatment),
             size = 8, 
             color = "black",
             stroke = 1, 
             show.legend = TRUE) +
  
  # Axis lines
  geom_hline(yintercept = 0, color = "gray52", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "gray52", linetype = "dashed") +
  
  
  # Color and shape codes
  
  scale_fill_manual(values = palette_CB, 
                    labels = c("c" = "Control", "w" = "Warmed-only", 
                               "p" = "Perturbed-only", "wp" = "Combined")) +
  
  scale_color_manual(values = palette_CB,
                     labels = c("c" = "Control", "w" = "Warmed-only", 
                                "p" = "Perturbed-only", "wp" = "Combined")) +
  
  scale_shape_manual(values = c("c"  = 21,  # Circle
                                "w"  = 22,  # Square
                                "p"  = 24,  # Triangle
                                "wp" = 23), # Rhombus
                     labels = c("c" = "Control", "w" = "Warmed-only", 
                                "p" = "Perturbed-only", "wp" = "Combined")) +
  
  labs(
    subtitle = paste0("Stress = ", round(nmds_bc_treatment$stress, 3)),
    x = "NMDS1", y = "NMDS2", 
    fill = NULL, color = NULL, shape = NULL) +
  theme1


}


print(list_sorensen)
#View(sorensen_df)
gg_sp_comp_dynamics
ggsave("results/NMDS_sp_composition_dyn.png", plot = gg_sp_comp_dynamics, dpi = 600)

gg_sp_comp_treatment
ggsave("results/NMDS_sp_composition_agg.png", plot = gg_sp_comp_treatment, dpi = 600)













########### 1.2. STATISTICAL ANALYSIS: PERMANOVA ##########################

adonis_treatment <- adonis2(
  distance_matrix_treatment_bc ~ treatment,  # puedes agregar más variables si quieres
  data = sp_wide_treatment,                 # debe tener las variables explicativas
  permutations = 999,                      # número de permutaciones
  method = "bray"
)

# Mostrar resultados
print(adonis_treatment)
 # R2 reflects how much variance the treatment explain
adonis_treatment$`Pr(>F)`[1]
adonis_treatment$R2[1]

# Dispersion test. It is recommended that dispersion is similar within treatments
# aka, dispersion test is non-significant
bd <- betadisper(distance_matrix_treatment_bc, sp_wide_treatment$treatment)

test_dispersion <- permutest(bd) 
test_dispersion
test_dispersion[[1]]$`Pr(>F)`[1] #p-value is bigger than 0.05. 
plot(bd)
boxplot(bd)

permutest(bd, pairwise = TRUE)
TukeyHSD(bd)  # Differences of dispersion across treatment comparisons 

library(pairwiseAdonis)

# Ejecutamos las comparaciones por pares
pw_adonis <- pairwise.adonis(
  x           = distance_matrix_treatment_bc,                 # tu matriz de distancias Hellinger–Bray
  factors     = sp_wide_treatment$treatment,  # factor con los cuatro tratamientos
  perm        = 999,                         # número de permutaciones
  p.adjust.m  = "BH"                         # corrección de p por Benjamini–Hochberg
)

print(pw_adonis)




# Data for colonizers

species_ab_plot <- flora_abrich |> 
  select(date, code, sampling, plot, treatment, family, genus_level, species_level, abundance_s)

species_ab_plot <- species_ab_plot |> 
  mutate(id = paste0(treatment, "/", sampling, "/", plot))

species_ab_plot <- species_ab_plot |> 
  filter(!(sampling == "1" & treatment %in% c("p", "wp")))

sp_wide_plot <- species_ab_plot |>
  pivot_wider(id_cols = c(sampling, date, treatment, plot),
              names_from = code,
              values_from = abundance_s,
              values_fill = list(abundance = 0)) |>
  mutate(across(everything(), ~ replace_na(., 0)))

sp_wide_plot |> write.csv("data/processed_data/relative_abudance_species_time.csv")








###### Integratting Abundance matrix at sampling level within NMDS at treatment level

{
  sp_wide_sampling <- species_ab_sampling |>
    pivot_wider(id_cols = c(sampling, date, treatment, id),
                names_from = code,
                values_from = abundance,
                values_fill = list(abundance = 0)) |> 
    ungroup() |> 
    as_tibble()
  
  abundance_data_sampling <- sp_wide_sampling |> 
    select(-treatment, -sampling, -date, -id)
  
  # Optional step in which I transform abundance data sampling in Hellinger abundances
  abundance_data_sampling_hellinger <- vegan::decostand(abundance_data_sampling, method = "hellinger")
  
  abundance_list <- list(abundance_data_sampling, abundance_data_sampling_hellinger)
  
  distance_matrix_sampling_bc <- vegan::vegdist(as.data.frame(abundance_list[[i]]), method = "bray")
  
  
  # 1. Alinear las matrices para que tengan exactamente las mismas especies
  all_species <- unique(c(colnames(abudance_data_treatment), colnames(abundance_data_sampling)))
  
  align_columns <- function(df, cols) {
    missing_cols <- setdiff(cols, colnames(df))
    if(length(missing_cols) > 0) {
      df[missing_cols] <- 0
    }
    return(df[, cols]) # Ordena las columnas de forma consistente
  }
  
  abudance_data_treatment_aligned <- align_columns(abudance_data_treatment, all_species)
  abundance_data_sampling_aligned <- align_columns(abundance_data_sampling, all_species)
  
  # (Opcional) Si decides aplicar la transformación Hellinger, hazlo aquí para ambas matrices:
  # abudance_data_treatment_aligned <- vegan::decostand(abudance_data_treatment_aligned, method = "hellinger")
  # abundance_data_sampling_aligned <- vegan::decostand(abundance_data_sampling_aligned, method = "hellinger")
  
  
  # 2. Ejecutar el NMDS base usando la matriz de tratamientos alineada
  distance_matrix_treatment_bc <- vegan::vegdist(abudance_data_treatment_aligned, method = "bray")
  nmds_bc_treatment <- vegan::metaMDS(distance_matrix_treatment_bc, k = 2, trymax = 250, maxit = 999)
  
  
  # -------------------------------------------------------------------------
  # PASO 1: NMDS BASE (Espacial - 16 parcelas)
  # -------------------------------------------------------------------------
  distance_matrix_treatment_bc <- vegan::vegdist(abudance_data_treatment_aligned, method = "bray")
  nmds_bc_treatment <- vegan::metaMDS(distance_matrix_treatment_bc, k = 2, trymax = 250, maxit = 999)
  
  # Dataframe con los 16 puntos grandes y base para las elipses
  nmds_df_treatment <- data.frame(
    NMDS1 = nmds_bc_treatment$points[, 1],
    NMDS2 = nmds_bc_treatment$points[, 2],
    treatment = sp_wide_treatment$treatment,
    id_treat_plot = sp_wide_treatment$id_treat_plot
  )
  
  
  # -------------------------------------------------------------------------
  # PASO 2: PROYECCIÓN TEMPORAL (Temporal - 82 puntos)
  # -------------------------------------------------------------------------
  # A. Calcular coordenadas de las especies cruzando con el NMDS de 16 parcelas
  spp_scores <- as.matrix(vegan::wascores(nmds_bc_treatment$points, abudance_data_treatment_aligned))
  
  # B. Abundancias relativas de los 82 muestreos
  new_abund <- data.matrix(abundance_data_sampling_aligned)
  new_abund_rel <- sweep(new_abund, 1, rowSums(new_abund, na.rm = TRUE), "/")
  new_abund_rel[is.na(new_abund_rel)] <- 0 
  
  # C. Multiplicación matricial (proyección inicial con encogimiento)
  sampling_projected_scores <- new_abund_rel %*% spp_scores
  
  # D. CORRECCIÓN DE ENCOGIMIENTO (Factor de expansión)
  factor_expansion <- sd(nmds_bc_treatment$points) / sd(sampling_projected_scores)
  sampling_projected_scores <- sampling_projected_scores * factor_expansion
  
  # E. Dataframe de los 82 puntos proyectados (líneas y puntos pequeños)
  nmds_df_sampling_projected <- data.frame(
    NMDS1 = sampling_projected_scores[, 1],
    NMDS2 = sampling_projected_scores[, 2],
    treatment = sp_wide_sampling$treatment,
    sampling = sp_wide_sampling$sampling,
    date = sp_wide_sampling$date
  ) |> 
    dplyr::arrange(sampling) # Vital para que las líneas temporales tengan sentido
  
  
  # -------------------------------------------------------------------------
  # PASO 3: GRAFICAR (Integrando ambas escalas)
  # -------------------------------------------------------------------------
  ggplot() +
    
    # --- CAPA 1: ELIPSES ESPACIALES (16 parcelas) ---
    stat_ellipse(data = nmds_df_treatment, 
                 aes(x = NMDS1, y = NMDS2, fill = treatment, color = treatment),
                 geom = "polygon", alpha = 0.12, show.legend = FALSE, level = 0.95) + 
    
    # --- CAPA 2: TRAYECTORIAS TEMPORALES (4 líneas base = treatment) ---
    geom_path(data = nmds_df_sampling_projected,
              aes(x = NMDS1, y = NMDS2, color = treatment, group = treatment),
              linewidth = 0.5, alpha = 0.3) +  
    
    # --- CAPA 3: PUNTOS TEMPORALES (Los 82 muestreos) ---
    geom_point(data = nmds_df_sampling_projected, 
               aes(x = NMDS1, y = NMDS2, color = treatment),
               size = 1.5, alpha = 0.6) +
    
    # --- CAPA 4: ETIQUETAS DE TIEMPO (Solo algunos samplings) ---
    geom_text_repel(data = nmds_df_sampling_projected,
                    aes(x = NMDS1, y = NMDS2, color = treatment, 
                        label = case_when(
                          sampling %in% c(0, 2, 10, 15, 20) ~ as.character(sampling), 
                          TRUE ~ NA_character_
                        )),
                    max.overlaps = 12, size = 3.5, show.legend = FALSE) +
    
    # --- CAPA 5: CENTROIDES REALES (Las 16 parcelas del experimento) ---
    geom_point(data = nmds_df_treatment, 
               aes(x = NMDS1, y = NMDS2, fill = treatment), 
               shape = 21, size = 5, color = "black", stroke = 1, show.legend = TRUE) +
    
    # --- REFERENCIAS Y ESTÉTICA ---
    geom_hline(yintercept = 0, color = "gray52", linetype = "dashed") +
    geom_vline(xintercept = 0, color = "gray52", linetype = "dashed") +
    
    scale_color_manual(values = palette_CB, 
                       labels = c("c" = "Control", "w" = "Warmed-only", 
                                  "p" = "Perturbed-only", "wp" = "Combined")) +
    scale_fill_manual(values = palette_CB, guide = "none") +
    
    labs(
      subtitle = paste0("Stress = ", round(nmds_bc_treatment$stress, 3)),
      x = "NMDS1", y = "NMDS2", color = "Tratamiento"
    ) +
    theme1
  
}




