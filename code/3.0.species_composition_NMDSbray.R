


rm(list = ls(all.names = TRUE))
pacman::p_load(dplyr, tidyverse, DT, viridis, ggrepel, codyn, vegan, eulerr, ggplot2, ggthemes, ggpubr, ggforce )#


flora_abrich <- read.csv("data/processed_data/flora_abrich.csv")



source("code/palettes_labels.R")
palette <- palette_CB
labels <- labels1



species_ab_sampling <- flora_abrich %>% 
  group_by(code, date, sampling, treatment) %>% 
  summarise(abundance = mean(abundance_s, na.rm = T))

species_ab_sampling <- species_ab_sampling %>% 
  mutate(id = paste0(as.character(treatment), "/" , as.character(sampling)))

species_ab_sampling <- species_ab_sampling %>% 
  filter(!(sampling == "1" & treatment %in% c("p", "wp")))




####### DIFFERENCES AT SAMPLING x TREATMENT LEVEL #############

treats <- unique(flora_abrich$treatment)

list1 <- list()
gglist1 <- list()
gglist2 <- list()
count = 0

for(i in 1:length(treats)){
  
  count = count + 1
  
  list1[[count]] <- subset(species_ab_sampling, treatment == treats[i])
  
  sp_wide <- list1[[count]] %>%
    pivot_wider(id_cols = sampling,
                names_from = code,
                values_from = abundance,
                values_fill = list(abundance = 0)) %>% 
    column_to_rownames(var = "sampling") %>% 
    arrange(as.numeric(rownames(.)))
  
  
  #Relative abundance
  
  # Perform NMDS using Bray-Curtis distance
  nmds_res <- metaMDS(sp_wide, distance = "bray", k = 2, trymax = 250, maxit = 999) 
  
  # Extract NMDS sample scores
  nmds_samples <- as.data.frame(scores(nmds_res, display = "sites"))
  
  # Extract NMDS species scores (optional)
  nmds_species <- as.data.frame(scores(nmds_res, display = "species"))
  
  gglist1[[count]] <- ggplot() +
    geom_text_repel(data = nmds_species %>% 
                      rownames_to_column(var = "sp"),
                    aes(x = NMDS1, y = NMDS2, label = sp),
                    color = "grey",
                    max.overlaps = 30) +
    geom_point(data = nmds_samples %>% 
                 rownames_to_column(var = "sampling"),
               aes(x = NMDS1, y = NMDS2),
               size = 1.5) +
    geom_text_repel(data = nmds_samples %>% 
                      rownames_to_column(var = "sampling"),
                    aes(x = NMDS1, y = NMDS2, label = sampling),
                    max.overlaps = 9) +
    geom_path(data = nmds_samples %>% 
                rownames_to_column(var = "sampling"),
              aes(x = NMDS1, y = NMDS2)) +
    geom_hline(aes(yintercept = 0), color = "gray52", linetype = "dashed") +
    geom_vline(aes(xintercept = 0), color = "gray52", linetype = "dashed") +
    labs(title = paste("NMDS using Bray-Curtis:", treats[i], sep = " "),
         subtitle = paste0("Stress = ", round(nmds_res$stress, 3)),
         x = "NMDS1",
         y = "NMDS2")
  
  
}


ggarrange(
  gglist1[[2]],
  gglist1[[1]],
  gglist1[[3]],
  gglist1[[4]], 
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
  
  sp_wide <- species_ab_sampling %>% 
    filter(sampling == samps[i]) %>% 
    pivot_wider(id_cols = c(sampling, date, treatment, id),
                names_from = code,
                values_from = abundance,
                values_fill = list(abundance = 0)) %>% 
    ungroup() %>% 
    as_tibble()
  
  
  abundance_data <- sp_wide %>%
    select(-treatment, -sampling, -date, -id)
  
    sorensen <- vegdist(abundance_data, method = "bray", binary = TRUE)
    sorensen <- as.matrix(sorensen)  
    rownames(sorensen) <- sp_wide$id
    colnames(sorensen) <- sp_wide$id
    sorensen[upper.tri(sorensen)] <- NA
    
    sorensen_df <- sorensen %>% 
      as.data.frame() %>%
      rownames_to_column(var = "row_name") %>%
      pivot_longer(-row_name, names_to = "col_name", values_to = "value") %>% 
      filter(!is.na(value))
  
  list_sorensen[[i]] <- sorensen
  list_sorensen_df[[i]] <- sorensen_df
  
  
}


print(list_sorensen)


sorensen_df <- bind_rows(list_sorensen_df) %>% 
  filter(!value == "0") %>%
  separate(col_name, into = c("treatment_x", "sampling_x"), sep = "/") %>%
  separate(row_name, into = c("treatment_y", "sampling_y"), sep = "/") %>% 
  select(-sampling_x) %>% 
  rename(sampling = sampling_y) %>% 
  mutate(sampling = factor(as.numeric(sampling), levels = sort(unique(as.numeric(sampling))))) %>% 
  arrange(sampling) %>% 
  mutate(comparison = paste0(treatment_x, "-", treatment_y)) %>% 
  select(-treatment_x, -treatment_y) %>% 
  mutate(comparison = ifelse(comparison == "p-c", "c-p", comparison),
         comparison = ifelse(comparison == "w-c", "c-w", comparison),
         comparison = ifelse(comparison == "wp-c", "c-wp", comparison))





##############################################################################
#### 1. SPECIES COMPOSITION DIFFERENCES AT SAMPLING LEVEL #######################
##############################################################################


########### 1. VISUALIZATION: NMDS ##########################

i = 1
# 1: no abundance transformation
# 2: Hellinger transformation of abundance


{
  sp_wide_sampling <- species_ab_sampling %>%
  pivot_wider(id_cols = c(sampling, date, treatment, id),
              names_from = code,
              values_from = abundance,
              values_fill = list(abundance = 0)) %>% 
  ungroup() %>% 
  as_tibble()

# create a distance matrix 
abundance_data_sampling <- sp_wide_sampling %>% 
  select(-treatment, -sampling, -date, -id)

# Optional step in which I transform abundance data sampling in Hellinger abundances
abundance_data_sampling_hellinger <- vegan::decostand(abundance_data_sampling, method = "hellinger")

abundance_list <- list(abundance_data_sampling, abundance_data_sampling_hellinger)
# Compute Bray-Curtis distance matrix

distance_matrix_sampling_bc <- vegan::vegdist(as.data.frame(abundance_list[[i]]), method = "bray")

# Run NMDS (2 dimensions, 100 tries)
nmds_bc_sampling <- metaMDS(distance_matrix_sampling_bc, k = 2, trymax = 250, maxit = 999)

# Extract NMDS coordinates
nmds_df_sampling <- data.frame(
  NMDS1 = nmds_bc_sampling$points[, 1],
  NMDS2 = nmds_bc_sampling$points[, 2],
  treatment = sp_wide_sampling$treatment,
  sampling = sp_wide_sampling$sampling,
  date = sp_wide_sampling$date
)

# Arrange by sampling order
nmds_df_sampling <- nmds_df_sampling %>% arrange(sampling)


# Species arrows visualization

set.seed(123)  # reproducibility for envfit permutations
fit <- vegan::envfit(nmds_bc_sampling, as.data.frame(abundance_list[[i]]), permutations = 999)

# Extract species scores (vectors) 

sp_scores <- as.data.frame(fit$vectors$arrows) %>%
  mutate(p = fit$vectors$pvals,
         R2 = fit$vectors$r,
         species = rownames(.)) %>%
  filter(p < 0.05, R2 > 0.15) %>%    # Filter significant and high correlation scores                             
  mutate(NMDS1 = NMDS1 * R2 * 1.5,     # Scaling arrows
         NMDS2 = NMDS2 * R2 * 1.5)     # Scaling arrows



# Plot NMDS results using ggplot

ggnmds_alltreatments <- 
  ggplot(nmds_df_sampling, aes(x = NMDS1, y = NMDS2, color = treatment)) +
    
    
    stat_ellipse(geom = "polygon", aes(fill = treatment),
                 alpha = 0.12, show.legend = FALSE, level = 0.68) + 
    
    geom_point(size = 2, show.legend = T) +
    
    geom_text_repel(aes(label = sampling),
                    max.overlaps = 8,
                    size = 4.5,
                    show.legend = F) +
  
    geom_hline(yintercept = 0, color = "gray52", linetype = "dashed") +
    
    geom_path(aes(group = treatment), linewidth = 0.5, alpha = 0.2) +
    
    geom_vline(xintercept = 0, color = "gray52", linetype = "dashed") +
  
    #geom_segment(data = sp_scores,
    #             aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
    #             arrow = arrow(length = unit(0.1, "cm")),
    #             colour = "black", 
    #             alpha = 0.5) +
    
    #geom_text(data = sp_scores,
    #          aes(x = NMDS1, y = NMDS2, label = species),
    #          hjust = 0.5, vjust = -0.3, size = 3.5,
    #          colour = "black",
    #          alpha = 0.5) +
   
    scale_color_manual(values = palette_CB, labels = labels, guide = "legend") +
    
    scale_fill_manual(values = palette_CB, guide = "none" ) +
    
    scale_shape_manual(values = point_shapes, guide = "none") +
    
    labs(
      #title = "NMDS Bray-Curtis: mean abundance of species at sampling level",
         subtitle = paste0("Stress = ", round(nmds_bc_sampling$stress, 3)),
         x = "NMDS1", y = "NMDS2", color = " ") +
    theme1

print(ggnmds_alltreatments) # Supplementary Fig. 1

}







########### 1.2. STATISTICAL ANALYSIS: PERMANOVA ##########################

adonis_sampling <- adonis2(
  distance_matrix_sampling_bc ~ treatment,  # puedes agregar más variables si quieres
  data = sp_wide_sampling,                 # debe tener las variables explicativas
  permutations = 999,                      # número de permutaciones
  method = "bray"
)

# Mostrar resultados
print(adonis_sampling)
# Hay un efecto significativo del tratamiento sobre la composición de especies (p = 0.001). 
# El tratamiento explica aproximadamente el 33.7% de la variación en la composición.

bd <- betadisper(distance_matrix_sampling_bc, sp_wide_sampling$treatment)
anova(bd)
# El resultado de ANOVA para las dispersiónes dentro de grupos (tratamientos) es significativo (p = 0.0004).
# Esto significa que la variabilidad o dispersión dentro de al menos un grupo es diferente respecto a otros grupos.
permutest(bd) 
plot(bd)
boxplot(bd)

permutest(bd, pairwise = TRUE)
TukeyHSD(bd)

library(pairwiseAdonis)

# Ejecutamos las comparaciones por pares
pw_adonis <- pairwise.adonis(
  x           = distance_matrix_sampling_bc,                 # tu matriz de distancias Hellinger–Bray
  factors     = sp_wide_sampling$treatment,  # factor con los cuatro tratamientos
  perm        = 999,                         # número de permutaciones
  p.adjust.m  = "BH"                         # corrección de p por Benjamini–Hochberg
)

print(pw_adonis)







# Data for colonizers

species_ab_plot <- flora_abrich %>% 
  select(date, code, sampling, plot, treatment, family, genus_level, species_level, abundance_s)

species_ab_plot <- species_ab_plot %>% 
  mutate(id = paste0(treatment, "/", sampling, "/", plot))

species_ab_plot <- species_ab_plot %>% 
  filter(!(sampling == "1" & treatment %in% c("p", "wp")))

sp_wide_plot <- species_ab_plot %>%
  pivot_wider(id_cols = c(sampling, date, treatment, plot),
              names_from = code,
              values_from = abundance_s,
              values_fill = list(abundance = 0)) %>%
  mutate(across(everything(), ~ replace_na(., 0)))

sp_wide_plot %>% write.csv("data/processed_data/relative_abudance_species_time.csv")


