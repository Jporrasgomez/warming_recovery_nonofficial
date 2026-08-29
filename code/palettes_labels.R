



# THEMES




gg_RR_theme <- theme(
  # Panel & grid
  
  panel.grid        = element_blank(),
  panel.background = element_rect(fill = "white", color = "black"),
  
  # Text (global)
  text              = element_text(size = 14),
  
  # Facets / strips
  strip.background  = element_blank(),
  strip.text        = element_text(face = "bold", size = 14),
  strip.text.x      = element_blank(),
  strip.text.y      = element_blank(),
  
  # Axes
  axis.text.y       = element_text(
    angle = 90,
    hjust = 0.5,
    face  = "plain",
    size  = 12
  ),
  axis.text.x       = element_text(
    angle = 45,
    hjust = 1,
    face  = "plain",
    size  = 12
  ),
  
  # Legend
  legend.position   = "bottom",
  legend.text       = element_text(size = 14, face = "plain"), 
  legend.key = element_rect(fill = "white", colour = NA)
)



theme1 <- 
  theme(
    # Panel & grid
    panel.grid        = element_blank(),
    panel.background = element_rect(fill = "white", color = "black"),
    
    # Text (global)
    text              = element_text(size = 14),
    
    # Facets / strips
    strip.background  = element_blank(),
    strip.text        = element_text(face = "bold", size = 14),
    
    # Axes
    axis.text.y       = element_text(
      angle = 90,
      hjust = 0.5,
      face  = "plain",
      size  = 14
    ),
    axis.text.x       = element_text(
      angle = 0,
      hjust = 0.5,
      face  = "plain",
      size  = 14
    ),
    
    # Legend
    legend.position   = "bottom",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text       = element_text(size = 14, face = "plain"),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_rect(fill = "white", colour = NA)
  )



theme2 <- 
theme(
  legend.position = "bottom",
  axis.text.x     = element_blank(),
  axis.ticks.x    = element_blank(),
  panel.grid = element_blank(),
  strip.background = element_blank(),
  strip.text = element_text(face = "bold"),
  text = element_text(size = 15)
)

theme3 <- 
  theme(
    legend.position = "bottom",
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    text = element_text(size = 15)
  )


## Variable names

limits_variables <- c("richness",
                      "abundance",
                      "Y_zipf",
                      #"biomass",
                      #"biomass012",
                      "biomass_lm_plot",
                      "NMDS1",
                      "NMDS2",
                      "PC1",
                      "PC2")

labels_variables <- c("richness" = "Richness",
                      "abundance" = "Cover",
                      "Y_zipf" = "Evenness",
                      #"biomass" = "Biomass",
                      #"biomass012" = "Biomass",
                      "biomass_lm_plot" = "Biomass", 
                      "NMDS1" = "SC1",
                      "NMDS2" = "SC2",
                      "PC1" = "FT1", 
                      "PC2" = "FT2")




# SENSORS 

palette_sensor <- c("t_top" = "red", "t_bottom" = "purple", "t_ground" = "orange", "vwc" = "blue3")

palette_OTC = c("control" = "#1FBDC7", "otc"     = "#EA6E13")

labels_OTC = c("control" = "Ambient", "otc" = "OTC")

labels_variables_sensor <- c("t_top" = "Temp. 40 cm", "t_bottom"= "Temp. 2 cm", "t_ground" = "Temp. -6 cm", "vwc" = "VWC")





# TREATMENTS

c_CB = "#12D08C" # "CB" color blind
w_CB = "#E4345A"
p_CB = "#1CC4FD"
wp_CB = "#842F8A"



labels1 <- c("c" = "Control", "w" = "Warming", "p" = "Perturbation", "wp" = "Combined")


palette_CB <- c("c" = c_CB, "w" = w_CB, "p" = p_CB, "wp" = wp_CB)

## Response ratio
labels_RR2 <- c("w_vs_c" = "Warming effect on assembly",
                "p_vs_c" = "Recovery of ambient communities",
                "wp_vs_c" = "Combined effect on assembly")

palette_RR_CB <- c("w_vs_c" = w_CB, "p_vs_c" = p_CB, "wp_vs_c" = wp_CB)


palette_RR_wp <- c(wp_vs_p = "#903996", wp_vs_w = "#903996")

labels_RR_wp <- c("wp_vs_w" = "Global Change vs Warming", "wp_vs_p" = "Warming effect on recovery")

labels_RR_wp2 <- c("wp_vs_w" = "Perturbation effect", "wp_vs_p" = "Warming effect on recovery")
labels_RR_wp3 <- c("wp_vs_w" = "Global Change vs Warming", "wp_vs_p" = "Combined vs Perturbation")




palette3 <- c("c" = "#1F8A8C", "w" = "#D93232", "p" = "#F4A300", "wp" = "#3A3A3A")


palette8 <- c("c" = "#00C4A7", "w" = "#A238A2", "p" = "#F4A300", "wp" = "#3A3A3A")




point_shapes <- c("c" = 16, "w" = 17, "p" = 15, "wp" = 18)




