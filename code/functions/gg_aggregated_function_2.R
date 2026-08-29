




ggagg2 <- function(data, palette, labels, colorline, limitvar, labelvar, breaks_axix_y){
  
 
  data <- data %>% 
   mutate(
      eff_descriptor = factor(eff_descriptor),
      x_jit = (as.integer(eff_descriptor) - 2) 
    )
  
  ggplot(data, aes(
    x = eff_descriptor,                 # centrado en 0 + pequeño desplazamiento
    y = eff_value,
    color = eff_descriptor
  )) +
    facet_grid(rows = vars(variable), scales = "free_y", switch = "y") +
    
    geom_hline(yintercept = 0, linetype = "dashed",
               color = colorline, linewidth = 0.5) +
    
    
    geom_linerange(aes(ymin = lower_limit, ymax = upper_limit),
                   linewidth = 1, alpha = 1) +
    
    geom_point(size = 6) +
    
    geom_text(aes(
      y = ifelse(eff_value < 0, lower_limit - scale, upper_limit + scale),
      label = ifelse(null_effect == "NO", "*", NA_character_)
    ),
    show.legend = FALSE,
    size = 10) +
    
    scale_color_manual(values = palette, labels = labels) +
    
   # Chat gpt help
    #scale_x_continuous(limits = c(-1.8 , 1.8 ), expand = expansion(mult = 0.1)) +
    
    scale_y_continuous( breaks = scales::pretty_breaks(n = breaks_axix_y), 
                        # añade margen relativo por abajo y por arriba (5% y 25% como ejemplo)
                        expand = expansion(mult = c(0.05, 0.1)) ) +
    
    
    labs(x = NULL, y = NULL, color = NULL) +
    
    gg_RR_theme +
    theme(
      strip.background   = element_blank(),
      strip.placement    = "outside",
      strip.text         = element_text(face = "bold", size = 12),
      strip.text.y.left  = element_text(face = "bold"),
      axis.text.y        = element_text(angle = 90, hjust = 0.5, face = "plain", size = 12),
      axis.text.x        = element_blank(),
      axis.ticks.x       = element_blank(),
      legend.position    = "bottom",
      legend.text        = element_text(size = 14, face = "plain")
    )
}




