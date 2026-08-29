



ggdyn2 <- function(data, palette, labels, colorline, position, asterisk, caps, breaks_axix_y){
  

pos_dod_c_dyn <- position_dodge(width = 0.5)

levs  <- levels(data$date_label_noyear)
v_perturbation  <- match("May 04", levs) + 0.5 
v_year  <- match("Nov 13", levs) + 0.5 
  
  plot <- 
  ggplot(data, aes(x = date_label_noyear, y = eff_value,
            group = eff_descriptor, color = eff_descriptor)) +
  
  facet_wrap(~ variable, scales = "free_y", ncol = 1, nrow = 9,
             labeller = labeller(eff_descriptor = as_labeller(labels))
  ) +
  
  scale_x_discrete(drop = FALSE) +
    
  #scale_x_discrete( drop = FALSE, labels = data$date_label) + 
  
  geom_hline(yintercept = 0, linetype = "dashed", color = colorline, linewidth = 0.5) +
  
  geom_vline(xintercept = v_perturbation, linetype = "dashed", color = "grey40", linewidth = 0.5) +
    
  geom_vline(xintercept = v_year, linetype = "solid", color = "black", linewidth = 0.2) +
  
  geom_linerange(aes(ymin = lower_limit, ymax = upper_limit),
                 position = pos_dod_c_dyn, alpha = 1, linewidth = 0.8, 
                 show.legend = FALSE) +
  
  geom_point(position = pos_dod_c_dyn, size = 3.5, alpha = 1) +
  
  geom_line(position = pos_dod_c_dyn, linewidth = 0.5) + 
    
    geom_text(aes(
      x = date_label_noyear,
      y = ifelse(eff_value < 0,
                 lower_limit - asterisk * scale,
                 upper_limit + asterisk * scale),
      label = ifelse(null_effect == "NO", "*", NA_character_),
      color = eff_descriptor
    ),
    position      = position,
    inherit.aes   = FALSE,
    size          = 6, 
    show.legend   = FALSE
    ) +
    
    scale_y_continuous(
      breaks      = scales::pretty_breaks(n = breaks_axix_y),
      minor_breaks = NULL,
      expand = expansion(mult = c(0.1, 0.1))
    ) +
    
    
    scale_color_manual(values = palette, labels = labels) +
    
    labs(y = NULL, x = NULL, color = NULL) +
    
    gg_RR_theme 
  
}








