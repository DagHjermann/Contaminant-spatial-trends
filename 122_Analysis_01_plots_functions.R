
## Function for ploting GAMM results for each group ----     

# * Parameters with very low mean concentration not shown (mean_value > minvalue)
# * Only graphing parameters with a signinficant geographic pattern (p < 0.05) 
# * Plots 
#    gg1, a line plot (visreg without residuals)
#    gg2, a "fixed" map

plotgamm_group <- function(group = "Organobromines", minvalue = -0.8, 
                           tissue = "Lever",
                           df_map = mapdata,
                           report = FALSE) {
  
  # Get positions for km's to show:
  points <- c(0, 500, 1000, 1500, 2000, 2500, 2685) %>%   
    map_dfr(~get_point_on_coastline(.))
  
  if (!is.na(group)){
    pick_data1 <-plotvalues_gamm_med %>%
      # Pick data
      filter(grepl(group, Substance.Group,ignore.case = TRUE) &
               TISSUE_NAME %in% tissue
      )
    
  } else {
    pick_data1 <- plotvalues_gamm_med %>%
      filter(is.na(Substance.Group) &
               TISSUE_NAME %in% tissue
      )
    
  }
  
  pick_data2 <- pick_data1 %>% 
    # Add 'Position_p' (p-values) to data
    left_join(ttable_gamm %>% select(PARAM, TISSUE_NAME, Position_p), 
              by = c("PARAM", "TISSUE_NAME")) %>%
    mutate(Significant = ifelse(Position_p < 0.05, "Significant", "Not significant")) %>%
    # Find mean value (next 3 lines)
    group_by(PARAM) %>%
    mutate(mean_value = mean(visregFit)) %>%
    ungroup()
  pick_data3 <- pick_data2 %>% 
    # Pick only the most common parameters (next 4 lines)
    filter(mean_value > minvalue & Significant == "Significant") %>%
    # Set factors in order following mean_value
    mutate(PARAM = forcats::fct_reorder(PARAM, mean_value, .desc = TRUE))
  
  if (report){
    cat("===================================================")
    cat(group, "\n")
    cat("---------------------------------------------------")
    print(
      xtabs(~(mean_value > minvalue) + Significant + PARAM, pick_data2)
    )
  }
  non_significant_params <- pick_data2 %>%
    filter(Significant == "Not significant") %>%
    pull(PARAM) %>% 
    unique()
  
  
  plottitle <- case_when(
    tissue == "Lever" ~ paste(group, "in cod liver"),
    tissue == "Muskel" ~ paste(group, "in cod muscle")
  )
  gg1 <- ggplot(pick_data3, aes(Dist_along_coast, visregFit, color = PARAM)) + 
    geom_line(size = rel(1)) + 
    geom_vline(data = points, 
               aes(xintercept = distance),
               linetype = "dashed") +  
    annotate("vline", xintercept = points$distance) +
    labs(x = "Distance along coast", y = "log(Concentration in cod liver)") +
    theme_bw()  +
    labs(title = plottitle)
  
  if (length(non_significant_params) > 0)
    gg1 <- gg1 +
    labs(subtitle = paste("Non-significant parameters (not in graph):", 
                          paste(non_significant_params, collapse = ", ")
    ))
  
  
  # Direction of text labels:
  points$Text_direction <- "West"
  points$Text_direction[c(1,6,7)] <- "East"
  
  # PLot
  gg2 <- ggplot(df_map, aes(x, y)) +
    geom_path() +
    coord_fixed() +
    geom_path(data = coast, color = "red") +
    geom_point(data = points, color = "blue") +
    geom_text(data = points %>% filter(Text_direction == "West"), 
              aes(x = x - 20000, label = paste(distance, "km")), 
              color = "blue", hjust = 1, size = rel(3)) +
    geom_text(data = points %>% filter(Text_direction == "East"), 
              aes(x = x + 20000, label = paste(distance, "km")), 
              color = "blue", hjust = 0, size = rel(3)) +
    expand_limits(x = c(min(df_map$x, na.rm = TRUE) - 150000,
                        max(df_map$x, na.rm = TRUE) + 200000)) +
    theme_minimal() +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank()
    )
  
  cowplot::plot_grid(gg1, gg2, nrow = 1, rel_widths = c(2, 1))
  
}

#
# plot_observations
# For plotting a single parameter
#
plot_observations <- function(param = "BDE47", tissue = "Lever", 
                              data_plot1 = dat2,
                              data_analysis = dat2_notimpacted,
                              tile = TRUE, nrow = 2){
  
  title <- case_when(
    tissue == "Lever" ~ paste0(param, " (liver)"),
    tissue == "Muskel" ~ paste0(param, " (muscle)")
  )
  
  
  if (tile == FALSE){
    gg1 <- data_plot1 %>%
      filter(PARAM %in% param & TISSUE_NAME %in% tissue) %>%
      ggplot(aes(Dist_along_coast, VALUE_WW)) +
      geom_point(alpha = 0.3) +
      stat_summary(fun = "median", colour = "red3", size = 5, geom = "point",
                   shape = 4) +
      labs(title = title) +
      scale_y_log10()
    
  } else {
    df <- data_plot1 %>%
      filter(PARAM %in% param & TISSUE_NAME %in% tissue & !is.na(VALUE_WW)) %>%
      group_by(PARAM, TISSUE_NAME, STATION_CODE, Dist_along_coast, MYEAR) %>%
      summarise(VALUE_WW = median(VALUE_WW), .groups = "drop") %>%
      mutate(
        STATION_CODE = forcats::fct_reorder(STATION_CODE, Dist_along_coast),
        VALUE_WW_c = cut(
          VALUE_WW, 
          breaks = c(0.001,0.003,0.01,0.03,0.1,0.3,1,3,10,30,
                     100,300,1000,3000,10000),
          include.lowest = TRUE)
      )
    
    gg1 <- ggplot(df, aes(STATION_CODE, MYEAR, fill = VALUE_WW_c)) +
      geom_tile() +
      viridis::scale_fill_viridis("WW conc.", discrete = TRUE) +
      scale_y_reverse() +
      labs(title = title)
    
  }
  
  print(gg1)

  modeldata <- try(
    gamm_one_par_one_repl(param, tissue, dat2_notimpacted, return = "model_and_data")
  )

  if (class(modeldata) != "try-error"){
    
    par(mfrow = c(1,2), mar = c(4,5,1,1), oma = c(0,0,2,0))
    visreg(modeldata$mod$gam)
    mtext("Length- and year-adjusted numbers - time (left) and position (right) effects", outer = TRUE)
    
    cat("Significance of year and position effects: \n")
    summary(modeldata$mod$gam)$s.table
    
  }
  
  # gg2 <- ggplot(plotdata$fit, aes(Dist_along_coast, visregFit)) +
  #   geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), fill = "grey70") +
  #   geom_line() +
  #   geom_point(data = plotdata$res, aes(y = visregRes), color = "black", alpha = 0.3) +
  #   labs(title = "Length- and year-adjusted numbers")
  
  # cowplot::plot_grid(gg1, gg2, nrow = nrow)
  
}





