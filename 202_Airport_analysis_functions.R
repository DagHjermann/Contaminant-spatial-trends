
#
# Make new variable 'VALUE_WW_r' with random values between
# 0.5*LOQ and LOQ for data undcer LOQ
#
# NOTE: this is a midfied version using either 
# - use_sample_LOQ = TRUE for sample-level LOQ
# - use_sample_LOQ = FALSE for year/tissue-level LOQ (given by LOQ variable, hich is computed in
#   Jupyterhub )
# 
#
# Returns a new data set (the difference with input data set is the addition of a 'VALUE_WW_r' variable)
#
add_random_data <- function(data, use_sample_LOQ = TRUE){
  
  # Set todata frame (not tibble)
  data <- as.data.frame(data)
  
  # Under LOQs
  sel <- !is.na(data$FLAG1)
  
  # Test that runif works for vactors of min, max:
  # runif(2, min = c(1,50), max = c(10,100))
  data$VALUE_WW_r <- data$VALUE_WW
  
  if (use_sample_LOQ){
    data$VALUE_WW_r[sel] <- runif(
      n = sum(sel), 
      min = data$VALUE_WW[sel]/2,
      max = data$VALUE_WW[sel]
    ) %>% round(4)
    
  } else {
    data$VALUE_WW_r[sel] <- runif(
      n = sum(sel), 
      min = data$LOQ[sel]/2,
      max = data$VALUE_WW[sel]
    ) %>% round(4)
    
  }
  
  data
  
}


# Formats p-values
format_p <- function(p){
  case_when(
    p < 0.001 ~ "< 0.001",
    p < 0.01 ~ round(p, 3) %>% as.character(),
    p <= 1 ~ round(p, 2) %>% as.character()
  )
}



# Get the n_closest airports to station no. 'i' in data frame 'data'
# Output: data frame with n_closest
get_closest <- function(i, data, airports, n_closest = 3){
  df <- as.data.frame(data)[i,]
  dx <- df$x - airports$x
  dy <- df$y - airports$y
  dist <- sqrt((dx^2) + (dy^2))/1000
  i_order <- order(dist)
  i_min <- i_order[1:n_closest]
  tibble(i = 1:n_closest,
         STATION_CODE = df$STATION_CODE, 
         x = df$x,
         y = df$y,
         Distance = dist[i_min],
         Airport = dat_airports[i_min,] %>% pull(LUFTHAVN),
         Airport_x = dat_airports[i_min,] %>% pull(x),
         Airport_y = dat_airports[i_min,] %>% pull(y)
  ) 
}

# TEST
if (FALSE)
  get_closest(1, dat_stations, dat_airports)


#
# As get_closest, but plots the n_closest airports   
#
plot_closest <- function(..., range = 200, print_plot = FALSE){
  
  
  df <- get_closest(...) %>%
    mutate(Label = paste(i, "-", stringr::str_extract(Airport, "([^[[:blank:]]]+)")))
  
  ran <- range*1000
  
  # plot UTM
  gg <- ggplot() +
    geom_polygon(data = data_norway, aes(X, Y, group = interaction(L1, L2)), fill = "grey70") +
    geom_point(data = df, aes(x = Airport_x, y = Airport_y)) +
    geom_text(data = df, aes(x = Airport_x + 5000, y = Airport_y, label = Label), hjust = 0) +
    geom_point(data = df[1,], aes(x=x, y=y), color = "red2", size = rel(3)) +
    coord_fixed(
      xlim = df[1,]$x + c(-ran/2, ran/2),
      ylim = df[1,]$y + c(-ran/2, ran/2)
    ) +
    theme_bw() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      panel.grid = element_blank()
    ) +
    labs(title = df[1,]$STATION_CODE)
  
  if (print_plot)
    print(gg)
  
  invisible(gg)
}

if (FALSE)
  plot_closest(1, dat_stations, dat_airports, print_plot = TRUE)

#
#
#

get_estimates_parallel <- function(param, data, 
                                   distance_data, tissue = "Lever", reference_year = 2019){
  
  df_regr <- data %>%
    filter(PARAM %in% param & TISSUE_NAME %in% tissue) 
  
  mod0 <- lm(log10(VALUE_WW) ~ MYEAR, data = df_regr)
  mod1 <- lm(log10(VALUE_WW) ~ MYEAR + STATION_CODE, data = df_regr)
  
  anova = anova(mod0, mod1)
  
  # summary(mod)
  
  df_estimates <- tibble(
    MYEAR = reference_year,
    STATION_CODE = unique(df_regr$STATION_CODE)
  )
  
  pred <- predict(mod1, df_estimates, se.fit = TRUE)
  
  df_estimates <- df_estimates %>%
    mutate(
      log_VALUE_WW = pred$fit,
      se = pred$se.fit
    ) %>%
    safe_left_join(distance_data, by = "STATION_CODE", 
                   na_matches = "never", 
                   check = "bCV")
  
  list(anova=anova, estimates = df_estimates)

  }



#
# Regression assuming each station has different slope (mod2 below)
# Estimate

get_estimates_max <- function(param, data, 
                                distance_data, tissue = "Lever", reference_year = 2009){
  
  df_regr <- data %>%
    filter(PARAM %in% param & TISSUE_NAME %in% tissue) 
  
  mod0 <- lm(log10(VALUE_WW) ~ MYEAR, data = df_regr)
  mod1 <- lm(log10(VALUE_WW) ~ MYEAR + STATION_CODE, data = df_regr)
  mod2 <- lm(log10(VALUE_WW) ~ MYEAR*STATION_CODE, data = df_regr)
  
  anova = anova(mod0, mod1, mod2)
  
  # summary(mod)
  
  # Get one line per year for each station
  df_estimates <- df_regr %>%
    distinct(STATION_CODE, MYEAR)

  # PRedict for all of these
  pred <- predict(mod2, df_estimates, se.fit = TRUE)
  
  df_estimates <- df_estimates %>%
    mutate(
      log_VALUE_WW = pred$fit,
      se = pred$se.fit
    ) %>%
    # Next three lines: for each station, pick the line with highest log_VALUE_WW
    group_by(STATION_CODE) %>%
    mutate(log_VALUE_WW_max = max(log_VALUE_WW)) %>%
    filter(log_VALUE_WW_max == log_VALUE_WW) %>%
    safe_left_join(distance_data, by = "STATION_CODE", 
                   na_matches = "never", 
                   check = "bCV")
  
  list(anova=anova, estimates = df_estimates)
  
}

# TEST
# result <- get_estimates_max(param = "PFOS", dat, df_distance)

plot_estimates <- function(param, ..., analysis = "parallel", distance_measure){
  if (distance_measure == "Distance1"){
    label = "distance to closest airport"
  } else if (distance_measure == "Distance2"){
    label = "distance to closest 'upstream' airport"
  } else if (distance_measure == "Dist_along_coast"){
    label = "distance along coast (S -> N)"
  } else if (distance_measure == "Dist_along_coast_Dist2"){
    label = "distance along coast (S -> N)"
  }
  
  if (analysis == "parallel"){
    result <- get_estimates_parallel(param = param, ...)
  } else if (analysis == "max") {
    result <- get_estimates_max(param = param, ...)
  }
  if (distance_measure == "Distance1"){
    mod <- lm(log_VALUE_WW ~ Distance1, data = result$estimates)
    est_effect <- summary(mod)$coef[2, "Estimate"] %>% round(3)
    p_effect <- summary(mod)$coef[2, "Pr(>|t|)"] %>% round(3)
    ggplot(result$estimates, aes(Distance1, log_VALUE_WW)) +
      geom_point() +
      geom_text(aes(x = Distance1 + 1, label = STATION_CODE), hjust = 0) +
      geom_smooth(method = "lm", formula = 'y ~ x') +
      annotate("text", x = -Inf, y = Inf, label = paste("Slope =", est_effect, ", P =", p_effect), 
               hjust = -0.3, vjust = 1.3, size = 4) +
      labs(title = paste(param, "-", label)) +
      theme_bw()
  } else if (distance_measure == "Distance2"){
    mod <- lm(log_VALUE_WW ~ Distance2, data = result$estimates)
    est_effect <- summary(mod)$coef[2, "Estimate"] %>% round(3)
    p_effect <- summary(mod)$coef[2, "Pr(>|t|)"] %>% round(3)
    ggplot(result$estimates, aes(Distance2, log_VALUE_WW)) +
      geom_point() +
      geom_text(aes(x = Distance2 + 1, label = STATION_CODE), hjust = 0) +
      geom_smooth(method = "lm", formula = 'y ~ x') +
      annotate("text", x = -Inf, y = Inf, label = paste("Slope =", est_effect, ", P =", p_effect), 
               hjust = -0.1, vjust = 1.2, size = 5) +
      labs(title = paste(param, "-", label)) +
      theme_bw()
  } else if (distance_measure == "Dist_along_coast"){
  mod <- lm(log_VALUE_WW ~ Dist_along_coast, data = result$estimates)
  est_effect <- summary(mod)$coef[2, "Estimate"] %>% round(3)
  p_effect <- summary(mod)$coef[2, "Pr(>|t|)"] %>% round(3)
  ggplot(result$estimates, aes(Dist_along_coast, log_VALUE_WW)) +
    geom_point() +
    geom_text(aes(x = Dist_along_coast + 40, label = STATION_CODE), hjust = 0) +
    geom_smooth(method = "lm", formula = 'y ~ x') +
    annotate("text", x = -Inf, y = Inf, label = paste("Slope =", est_effect, ", P =", p_effect), 
             hjust = -0.1, vjust = 1.2, size = 5) +
    labs(title = paste(param, "-", label)) +
    theme_bw()
  } else if (distance_measure == "Dist_along_coast_Dist2"){
  mod <- lm(log_VALUE_WW ~ Dist_along_coast + Distance2, data = result$estimates)
  est_effect1 <- summary(mod)$coef[2, "Estimate"] %>% round(3)
  p_effect1 <- summary(mod)$coef[2, "Pr(>|t|)"] %>% round(3)
  est_effect2 <- summary(mod)$coef[3, "Estimate"] %>% round(3)
  p_effect2 <- summary(mod)$coef[3, "Pr(>|t|)"] %>% round(3)
  par(mfrow = c(1,2), mar = c(5,5,2,1), oma = c(0,0,2,0))
  visreg::visreg(mod, "Dist_along_coast", points=list(cex=1.2, pch=1), 
                 main = paste("Slope =", est_effect1, ", P =", p_effect1))
  visreg::visreg(mod, "Distance2", points=list(cex=1.2, pch=1), 
                 main = paste("Slope =", est_effect2, ", P =", p_effect2))
  mtext(param, cex = 1.5, line = 1, outer = TRUE)
  mtext(param, cex = 1.5, line = 1, outer = TRUE)
  }
}


# TEST
if (FALSE){
  # debugonce(plot_estimates)
  plot_estimates(param = "PFBS", dat, df_distance, 
                 analysis = "max",
                 distance_measure = "Distance1")
  plot_estimates(param = "PFBS", dat, df_distance, 
                 analysis = "max",
                 distance_measure = "Dist_along_coast_Dist2")
}


# Like ramge, but expand the tange by a fraction down and up 
# Used in plots to avoid "clipping" of text
rangex <- function(x, down = 0, up = 0.1){
  ran <- range(x, na.rm = TRUE)
  diff <- diff(ran)
  ran + c(-down*diff, up*diff)
}
if (FALSE)
  rangex(10:20)