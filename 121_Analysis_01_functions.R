
#
# OVERVIEW ----
#
# LINEAR ANALYSIS - part 6
# analysis_one_par - analysis for one parameter/tissue
# - works by running 
#     analysis_one_par_one_repl 
#   several times (returns a one-line data frame each time)
# - which again does
#   1. selection of data by parameter and tissue  
#   2. add_random_data() for making random data under LOQ
#   3. add_adjusted_value_F() or add_adjusted_value_LF()
#      for adjusting for fat (F) and possibly also length (LF)
#   4. run_gls() for the actual analysis (by MYEAR + Dist_along_coast)  


# GAMM - part 7
# The main analysis (resulting in 'gamm_list') uses 'gamm_one_par_one_repl'
# 'gamm_one_par_one_repl':
#  - selects data, adds add_random_data
#  - then runs 'run_gamm' or 'run_gamm_int' (with/without interactions)
# 'run_gamm' runs 'run_gamm_model' which returns model_and_data, and then sends this to 
#    'run_gamm_extract_results' (the version used in the loops) or  
#    'run_gamm_extract_plotdata'
# 'run_gamm_model':  selects data, runs gamm, and returns both data and model
# 'run_gamm_extract_results': Get t-table from model as well as the visreg data for the model fit
# 'run_gamm_extract_plotdata': Get all the visreg data for the model 

#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# Utility functions ----
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

#
#  add_random_data  
#  - Makes new variable 'VALUE_WW_r' with random values for data under LOQ  
# - Make new variable 'VALUE_WW_r' with random values between
#  0.5*LOQ and LOQ for data undcer LOQ
#
# Returns data set
#
add_random_data <- function(data){
  
  # Set todata frame (not tibble)
  data <- as.data.frame(data)
  
  # Under LOQs
  sel <- !is.na(data$FLAG1)
  
  # Test that runif works for vactors of min, max:
  # runif(2, min = c(1,50), max = c(10,100))
  data$VALUE_WW_r <- data$VALUE_WW
  data$VALUE_WW_r[sel] <- runif(
    n = sum(sel), 
    min = data$LOQ[sel]/2,
    max = data$VALUE_WW[sel]
  )
  
  data
  
}

#
# Formats p-values
#
format_p <- function(p){
  case_when(
    p < 0.001 ~ "< 0.001",
    p < 0.01 ~ round(p, 3) %>% as.character(),
    p <= 1 ~ round(p, 2) %>% as.character()
  )
}
# format_p(0.0001)
# format_p(0.002334324)
# format_p(0.02334324)

sort_by_number <- function(x){
  number <- stringr::str_extract(x, "[0-9]+") %>% as.numeric()
  number[is.na(number)] <- 0
  order <- order(number)
  x[order]
}

if (FALSE){
  # Test
  sort_by_number(c("CB118", "CB28"))
}


move_to_back <- function(x, pattern){
  sel <- grepl(pattern, x)
  c(x[!sel], x[sel])
}

if (FALSE){
  # Test
  move_to_back(LETTERS, "B")
}

move_to_front <- function(x, pattern){
  sel <- grepl(pattern, x)
  c(x[sel], x[!sel])
}
if (FALSE){
  # Test
  move_to_front(LETTERS, "B")
}


#
# Define `df_parameter_groups` 
#
# data_check is your raw data (you want to check if all of them hava a parameter group)
# 

get_df_parameter_groups <- function(filename = "Input_data/47_df_par.csv", 
                                    data_check = dat, check = FALSE){
  
  group_order <- c("Metals and metalloids", "Chlorobiphenyls", 
                   "Polycyclic aromatic hydrocarbons (PAHs)",
                   "Organobromines", "Organochlorines (general)", "Organofluorines",
                   "Phosphorus flame retardant (PFR)",
                   "Phenols/chlorophenols", "Bisphenols",
                   "Chlorinated paraffins", "Dichloro-diphenyl-trichloroethane (DDTs)", 
                   "Hexachlorocyclohexanes", 
                   "Biological effects: molecular/biochemical/cellular/assays", 
                   "Organic esters", "Isotopes", "Cyclodienes", "Dioxins", "Biomarkers", 
                   "Phthalates", "Organo-metallic compounds", "Major inorganic constituents", 
                   "Triazines", "Siloxanes", "Chlorinated flame retardants",
                   "Others")
  
  df_parameter_groups <- read.csv2(filename, stringsAsFactors = FALSE)
  
  #
  # Extra metadata
  #
  df <- textConnection("
Parameter.Code;Substance.Group
MCCP inkl. LOQ;Chlorinated paraffins
MCCP eksl. LOQ;Chlorinated paraffins
SCCP eksl. LOQ;Chlorinated paraffins
KPAH;Polycyclic aromatic hydrocarbons (PAHs)
PAH16;Polycyclic aromatic hydrocarbons (PAHs)
Sum HBCD;Organobromines
% C;Others
% N;Others
BDE156;Organobromines
BDE17;Organobromines
BDE184;Organobromines
BDE191;Organobromines
BDE197;Organobromines
BDE206;Organobromines
BDE207;Organobromines
BPA;Bisphenols
C/N;Others
Delta15N;Others
DRYWT%;Others
Fett;Others
Oktaklorstyren (OCS);Organochlorines (general)
Aldrin;Organochlorines (general)
Endrin;Organochlorines (general)
Dieldrin;Organochlorines (general)
Mirex;Organochlorines (general)
alfa-Klordan (cis);Organochlorines (general)
gamma-Klordan (trans);Organochlorines (general)
Heptaklor;Organochlorines (general)
Heptaklor epoksid;Organochlorines (general)
;Organochlorines (general)
;Organochlorines (general)
Pentaklorbenzen (QCB);Chlorobiphenyls
") %>% 
    read.csv2(stringsAsFactors = FALSE) %>%
    mutate(Parameter.Name = Parameter.Code)
  # df
  
  df_parameter_groups <- bind_rows(df_parameter_groups, df)
  
  sel <- df_parameter_groups$Substance.Group == ""
  df_parameter_groups$Substance.Group[sel] <- "Others"
  
  if (check){
    # Check that 'Parameter.Code' covers all PARAM in the data
    # (not completey but good enough; check later using e.g. 'par  <- series_param_tissue$PARAM')
    par  <- unique(data_check$PARAM)
    par[!par %in% df_parameter_groups$Parameter.Code]
    
    # Check that 'group_order' covers all groups
    groups <- df_parameter_groups$Substance.Group %>% unique()
    groups[!groups %in% group_order] %>% dput()
  }
  
  # Order groups  
  df_parameter_groups <- df_parameter_groups  %>% # filter(grepl("PAH", PARAM))
    mutate(Substance.Group = factor(Substance.Group, levels = group_order))
  
  df_parameter_groups
  
}



#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# Adding adjusted value (adjusted for length + fat or fat only) ----    
#
# add_adjusted_value_LF
# add_adjusted_value_F
#
# Functions are used for 'data', where 'data' is the raw data filtered so it contains only one parameter/tissue
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o


add_adjusted_value_LF <- function(data, 
                                  LNMEA_fixed = 500, FAT_PERC_fixed = 40,
                                  response_variable = "log_CONC",
                                  plot = FALSE){
  
  # Rename given 'response_variable' to 'Response'
  sel <- names(data) == response_variable
  if (sum(sel) != 1)
    stop(paste("variable", response_variable, "not found!"))
  names(data)[sel] <- "Response" 
  
  # Select complete records only
  sel <- complete.cases(data %>% select(Response, LNMEA, FAT_PERC, STATION_CODE, MYEAR))
  data <- data[sel,]
  
  # Number of stations*years
  no_stationyears <- paste(data$STATION_CODE, data$MYEAR, sep = "_") %>% 
    unique() %>% length()
  
  if (no_stationyears > 1){   # Check if several years for this station; if not, no need to adjust for year
    data$fSTATION_CODE <- factor(data$STATION_CODE)
    data$fMYEAR <- factor(data$MYEAR)
    mod <- lm(Response  ~ LNMEA + FAT_PERC + fSTATION_CODE*fMYEAR, data = data)
  } else {
    mod <- lm(Response  ~ LNMEA + FAT_PERC, data = data)
  }
  
  if (plot){
    par(mfrow = c(1,4), mar = c(4,5,2,1))
    visreg::visreg(mod)
  }
  
  
  # summary(mod)
  
  pred_data <- data
  n <- nrow(data)
  pred_data$LNMEA <- rep(LNMEA_fixed, n)
  pred_data$FAT_PERC <- rep(FAT_PERC_fixed, n)
  pred <- predict(mod, pred_data)
  data$Response_adj <- pred + residuals(mod)
  
  
  # Change names back
  sel <- names(data) == "Response"
  names(data)[sel] <- response_variable
  sel <- names(data) == "Response_adj"
  names(data)[sel] <- paste0(response_variable, "_adj")
  
  # remove fMYEAR, if needed
  if (no_stationyears > 1){
    result <- as.data.frame(data %>% select(-fMYEAR,-fSTATION_CODE))
  } else {
    result <- as.data.frame(data)
  }
  invisible(result)
  
}

if (FALSE){
  
  param <- "NI"
  tissue <- "Lever"

  dat_param <-  dat2 %>%
    filter(PARAM %in% param & TISSUE_NAME %in% tissue) %>% 
    add_random_data() %>%
    mutate(log_CONC = log10(VALUE_WW)) %>%
    select(PARAM, STATION_CODE, TISSUE_NAME, MSTAT, log_CONC, MYEAR, Dist_along_coast, LNMEA, FAT_PERC, TL, TL_mean) # %>%
  
  dim(dat_param)
  dat_param2 <- add_adjusted_value_LF(dat_param, LNMEA_fixed = 500, FAT_PERC_fixed = 40, plot = FALSE)
  cat("\n")
  dim(dat_param2)

  }

#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
# Adjusting for length only    
# Used for muscle
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

add_adjusted_value_L <- function(data, 
                                  LNMEA_fixed = 500, 
                                  response_variable = "log_CONC",
                                  plot = FALSE){
  
  # Rename given 'response_variable' to 'Response'
  sel <- names(data) == response_variable
  if (sum(sel) != 1)
    stop(paste("variable", response_variable, "not found!"))
  names(data)[sel] <- "Response" 
  
  # Select complete records only
  sel <- complete.cases(data %>% select(Response, LNMEA,  STATION_CODE, MYEAR))
  data <- data[sel,]
  
  # Number of stations*years
  no_stationyears <- paste(data$STATION_CODE, data$MYEAR, sep = "_") %>% 
    unique() %>% length()
  
  if (no_stationyears > 1){   # Check if several years for this station; if not, no need to adjust for year
    data$fSTATION_CODE <- factor(data$STATION_CODE)
    data$fMYEAR <- factor(data$MYEAR)
    mod <- lm(Response  ~ LNMEA + fSTATION_CODE*fMYEAR, data = data)
  } else {
    mod <- lm(Response  ~ LNMEA, data = data)
  }
  
  if (plot){
    par(mfrow = c(1,4), mar = c(4,5,2,1))
    visreg::visreg(mod)
  }
  
  
  # summary(mod)
  
  pred_data <- data
  n <- nrow(data)
  pred_data$LNMEA <- rep(LNMEA_fixed, n)
  # pred_data$FAT_PERC <- rep(FAT_PERC_fixed, n)
  pred <- predict(mod, pred_data)
  data$Response_adj <- pred + residuals(mod)
  
  
  # Change names back
  sel <- names(data) == "Response"
  names(data)[sel] <- response_variable
  sel <- names(data) == "Response_adj"
  names(data)[sel] <- paste0(response_variable, "_adj")
  
  # remove fMYEAR, if needed
  if (no_stationyears > 1){
    result <- as.data.frame(data %>% select(-fMYEAR,-fSTATION_CODE))
  } else {
    result <- as.data.frame(data)
  }
  invisible(result)
  
}

if (FALSE){
  
  param <- "NI"
  tissue <- "Lever"
  
  dat_param <-  dat2 %>%
    filter(PARAM %in% param & TISSUE_NAME %in% tissue) %>% 
    add_random_data() %>%
    mutate(log_CONC = log10(VALUE_WW)) %>%
    select(PARAM, STATION_CODE, TISSUE_NAME, MSTAT, log_CONC, MYEAR, Dist_along_coast, LNMEA, FAT_PERC, TL, TL_mean) # %>%
  
  dim(dat_param)
  dat_param2 <- add_adjusted_value_LF(dat_param, LNMEA_fixed = 500, FAT_PERC_fixed = 40, plot = FALSE)
  cat("\n")
  dim(dat_param2)
  
}


#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
# Adjusting for  for fat only    
# Used fr blue mussel, where we lack length for a lot of years  
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

add_adjusted_value_F <- function(data, 
                                 LNMEA_fixed = 500, FAT_PERC_fixed = 40,
                                 response_variable = "log_CONC",
                                 plot = FALSE){
  
  # Rename given 'response_variable' to 'Response'
  sel <- names(data) == response_variable
  if (sum(sel) != 1)
    stop(paste("variable", response_variable, "not found!"))
  names(data)[sel] <- "Response" 
  
  # Select complete records only
  sel <- complete.cases(data %>% select(Response, FAT_PERC, STATION_CODE, MYEAR))
  data <- data[sel,]
  
  # Number of stations*years
  no_stationyears <- paste(data$STATION_CODE, data$MYEAR, sep = "_") %>% 
    unique() %>% length()
  
  if (no_stationyears > 1){   # Check if several years for this station; if not, no need to adjust for year
    data$fSTATION_CODE <- factor(data$STATION_CODE)
    data$fMYEAR <- factor(data$MYEAR)
    mod <- lm(Response  ~ FAT_PERC + fSTATION_CODE*fMYEAR, data = data)
  } else {
    mod <- lm(Response  ~ FAT_PERC, data = data)
  }
  
  if (plot){
    par(mfrow = c(1,4), mar = c(4,5,2,1))
    visreg::visreg(mod)
  }
  
  
  # summary(mod)
  
  pred_data <- data
  n <- nrow(data)
  pred_data$FAT_PERC <- rep(FAT_PERC_fixed, n)
  pred <- predict(mod, pred_data)
  data$Response_adj <- pred + residuals(mod)
  
  
  # Change names back
  sel <- names(data) == "Response"
  names(data)[sel] <- response_variable
  sel <- names(data) == "Response_adj"
  names(data)[sel] <- paste0(response_variable, "_adj")
  
  # remove fMYEAR, if needed
  if (no_stationyears > 1){
    result <- as.data.frame(data %>% select(-fMYEAR,-fSTATION_CODE))
  } else {
    result <- as.data.frame(data)
  }
  invisible(result)
  
}

if (FALSE){
  #
  # Note: add_adjusted_value_LF() may remove some records
  #
  dim(dat_param)
  dat_param2 <- add_adjusted_value_F(dat_param, plot = FALSE)
  cat("\n")
  dim(dat_param2)
  
}


#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# Functions for linear analysis using gls  ----
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o


#
# Function 'run_gls'    
#
# - Runs gls for 5 models (combinations of MYEAR + Dist_along_coast)
# - Returns a one-line data frame
#
# Used in analysis_one_par_one_repl
#

run_gls <- function(data, response_variable = "log_CONC_adj"){
  
  # Rename given 'response_variable' to 'Response'
  sel <- names(data) == response_variable
  if (sum(sel) != 1)
    stop(paste("variable", response_variable, "not found!"))
  names(data)[sel] <- "Response" 
  
  
  data <- data %>%
    select(STATION_CODE, Response, MYEAR, Dist_along_coast) %>%
    filter(complete.cases(.))  # %>%
  # group_by(STATION_CODE, MYEAR, Dist_along_coast) %>%
  # summarize(log_CONC_adj = mean(log_CONC_adj))
  
  mod10 <- gls(Response ~ 1, data = data,
               correlation = corAR1(form = ~1 | STATION_CODE))
  mod1a <- gls(Response ~ MYEAR, data = data,
               correlation = corAR1(form = ~1 | STATION_CODE))
  mod1b <- gls(Response ~ Dist_along_coast, data = data,
               correlation = corAR1(form = ~1 | STATION_CODE))
  mod1c <- gls(Response ~ MYEAR + Dist_along_coast, data = data,
               correlation = corAR1(form = ~1 | STATION_CODE))
  mod1d <- gls(Response ~ MYEAR*Dist_along_coast, data = data,
               correlation = corAR1(form = ~1 | STATION_CODE))
  
  aic <- AIC(mod10, mod1a, mod1b, mod1c, mod1d)
  aic$dAIC <- aic$AIC - min(aic$AIC)
  
  best_model <- which.min(aic$AIC)
  
  # Resuts from model with MYEAR + Dist_along_coast, without interaction
  ttable_values <- summary(mod1c)$tTable[2:3, c(1,2,4)] %>% as.numeric()
  result <- c(best_model, 
              aic$dAIC[2:5],
              ttable_values) %>% 
    matrix(nrow = 1) %>% 
    data.frame(stringsAsFactors = FALSE)
  names(result) <- c("Best_model",
                     "dAIC_yr", "dAIC_dist", "dAIC_yr_dist", "dAIC_yr_x_dist",
                     "Year_est", "Position_est",
                     "Year_se", "Position_se",
                     "Year_p", "Position_p")
  result
}

# Used in section 6
model_string <- data.frame(
  modelnumber = 1:5,
  model = c("~", 
            "~ MYEAR",
            "~ Dist_along_coast",
            "~ MYEAR + Dist_along_coast",
            "~ MYEAR*Dist_along_coast")
)

if (FALSE){
  debugonce(run_gls)
  run_gls(dat_param2)
}

#
# analysis_one_par_one_repl
#
# Linear analysis by year and dist. along coast
# Function does analysis for *one replicate*
#
analysis_one_par_one_repl <- function(param, tissue, data){
  
  # Select data, addrandom values for data below LOC, and calculate the log concentration
  data_selected <- data %>%
    filter(PARAM %in% param & TISSUE_NAME %in% tissue & !grepl("F", STATION_CODE)) %>% 
    add_random_data() %>%
    mutate(log_CONC = log10(VALUE_WW),
           log_CONC_r = log10(VALUE_WW_r)) %>%
    select(STATION_CODE, TISSUE_NAME, MSTAT, log_CONC, log_CONC_r, MYEAR, Dist_along_coast, LNMEA, FAT_PERC, TL, TL_mean)
  
  
  # Add adjusted value
  if (tissue %in% "Whole soft body"){
    data_selected <- add_adjusted_value_F(data_selected,
                                          FAT_PERC_fixed = 40, 
                                          plot = FALSE, 
                                          response_variable = "log_CONC_r")
  } else {
    data_selected <- add_adjusted_value_LF(data_selected,
                                           LNMEA_fixed = 500, FAT_PERC_fixed = 40, 
                                           plot = FALSE, 
                                           response_variable = "log_CONC_r")
  }
  
  # Run analysis
  result <- run_gls(data_selected,  response_variable = "log_CONC_r_adj") 
  
  # Output
  data.frame(PARAM = param, TISSUE_NAME = tissue, result,
             stringsAsFactors = FALSE)
  
}

if (FALSE){
  # Testing 
  analysis_one_par_one_repl("NI", "Lever", dat2)
  analysis_one_par_one_repl("NI", "Whole soft body", dat2)
  analysis_one_par_one_repl("PFUdA", "Lever", dat2)
  # analysis_one_par_one_repl("CD", "Whole soft body", dat2)
  # debugonce(analysis_one_par_one_repl)
  # debugonce(add_adjusted_value_LF)
  # analysis_one_par_one_repl("HG", "Muskel", dat2)
}


#
# analysis_one_par
#
# Function for n replicates  
# Replicates needed since 'add_random_data' will give different results each time  


analysis_one_par <- function(param, tissue, data, nrepl = 5){
  
  seq(1, nrepl) %>%
    set_names() %>%
    map_dfr(
      ~analysis_one_par_one_repl(param = param, 
                                 tissue = tissue, 
                                 data = data),
      .id = "Repl")
  
}

if (FALSE){
  analysis_one_par("NI", "Lever", dat2, nrepl = 2)
  analysis_one_par("PFUdA", "Lever", dat2, nrepl = 5)
}




#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# Functions for GAMM models  ----
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

#
# Function 1
#
run_gamm_model <- function(data, response_variable = "log_CONC_adj"){
  
  # Rename given 'response_variable' to 'Response'
  sel <- names(data) == response_variable
  if (sum(sel) != 1)
    stop(paste("variable", response_variable, "not found!"))
  names(data)[sel] <- "Response" 
  
  data <- data %>%
    select(PARAM, TISSUE_NAME, STATION_CODE, Response, MYEAR, Dist_along_coast) %>%
    filter(complete.cases(.))  # %>%
  # group_by(STATION_CODE, MYEAR, Dist_along_coast) %>%
  # summarize(log_CONC_adj = mean(log_CONC_adj))
  
  n_years <- length(unique(data$MYEAR))
  n_stations <- length(unique(data$STATION_CODE))

  # For series with few years (n_years < 3), we cannot fit a smooth for year, so
  #   we need an "if" here 
  if (n_stations >= 9){
    if (n_years == 1) {
      # No MYEAR in model
      mod <- gamm(Response ~ s(Dist_along_coast, k = 8), 
                  data = data,
                  correlation = corAR1(form = ~1 | STATION_CODE))
    } else if (n_years <= 3) {
      # MYEAR as linear effect in model
      mod <- gamm(Response ~ MYEAR + s(Dist_along_coast, k = 8), 
                  data = data,
                  correlation = corAR1(form = ~1 | STATION_CODE))
    } else {
      # MYEAR as non-linear effect in model
      mod <- gamm(Response ~ s(MYEAR, k = 4) + s(Dist_along_coast, k = 8), 
                  data = data,
                  correlation = corAR1(form = ~1 | STATION_CODE))
    }
  } else {
    if (n_years == 1) {
      # No MYEAR in model
      mod <- gamm(Response ~ s(Dist_along_coast, k = 5), 
                  data = data,
                  correlation = corAR1(form = ~1 | STATION_CODE))
    } else if (n_years <= 3) {
      # MYEAR as linear effect in model
      mod <- gamm(Response ~ MYEAR + s(Dist_along_coast, k = 5), 
                  data = data,
                  correlation = corAR1(form = ~1 | STATION_CODE))
    } else {
      # MYEAR as non-linear effect in model
      mod <- gamm(Response ~ s(MYEAR, k = 4) + s(Dist_along_coast, k = 5), 
                  data = data,
                  correlation = corAR1(form = ~1 | STATION_CODE))
    }
  }
  
  
  list(mod = mod, data = data)
  
}

#
# Function 2
#
run_gamm_extract_results <- function(model, data, repl = 1){
  
  # Get parameter and tissue (will be saved in resulting data frames)
  parameter <- data$PARAM %>% unique() %>% paste(collapse = ",")
  tissue_name <- data$TISSUE_NAME %>% unique() %>% paste(collapse = ",")
  
  
  #
  # 1. "ttable" object (t-table)
  #
  # Get t-table from model (data frame of 1 row)
  s_table <- summary(model$gam)$s.table
  
  # For series with few years, we fit only smooth for distance, so
  #   we need an "if" here (see 'run_gamm_model')  
  if (nrow(s_table) == 2){
    # Normal case - smooths are fit for both time and distance 
    ttable_values <- s_table[,c(1,3,4)] %>% as.numeric()
  } else if (nrow(s_table) == 1){
    # 1-3 years - smooths are NA, ttable_values[1], only
    ttable_values <- s_table[,c(1,3,4)] %>% as.numeric()
    # Add NA values for the lacking time values 
    ttable_values <- c(NA, ttable_values[1],
                       NA, ttable_values[2],
                       NA, ttable_values[3])
  }
  # Add parameter and tissue to t values
  ttable <- c(parameter, tissue_name,
              repl,
              ttable_values) %>%
    matrix(nrow = 1) %>% 
    data.frame(stringsAsFactors = FALSE)
  names(ttable) <- c(
    "PARAM", "TISSUE_NAME", "Repl",
    "Year_edf", "Position_esdf",
    "Year_F", "Position_F",
    "Year_p", "Position_p")
  
  #
  # 2. "plotvalues" object 
  #
  # Get plot values (data frame of 100 rows)
  model$gam$data <- data   # needed, otherwise visreg doesn't work
  plotvalues <- visreg::visreg(
    model$gam, "Dist_along_coast", plot = FALSE)$fit
  
  # Add parameter and tissue to plot values
  plotvalues <- plotvalues %>%
    mutate(PARAM = parameter,
           TISSUE_NAME = tissue_name,
           Repl = repl) %>%
    select(PARAM, TISSUE_NAME, Repl, everything())
  
  # Final result = list of two data frames
  list(ttable=ttable, plotvalues=plotvalues)
  
}

#
# A variant of function 2 - returns all plot values   
#
run_gamm_extract_plotvalues <- function(model, data, repl = 1){
  
  # Get parameter and tissue (will be saved in resulting data frames)
  parameter <- data$PARAM %>% unique() %>% paste(collapse = ",")
  tissue_name <- data$TISSUE_NAME %>% unique() %>% paste(collapse = ",")
  
  
  #
  # 2. "plotvalues" object 
  #
  # Get plot values (data frame of 100 rows)
  model$gam$data <- data   # needed, otherwise visreg doesn't work
  plotvalues_list <- visreg::visreg(
    model$gam, "Dist_along_coast", plot = FALSE)
  
  plotvalues_list
  
}

#
# Function 3 - main function which calls the two others  
#
run_gamm <- function(data, response_variable = "log_CONC_adj", repl = 1,
                     return = "ttable_and_fit"){
  model_and_data <- run_gamm_model(data = data, 
                                   response_variable = response_variable)
  if (return == "ttable_and_fit"){
    result <-  run_gamm_extract_results(model = model_and_data$mod, 
                                        data = model_and_data$data, 
                                        repl = repl)
  } else if (return == "plotvalues"){
    result <-  run_gamm_extract_plotvalues(model = model_and_data$mod, 
                                        data = model_and_data$data, 
                                        repl = repl)
  } else if (return == "model_and_data"){
    model_and_data$mod$gam$data <- model_and_data$data 
    result <- model_and_data
  }
  result
  
}

#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# ... for testing ----
#
# Testing 'run_gamm_model' and 'run_gamm_extract_results'
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o

if (FALSE){
  
  # Test function 1
  x <- run_gamm_model(dat_param2)
  str(x, 1)
  str(x$mod, 1)
  
  # "Standard plot
  par(mfrow = c(1,2), mar = c(4,5,2,1))
  plot(x$mod$gam, res = TRUE)
  
  # visreg plot
  x$mod$gam$data <- x$data   # needed, otherwise visreg doesn't work
  par(mfrow = c(1,2), mar = c(4,5,2,1))
  visreg::visreg(x$mod$gam, "MYEAR")
  visreg::visreg(x$mod$gam, "Dist_along_coast")
  
  # Test function 2
  debugonce(run_gamm_extract_results)
  y <- run_gamm_extract_results(x$mod, x$data, repl = 1)
  str(y, 1)
  y$ttable
  
}


#
# Function for one replicate      
# This is actually the main version used  
#
# Can return 3 different tyoes of output, depending on 'return'
#   'ttable_and_fit' is the defaukt and the ine used for all parameters  
#

# interaction = FALSE is used in this part (section 7-9)  
# interaction = TRUE is used from section 10 onwards  
gamm_one_par_one_repl <- function(param, tissue, data, repl = 1, interaction = FALSE,
                                  return = "ttable_and_fit"){
  
  # Select data, addrandom values for data below LOC, and calculate the log concentration
  data_selected <- data %>%
    filter(PARAM %in% param & TISSUE_NAME %in% tissue & !grepl("F", STATION_CODE)) %>% 
    add_random_data() %>%
    mutate(log_CONC = log10(VALUE_WW),
           log_CONC_r = log10(VALUE_WW_r)) %>%
    select(PARAM, STATION_CODE, TISSUE_NAME, MSTAT, log_CONC, log_CONC_r, 
           MYEAR, Dist_along_coast, LNMEA, FAT_PERC, TL, TL_mean)
  
  # Add adjusted value
  if (tissue %in% "Whole soft body"){
    data_selected <- add_adjusted_value_F(
      data_selected,
      FAT_PERC_fixed = 40, 
      plot = FALSE, 
      response_variable = "log_CONC_r")
  } else if (tissue %in% "Muskel"){
    data_selected <- add_adjusted_value_L(
      data_selected,
      LNMEA_fixed = 500, 
      plot = FALSE, 
      response_variable = "log_CONC_r")
  } else {
    data_selected <- add_adjusted_value_LF(
      data_selected,
      LNMEA_fixed = 500, FAT_PERC_fixed = 40, 
      plot = FALSE, 
      response_variable = "log_CONC_r")
  }
  
  # Run analysis
  if (!interaction){
    result <- run_gamm(data_selected,  response_variable = "log_CONC_r_adj", repl = repl,
                       return = return) 
  } else {
    # run_gamm_int (GAMMwith interaction) is defined further down! (part 10)
    result <- run_gamm_int(data_selected,  response_variable = "log_CONC_r_adj", repl = repl,
                           return = return) 
  }
  # Output
  result
  
}


# As gamm_one_par_one_repl(), but returns '
#  - uses run_gamm_model() instaed of run_gamm() in the end
gamm_one_par_one_repl_model <- function(param, tissue, data, repl = 1, interaction = FALSE){
  
  # Select data, addrandom values for data below LOC, and calculate the log concentration
  data_selected <- data %>%
    filter(PARAM %in% param & TISSUE_NAME %in% tissue) %>% 
    add_random_data() %>%
    mutate(log_CONC = log10(VALUE_WW),
           log_CONC_r = log10(VALUE_WW_r)) %>%
    select(PARAM, STATION_CODE, TISSUE_NAME, MSTAT, log_CONC, log_CONC_r, 
           MYEAR, Dist_along_coast, LNMEA, FAT_PERC, TL, TL_mean)
  
  
  # Add adjusted value
  if (tissue %in% "Whole soft body"){
    data_selected <- add_adjusted_value_F(
      data_selected,
      FAT_PERC_fixed = 40, 
      plot = FALSE, 
      response_variable = "log_CONC_r")
  } else {
    data_selected <- add_adjusted_value_LF(
      data_selected,
      LNMEA_fixed = 500, FAT_PERC_fixed = 40, 
      plot = FALSE, 
      response_variable = "log_CONC_r")
  }
  
  # Run analysis
  if (!interaction){
    result <- run_gamm(data_selected,  response_variable = "log_CONC_r_adj", repl = repl) 
    model_and_data <- run_gamm_model(data = data_selected, 
                                     response_variable = "log_CONC_r_adj")
  } else {
    # run_gamm_int (GAMMwith interaction) is defined further down! (part 10)
    model_and_data <- run_gamm_model_int(data = data_selected, 
                                     response_variable = "log_CONC_r_adj")
  }
  # Output
  model_and_data
  
}


#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# Functions for GAMM with interaction  ----
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o


### Functions for GAMM with interactions  
# Functions named as functions above, with '_int' added, i.e.,  
# - run_gamm_model_int  
# - run_gamm_extract_results_int  
# - run_gamm_int

#
# Function 1
#
run_gamm_model_int <- function(data, response_variable = "log_CONC_adj"){
  
  # Rename given 'response_variable' to 'Response'
  sel <- names(data) == response_variable
  if (sum(sel) != 1)
    stop(paste("variable", response_variable, "not found!"))
  names(data)[sel] <- "Response" 
  
  data <- data %>%
    select(PARAM, TISSUE_NAME, STATION_CODE, Response, MYEAR, Dist_along_coast) %>%
    filter(complete.cases(.))  # %>%
  # group_by(STATION_CODE, MYEAR, Dist_along_coast) %>%
  # summarize(log_CONC_adj = mean(log_CONC_adj))
  
  n_years <- length(unique(data$MYEAR))
  
  # For series with few years (n_years < 3), we cannot fit a smooth for year, so
  #   we need an "if" here   
  if (n_years < 5) {
    mod1 <- NULL
    mod2 <- NULL
  } else {
    # MYEAR as non-linear effect in model
    mod1 <- gamm(
      Response ~ te(MYEAR, k = 4) + te(Dist_along_coast, k = 8), 
      data = data,
      correlation = corAR1(form = ~1 | STATION_CODE))
    mod2 <- gamm(
      Response ~ te(MYEAR, k = 4) + te(Dist_along_coast, k = 8) + te(Dist_along_coast, MYEAR, k = 8), 
      data = data,
      correlation = corAR1(form = ~1 | STATION_CODE))
  }
  
  list(mod1 = mod1, mod2 = mod2, data = data)
  
}

#
# Function 2
#
run_gamm_extract_results_int <- function(model1, model2, data, repl = 1){
  
  # Get parameter and tissue (will be saved in resulting data frames)
  parameter <- data$PARAM %>% unique() %>% paste(collapse = ",")
  tissue_name <- data$TISSUE_NAME %>% unique() %>% paste(collapse = ",")
  
  if (!is.null(model2)){
    #
    # 1. "anova" object
    #
    # Get t-table from model (data frame of 1 row)
    anova_table <- anova(model1$lme, model2$lme)
    
    anova <- data.frame(
      PARAM = parameter, TISSUE_NAME = tissue_name, Repl = repl,
      ddf = diff(anova_table$df),
      dAIC = diff(anova_table$AIC),
      p_value = anova_table$'p-value'[2],
      stringsAsFactors = FALSE
    )
    
    #
    # 2. "plotvalues" object 
    #
    
    year_range <- range(data$MYEAR)
    
    # 2 years if the data spans =< 10 years, 3 years if 11-20 years, etc. 
    fit_years <- seq(
      year_range[1], 
      year_range[2], 
      length = ceiling(diff(year_range)/10) + 1) %>%  floor()
    
    model2$gam$data <- data
    
    plotvalues <- fit_years %>% 
      map_dfr(
        ~visreg::visreg(model2$gam, "Dist_along_coast", cond = list(MYEAR = .), plot = FALSE)$fit
      ) %>%
      mutate(MYEAR = factor(MYEAR))
    
    plotvalues <- plotvalues %>%
      mutate(PARAM = parameter,
             TISSUE_NAME = tissue_name,
             Repl = repl) %>%
      select(PARAM, TISSUE_NAME, Repl, everything())
    
    
  } else {
    
    anova <- data.frame(
      PARAM = parameter, TISSUE_NAME = tissue_name, Repl = repl,
      ddf = NA,
      dAIC = NA,
      p_value = NA,
      stringsAsFactors = FALSE)
    
  }
  
  
  # Final result = list of two data frames
  list(anova=anova, plotvalues=plotvalues)
  
}


#
# Function 3 - main function which calls the two others  
#
run_gamm_int <- function(data, response_variable = "log_CONC_adj", repl = 1){
  model_and_data <- run_gamm_model_int(data = data, 
                                       response_variable = response_variable)
  run_gamm_extract_results_int(model1 = model_and_data$mod1, 
                               model2 = model_and_data$mod2, 
                               data = model_and_data$data, 
                               repl = repl)
  
}

#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# ... for testing ----
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o

if (FALSE){
  
  #o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
  # Test function 1
  #o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
  
  # debugonce(run_gamm_model_int)
  x <- run_gamm_model_int(dat_param2)
  str(x, 1)
  str(x$mod2, 1)
  anova(x$mod1$lme, x$mod2$lme)
  
  #
  # Get and plot fitted model, test 1
  #
  x$mod2$gam$data <- x$data
  x$mod2$gam$k <- 8
  df1 <- visreg::visreg(x$mod2$gam, "Dist_along_coast", cond = list(MYEAR = 1996), plot = FALSE)$fit
  df2 <- visreg::visreg(x$mod2$gam, "Dist_along_coast", cond = list(MYEAR = 2019), plot = FALSE)$fit
  ggplot(bind_rows(df1,df2), aes(Dist_along_coast, visregFit, group = MYEAR, color = MYEAR)) + geom_line()
  xx <- df2$visregFit - df1$visregFit
  plot(xx)
  
  #
  # Get and plot fitted model, test 2 (more general)
  #
  year_range <- range(x$data$MYEAR)
  
  # 2 years if the data spans =< 10 years, 3 years if 11-20 years, etc. 
  fit_years <- seq(
    year_range[1], 
    year_range[2], 
    length = ceiling(diff(year_range)/10)) %>%  floor()
  fit_years
  
  fits <- fit_years %>% 
    map_dfr(
      ~visreg::visreg(x$mod2$gam, "Dist_along_coast", cond = list(MYEAR = .), plot = FALSE)$fit
    ) %>%
    mutate(MYEAR = factor(MYEAR))
  ggplot(fits, aes(Dist_along_coast, visregFit, group = MYEAR, color = MYEAR)) + geom_line()
  
  
  #o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
  # Test function 2
  #o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
  
  # debugonce(run_gamm_extract_results_int)
  y <- run_gamm_extract_results_int(x$mod1, x$mod2, x$data, repl = 1)
  str(y, 1)
  y$anova
  ggplot(y$plotvalues, aes(Dist_along_coast, visregFit, group = MYEAR, color = MYEAR)) + geom_line()
  
  #o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
  # Test function 3
  #o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
  
  # debugonce(run_gamm_model_int)
  # debugonce(run_gamm_extract_results_int)
  z <- run_gamm_int(dat_param2, repl = 1)
  str(z, 1)
  ggplot(z$plotvalues, aes(Dist_along_coast, visregFit, group = MYEAR, color = MYEAR)) + geom_line()
  
  
}






