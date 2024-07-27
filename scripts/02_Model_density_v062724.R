# Tidyverse and ggpmisc
library(tidyverse)
suppressWarnings(library(ggpmisc)) 
# Data Import and Handling
library(vroom) 
library(readxl) 
# Statistical Modeling
library(mgcv)
library(gratia)
library(lmodel2) 
library(car)
# Progress Monitoring
library(progress)
# BBS data
library(bbsBayes2)
# Custom functions & data
source("scripts/00_Functions.R")
source("scripts/fit_model_with_fallback.R")
source("scripts/species_abundance_estimates.R")


# --- Repro ----

set.seed(1)


# ---- Setup ----

# Initial points
min_year <- 1966
detectability_factors <- read_xlsx("data/detectability-factors/detectability_factors.xlsx")

# Hash table for BCR names and numbers
bcrs <- data.frame(
  region = c("Atlantic Northern Forest","Northern Pacific Rainforest",
             "Appalachian Mountains", "Eastern Tallgrass Prairie",
             "Northern Rockies", "Shortgrass Prairie"),
  number = c(14, 5, 28, 22, 10, 18)) %>%
  mutate(bcr_code = paste0("BCR", number))

# Reproducibility
set.seed(1)


# ---- Select BCR region ----

bird_conservation_region <- "Atlantic Northern Forest"

# Get number ID
bcr_code <- bcrs %>% 
  filter(region == bird_conservation_region) %>% 
  pull(bcr_code)


# ---- Get species ----

# File name
BCR_file_suffix <- gsub(" ", "_", toupper(bird_conservation_region))
BCR_file <- paste("data/BCR-species/species_bcr_", BCR_file_suffix, ".csv", sep = "")

spp_initial <- vroom(BCR_file) %>% select(species)

sp_to_rm <- spp_initial %>% 
  mutate(species = gsub(" \\(all forms\\)", "", species)) %>%
  left_join(detectability_factors, by = "species") %>%
  filter(is.na(C_D)) %>%
  pull(species) %>%
  c(., "Eurasian Collared-Dove")

spp <- spp_initial %>% 
  filter(!(species %in% sp_to_rm)) %>% 
  pull(species)

n_spp <- length(spp)


# ---- Loop through species ----

pb <- progress_bar$new(
  format = "[:bar] :percent elapsed: :elapsed eta: :eta",
  total = n_spp, clear = FALSE)

options(warn = 1)
warnings_by_iteration <- list()

species_samples <- list()
species_estimates <- list()
species_parameter_variances <- list()
species_FirstLastDensities <- list()
for(i in 1:n_spp){
  
  # ---- Fetch the BBS data ----
  
  my_bcr <- bbs_strata %>%
    pluck("bcr") %>%
    filter(strata_name == bcr_code)
  
  s <- stratify(
    by = "bcr", 
    strata_custom = my_bcr,
    species = spp[i], 
    quiet = T)
  
  
  if(T){
    
    p <- prepare_data(s, quiet = T)
    
    
    # ---- Extract relevant BBS data -----
    
    # Observer data for route-year combos
    observer_df <- p$raw_data %>%
      select(year, route, observer, first_year, strat = strata_name)
    
    # Count data for route-year combos
    count_df <- s$birds_strata %>%
      select(species_total, route_data_id)
    
    # Route data for route-year combos
    route_df <- s$routes_strata %>%
      select(year, month, day, route,
             latitude, longitude, total_spp,
             start_temp, end_temp, temp_scale,
             start_time, end_time,
             start_wind, end_wind,
             start_sky, end_sky,
             route_data_id)
    
    
    # ---- Clean up relevant data ----
    
    df <- count_df %>%
      # Add route data to counts
      left_join(route_df, by = c("route_data_id")) %>%
      # Occasionally there is one row with duplicate data
      # So fix this by nesting by year and route
      # and keep only the first row
      filter(!is.na(year)) %>%
      filter(!is.na(route)) %>%
      nest(.by = c(year, route)) %>%
      mutate(data = map(data, ~ .x %>% slice(1))) %>%
      unnest(cols = c(data)) %>%
      # Add observer data to routes and counts
      left_join(x = ., y = observer_df, by = c("year","route")) %>%
      # Remove missing data
      na.omit() %>% # This removes instances where temp info is missing
      # Typically less than 30 observations
      # Convert bad time format to good, then to minutes of day
      mutate(start_time = sprintf("%04d", start_time),
             start_time = strptime(start_time, format = "%H%M"),
             start_time = hour(start_time) * 60 + minute(start_time)) %>%
      mutate(end_time = sprintf("%04d", end_time),
             end_time = strptime(end_time, format = "%H%M"),
             end_time = hour(end_time) * 60 + minute(end_time)) %>%
      # Convert mixed temp scale to all F, remove scale
      mutate(start_temp = if_else(temp_scale == "C",
                                  start_temp * 9/5 + 32, start_temp),
             temp_scale = if_else(temp_scale == "C",
                                  "F", temp_scale)) %>%
      select(-temp_scale) %>%
      # Date to yday
      mutate(date = make_date(year, month, day),
             yday = yday(date)) %>%
      mutate(yday = scale(yday)[,1]) %>%
      # Random effects to factors
      mutate(observer = as.factor(observer),
             route = as.factor(route),
             strat = as.factor(strat)) %>%
      # Arranging how I like the columns ordered
      select(observer, first_year, route, strat,
             year, month, day, yday,
             latitude, longitude, total_spp,
             start_temp, end_temp,
             start_time, end_time,
             start_wind, end_wind,
             start_sky, end_sky,
             species_total) %>%
      # Making a random effect for year
      mutate(year_factor = as.factor(as.character(year))) %>%
      mutate(year_strat_factor = as.factor(paste(year, strat, sep = "-"))) %>%
      mutate(duration = scale(end_time - start_time)[,1]) %>%
      mutate(year_full = year) %>%
      mutate(year = year - min_year) %>%
      mutate(od = row_number()) %>%
      mutate(od = as.factor(od))
    
    
    # ---- Modelling ----
    
    m <- fit_model_with_fallback(df)
    
    
    # ---- Calculate N0 and r ----
    
    # Get detectability-adjusted abundance
    det_facts_sp <- detectability_factors %>%
      filter(aou == unique(s$birds_strata$aou))
    
    C_P <- det_facts_sp$C_P
    C_D <- det_facts_sp$C_D
    DF <- ( (400/C_D)^2 ) * C_P
    
    area_original <- pi*(400^2) # in sq meters
    area_20ha <- 200000 # in sq meters
    
    # Number of samples
    n_samples <- 500
    
    # Get growth rate
    m_est_r <- coef(m)["year"]; names(m_est_r) <- "r"
    m_est_r_se <- summary(m, re.test = F)$se[2]; names(m_est_r_se) <- "se(r)"
    m_est_r_samples_link <- rnorm(n = n_samples, mean = m_est_r, sd = m_est_r_se)
    m_est_r_samples <- exp(m_est_r_samples_link)-1
    
    # Get count intercept
    m_est_y0 <- coef(m)["(Intercept)"]; names(m_est_y0) <- "y0"
    m_est_y0_se <- summary(m, re.test = F)$se[1]; names(m_est_y0_se) <- "se(y0)"
    N <- rnorm(n = n_samples, mean = m_est_y0, sd = m_est_y0_se)
    m_est_N0_samples <- log(area_20ha * (((exp(N)/50) * DF)/area_original))
    
    result <- data.frame(
      r = m_est_r_samples,
      N0 = m_est_N0_samples,
      species = spp[i]) %>%
      mutate(sample = row_number())
    
    # Save samples
    species_samples[[i]] <- result
    
    # Save point estimates
    species_estimates[[i]] <- result %>%
      group_by(species) %>%
      summarise(r_mean = mean(r),
                r_lwr = quantile(r, 0.025),
                r_upr = quantile(r, 0.975),
                N0_mean = mean(N0),
                N0_lwr = quantile(N0, 0.025),
                N0_upr = quantile(N0, 0.975)) %>%
      ungroup()
    
    
    # Variances
    variance_initial_abundance <- var(m_est_N0_samples)
    variance_growth_rate <- var(m_est_r_samples)
    
    species_parameter_variances[[i]] <- data.frame(
      var_initial_abundance = variance_initial_abundance,
      var_growth_rate = variance_growth_rate,
      species = spp[i])
    
    
    # ---- Calculate initial and final density ----
    
    species_FirstLastDensities[[i]] <- calculate_first_last_density20km(
                                          m = m, species = spp[i],
                                          DF = DF, 
                                          area_original = area_original, 
                                          area_20ha = area_20ha, 
                                          min_year = min_year)
    
  }
  
  # Capture and store warnings generated in this iteration
  iter_warnings <- warnings()
  if (length(iter_warnings) > 0) {
    # Store warnings in the list with the iteration number as the identifier
    warnings_by_iteration[[as.character(i)]] <- iter_warnings
  }
  
  # Reset warning state for the next iteration
  options(warn = 1)
  
  pb$tick()
  
}

# Save results
species_samples_df <- do.call(rbind, species_samples)
species_estimates_df <- do.call(rbind, species_estimates)
species_variances_df <- do.call(rbind, species_parameter_variances)
species_FirstLast_df <- do.call(rbind, species_FirstLastDensities)


# ---- Get resampled type 2 regression ----

# Progress bar
pb <- progress_bar$new(
  format = "[:bar] :percent elapsed: :elapsed eta: :eta",
  total = 100, clear = FALSE)

# Get samples from each regression
regr_samples <- list()
for(i in 1:100){
  
  # Select data
  df_i <- species_samples_df %>%
    group_by(species) %>%
    slice(i) %>%
    ungroup()
  
  # Fit model
  mod <- lmodel2(
    formula = r ~ N0,
    data = df_i,
    range.x = "interval", 
    range.y = "interval",
    nperm = 99)
  
  # Get estimated means
  means <- mod$regression.results %>%
    filter(Method == "RMA")
  mean_intercept <- means$Intercept
  mean_slope <- means$Slope
  
  # Get estimated intervals
  intervals <- mod$confidence.intervals %>%
    filter(Method == "RMA")
  upper_intercept <- intervals["97.5%-Intercept"]
  upper_slope <- intervals["97.5%-Slope"]
  
  # Get standard errors
  se_intercept <- as.numeric((upper_intercept - mean_intercept) / 1.96)
  se_slope <- as.numeric((upper_slope - mean_slope) / 1.96)
  
  # Type2Regr Intercept samples
  regr_intercept_samples <- rnorm(
    n = 100, 
    mean = mean_intercept, 
    sd = se_intercept)
  
  # Type2Regr Slope samples
  regr_slope_samples <- rnorm(
    n = 100, 
    mean = mean_slope, 
    sd = se_slope)
  
  # Type2Regr output
  results <- data.frame(
    intercept = regr_intercept_samples,
    slope = regr_slope_samples) %>%
    mutate(sample = i)
  
  # Save results
  regr_samples[[i]] <- results
  
  pb$tick()
  
}

regr_samples_df <- do.call(rbind, regr_samples)

hist(regr_samples_df$intercept, breaks = 50)
hist(regr_samples_df$slope, breaks = 50)


# ---- Calculating regressions ----

# Progress bar
pb <- progress_bar$new(
  format = "[:bar] :percent elapsed: :elapsed eta: :eta",
  total = nrow(regr_samples_df), clear = FALSE)

range_lwr <- min(species_estimates_df$N0_mean)
range_upr <- max(species_estimates_df$N0_mean)
x <- seq(range_lwr, range_upr, length.out = 100)

y_results <- list()
for(i in 1:nrow(regr_samples_df)){
  
  beta0 <- regr_samples_df$intercept[i]
  beta1 <- regr_samples_df$slope[i]
  y <- beta0 + beta1 * x
  
  y_df <- data.frame(y, x) %>%
    mutate(sample = i)
  
  y_results[[i]] <- y_df
  
  pb$tick()
}

y_results_df <- do.call(rbind, y_results)

final_relationship <- y_results_df %>%
  group_by(x) %>%
  summarise(mean = mean(y),
            lwr = quantile(y, 0.025),
            upr = quantile(y, 0.975)) %>%
  ungroup()

ggplot(final_relationship, aes(x = x, y = mean)) +
  
  geom_hline(yintercept = 0) +
  
  geom_errorbar(data = species_estimates_df,
                aes(x = N0_mean, y = r_mean, 
                    ymin = r_lwr, ymax = r_upr), 
                color = "gray70") +
  
  geom_errorbarh(data = species_estimates_df, 
                 aes(x = N0_mean, y = r_mean,  
                     xmin = N0_lwr, xmax = N0_upr), 
                 color = "gray70") +
  
  geom_point(data = species_estimates_df,
             aes(x = N0_mean, y = r_mean), 
             color = "gray70", size = 1) +
  
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "purple", alpha = 0.5) +
  geom_line(col = "purple", linewidth = 1.5) +
  
  scale_x_continuous(breaks = c(log(0.03), log(0.1), log(0.3), log(1), 
                                log(3),log(10),log(30),log(100)),
                     labels = c(0.03, 0.1, 0.3, 1,3,10,30,100)) +
  scale_y_continuous(labels = scales::percent_format()) +
  
  labs(x = "Initial abundance per 20 hectares ", 
       y = "Annual growth rate",
       subtitle = bird_conservation_region) +
  
  theme_light(16) +
  theme(aspect.ratio = 1,
        text = element_text(family = "Lato"),
        panel.grid.minor = element_blank())


# ---- Results ----

# Combine your data objects into a single list
combined_list <- list(
  FinalRelationship = final_relationship,
  SpeciesSamples = species_samples_df, 
  SpeciesEstimates = species_estimates_df,
  SpeciesVariances = species_variances_df,
  FirstLastDensities = species_FirstLast_df
)

# Save the combined list as an RDS file in your directory
BCR_output_file <- paste("data/BCR-results/results_bcr_", BCR_file_suffix, ".rds", sep = "")
saveRDS(combined_list, file = BCR_output_file)

