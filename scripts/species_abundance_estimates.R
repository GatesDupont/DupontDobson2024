require(tidyverse)
require(mgcv)
require(gratia)

# ---- Function ----

calculate_first_last_density20km <- function(m, species, DF, area_original, area_20ha, min_year){
  
  # ---- Recover average annual count estimates from model ----
  
  # Calculate average annual count estimates
  species_AvgAnnCt <- m %>%
    data_slice(
      year = unique(year),
      year_factor = levels(year_factor)) %>%
    filter(as.character(year + min_year) == year_factor) %>%
    fitted_values(object = m, data = ., scale = "link",
                  exclude = c("s(observer)","s(route)","s(duration)",
                              "s(yday)","te(start_temp,end_temp)", 
                              "start_wind", "first_year"))
  
  
  # ---- Simulate samples for average annual count estimates ----
  
  # Define the number of simulations
  n_sim <- 10000
  
  # Perform Monte Carlo simulation
  set.seed(1)  # for reproducibility
  simulations <- species_AvgAnnCt %>%
    filter(year %in% c(min(year), max(year))) %>%
    mutate(simulated_counts = map2(.fitted, .se, ~ rnorm(n_sim, mean = .x, sd = .y)))
  
  # Extract the first and last year
  first_year <- simulations %>%
    filter(year == min(year)) %>%
    pull(simulated_counts) %>%
    unlist()
  
  last_year <- simulations %>%
    filter(year == max(year)) %>%
    pull(simulated_counts) %>%
    unlist()
  
  # If a species is introduced or is extirpated, then count is zero
  if(min(simulations$year) > 0) first_year <- runif(n_sim, 0, 0)
  if(max(simulations$year) < 56) last_year <- runif(n_sim, 0, 0)
  
  
  # ---- Convert to density samples for average annual abundance ----
  
  calculate_20km_density <- function(x) area_20ha * (((exp(x)/50) * DF)/area_original)
  first_year_density20km <- data.frame(year = 0, density = calculate_20km_density(first_year)) 
  last_year_density20km <- data.frame(year = 56, density = calculate_20km_density(last_year))
  
  
  # ---- Ouptut ----
  
  # As dataframe object
  result <- rbind(
    first_year_density20km, 
    last_year_density20km) %>%
    mutate(species = species)
  
  return(result)
  
}