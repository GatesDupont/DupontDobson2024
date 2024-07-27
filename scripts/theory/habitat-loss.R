library(tidyverse)
library(lmodel2)
library(deSolve)
library(ggthemes)
library(ggpmisc)

set.seed(0)


# ---- MTE functions ----

r <- function(r_0, M){
  
  r_0 * (M^(-0.25))
  
}

K <- function(R, M){
  
  R * (M^(-0.75))
  
} 


# ---- Model -----

# Define the model for logistic growth of each species
lv_multi <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    dN <- numeric(n_species)
    
    for (i in 1:n_species) {
      
      N <- state[i]
      r <- parameters[["r"]][i]
      
      K_initial <- parameters[["K"]][i]
      K_decay_rate <- parameters[["K_decay_rate"]][i]
      K <- K_initial * exp(K_decay_rate * t)
      
      dN[i] <- r * N * ((K - N) / K)
    }
    
    return(list(dN))
  })
}


# ---- Population parameters ----

# Bird species
n_species <- 10000

# Mass in kilograms
Mass <- rlnorm(n = n_species, meanlog = log(0.043), sdlog = log(5.2))

# Intrinsic growth rate
r_values <- r(r_0 = 0.01, M = Mass)

# Carrying capacity
K_values <- K(R = 0.5, M = Mass)

# Initial abundance
N0_values <- runif(n_species, min = 0.01, max = 1) * K_values

# Decay rate for carrying capacity (2% per year)
K_decay_rate <- rnorm(n_species, mean = -0.04, sd = 0.02)
# hist(K_decay_rate)

# Hash table
spp_params_hash <- data.frame(
  Species_numeric = 1:n_species,
  initNfromK = N0_values/K_values
)

# ----- Organize parts of running model ----

# Create named vectors for initial state and parameters
initial_state <- c(N0_values)
parameters <- list(
  r = r_values, 
  K = K_values,
  K_decay_rate = K_decay_rate)

# Time points to solve the ODE
observed_times <- seq(0, 56, by = 0.1)


# ---- Run model ----

# Run the model and organize results
model_output <- ode(y = initial_state, times = observed_times, func = lv_multi, parms = parameters) %>%
  as.data.frame() %>%
  pivot_longer(-time, names_to = "Species", values_to = "Population") %>%
  mutate(Species = as.factor(Species))

if(F) {
  
  ggplot(model_output, aes(x = time, y = Population, color = Species)) +
    # facet_wrap(~Species, scales = "free") +
    geom_line() +
    scale_y_log10() +
    theme(legend.position = "none")
  
}


# ---- Analyze results for common vs rare ----

# Function to get initial abundance and growth rate for each species
get_ests <- function(y,x){
  
  m <- glm(y ~ x, family = gaussian(link = "log"))
  
  N0 <- exp(coef(m)[1])
  r <- exp(coef(m)[2])-1
  
  result <- c(N0, r)
  
  return(result)
  
}

# Calculate trends to get initial abundance and growth rate for each species
recovered_trends <- model_output %>%
  group_by(Species) %>%
  arrange(time) %>%
  summarise(N0_true = first(Population),
            N0_hat = get_ests(Population, time)[1],
            r_hat = get_ests(Population, time)[2]) %>%
  ungroup() %>%
  mutate(log_N0_hat = log(N0_hat)) %>%
  mutate(Species_numeric = as.numeric(as.character(Species))) %>%
  arrange(Species_numeric) %>%
  left_join(spp_params_hash, by = "Species_numeric")


# ---- Type II regr ----

m <- lmodel2(r_hat ~ log_N0_hat, 
             data = recovered_trends, 
             range.y = "interval", 
             range.x = "interval", 
             nperm = 99)

newdat <- data.frame(
  N0_hat = seq(min(recovered_trends$N0_hat), max(recovered_trends$N0_hat), length.out = 1000)) %>%
  mutate(log_N0_hat = log(N0_hat))

m_fit <- predict(
  object = m, newdata = newdat, 
  method = "RMA",  interval = "confidence",
  level = 0.95) %>%
  cbind(newdat)


# ---- Plot the results ----

# Plot species N0_hat and r_hat
ggplot(recovered_trends, aes(x = N0_hat, y = r_hat)) +
  geom_hline(yintercept = 0, col = 1) +
  geom_point(aes(color = initNfromK), size = 3, alpha = 1) +
  geom_ribbon(data = m_fit, aes(x = N0_hat, ymin = lwr, ymax = upr, y = fit), 
              fill = "black", alpha = 0.5) +
  geom_line(data = m_fit, aes(x = N0_hat, y = fit), color = "black", linewidth = 1.25) +
  labs(x = "Initial population density", 
       y = "Annual growth rate") +
  scale_x_log10(labels = scales::comma) +
  scale_y_continuous(labels = scales::percent_format(),
                     breaks = c(0.05,0.03,0.01,-0.01,-0.03,-0.05)) +
  scale_color_viridis_c("Initial distance from carrying capacity", option = "H") +
  theme_light(15) +
  guides(color = guide_colorbar(title.position = "bottom", title.hjust = 0.5, 
                                ticks.colour = NA)) +
  theme(aspect.ratio = 1,
        text = element_text(family = "Roboto"),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(angle = 0, vjust = 1),
        legend.key.height = unit(0.3, "cm"),  # Thinner bar
        legend.key.width = unit(2, "cm"))

