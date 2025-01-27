# Biodegradation rate calaculations

# Packages and libraries --------------------------------------------------
# Install packages
install.packages("tidyverse")
install.packages("ggplot2")

# Load libraries
{
  library(ggplot2)
  library(tidyverse)
}

# Read data ---------------------------------------------------------------
df <- read.csv("Data/kbCalculations.csv")
df <- df %>%
  mutate(PCB = as.factor(PCB))

# Define a function to fit the one-phase decay model
fit_decay_model <- function(data) {
  # Initial guesses for parameters
  initial_C0 <- max(data$concentration, na.rm = TRUE)
  initial_k <- 0.1  # small positive value as a starting point
  
  # Fit the model
  nls(concentration ~ C0 * exp(-k * time),
      data = data,
      start = list(C0 = initial_C0, k = initial_k))
}

# Fit model for each PCB
models <- df %>%
  group_by(PCB) %>%
  nest() %>%
  mutate(model = map(data, fit_decay_model))

# Extract the decay constant 'k' for each PCB congener
decay_constants <- models %>%
  mutate(k_value = map_dbl(model, ~ coef(.x)["k"])) %>%
  select(PCB, k_value)

# Calculate R^2 for each fitted model
model_metrics <- models %>%
  mutate(
    k_value = map_dbl(model, ~ coef(.x)["k"]),  # Extract k
    R2 = map2_dbl(data, model, ~ {
      # Predicted values
      preds <- predict(.y, newdata = .x)
      
      # Calculate SS_res and SS_tot
      SS_res <- sum((.x$concentration - preds)^2)
      SS_tot <- sum((.x$concentration - mean(.x$concentration))^2)
      
      # Calculate R-squared
      1 - (SS_res / SS_tot)
    })
  ) %>%
  select(PCB, k_value, R2)

# Export results
write.csv(model_metrics, file = "Output/Data/General/kb.csv")

# Extract coefficients and create fitted values
fitted_data <- models %>%
  mutate(coef = map(model, coef)) %>%
  unnest(data) %>%
  mutate(fitted_concentration = map2_dbl(time, coef, ~ .y["C0"] * exp(-.y["k"] * .x)))

# Plot observed vs fitted data
p.kb <- ggplot(fitted_data, aes(x = time, y = concentration, color = PCB)) +
  geom_point() +
  geom_line(aes(y = fitted_concentration), linetype = "dashed") +
  labs(title = "One-Phase Decay Model for Each PCB",
       x = "Time", y = "Concentration") +
  theme_bw() +
  scale_color_discrete(name = "PCB")

# see plot
p.kb

# Save plot in folder
ggsave("Output/Plots/General/kb.png", plot = p.kb, width = 10,
       height = 8, dpi = 500)
