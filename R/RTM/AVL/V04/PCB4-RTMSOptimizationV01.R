# Define function to run the model and calculate R-squared
optimize_model <- function(cf, kb, ko, observed_data, t.1) {
  
  # Set the initial conditions and parameters
  Ct <- 630.2023 * cf  # ng/g PCB 4 sediment concentration
  K <- 0.03 * 10^(0.94 * log10(10^(4.65)) + 0.42) # L/kg sediment-water equilibrium partition coefficient
  Cwi <- Ct * 0.1 * 1000 / (1 + 0.1 * K)  # PCB concentration in porewater
  
  cinit <- c(Cw = Cwi, mf = 0, Ca = 0, mpuf = 0)
  parms <- list(kb = kb, ko = ko, ro = 0.00025)  # Parameters
  
  # Run the model
  out.1 <- ode(y = cinit, times = t.1, func = rtm.PCB4, parms = parms)
  
  # Convert model results to tibble and merge with observed data
  model_results <- as_tibble(out.1) %>%
    mutate(time = as.numeric(time)) %>%
    rename(mf = "mf", mpuf = "mpuf") %>%
    select(time, mf, mpuf)
  
  comparison_data <- model_results %>%
    left_join(observed_data, by = "time")
  
  # Calculate the averages of mf and mpuf
  grouped_comparison <- comparison_data %>%
    group_by(time) %>%
    summarise(
      avg_mf_model = mean(mf, na.rm = TRUE),
      avg_mf_observed = mean(mf_Control, na.rm = TRUE)
    )
  
  # Calculate R-squared for mf
  mf_r2_value <- mf_r2(grouped_comparison$avg_mf_model, grouped_comparison$avg_mf_observed)
  
  return(mf_r2_value)
}

# Create a grid of values for cf, kb, and ko
cf_values <- seq(2, 6, by = 0.5)  # Adjust the range as needed
kb_values <- seq(0.1, 6, by = 0.5)  # Adjust the range as needed
ko_values <- seq(1, 100, by = 20)  # Adjust the range as needed

# Initialize variables to store the best parameters and R^2 value
best_r2 <- -Inf
best_params <- list()

# Loop through all combinations of cf, kb, and ko
for (cf in cf_values) {
  for (kb in kb_values) {
    for (ko in ko_values) {
      
      # Calculate R-squared for the current parameter combination
      r2_value <- optimize_model(cf, kb, ko, observed_data, t.1)
      
      # Update the best parameters if the current R^2 is higher
      if (r2_value > best_r2) {
        best_r2 <- r2_value
        best_params <- list(cf = cf, kb = kb, ko = ko)
      }
    }
  }
}

# Output the best parameters and R-squared value
print(paste("Best R-squared: ", best_r2))
print(paste("Optimal cf: ", best_params$cf))
print(paste("Optimal kb: ", best_params$kb))
print(paste("Optimal ko: ", best_params$ko))


