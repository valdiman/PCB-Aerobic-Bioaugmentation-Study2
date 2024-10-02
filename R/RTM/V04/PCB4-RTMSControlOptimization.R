
# PCB4-RTMSControlV05 needs to be run before

# Objective function to minimize (maximize combined R-squared for mf and mpuf)
objective_function <- function(parms_to_optimize) {
  # Extract ka and kd from parameters to optimize
  ka_opt <- parms_to_optimize[1]
  kd_opt <- parms_to_optimize[2]
  
  # Update parameters
  parms <- list(ka = ka_opt, kd = kd_opt)
  
  # Run the ODE model with the new parameters
  out_opt <- ode(y = cinit, times = t.1, func = rtm.PCB4, parms = parms)
  
  # Convert model results to tibble and merge with observed data
  model_results_opt <- as_tibble(out_opt) %>%
    mutate(time = as.numeric(time)) %>%
    rename(mf = "mf", mpuf = "mpuf") %>%
    select(time, mf, mpuf)
  
  comparison_data_opt <- model_results_opt %>%
    left_join(observed_data, by = "time")
  
  # Calculate R-squared values for mf and mpuf
  mf_r2_value_opt <- mf_r2(comparison_data_opt$mf, comparison_data_opt$mf_Control)
  mpuf_r2_value_opt <- mf_r2(comparison_data_opt$mpuf, comparison_data_opt$mpuf_Control)
  
  # Combine R-squared values (e.g., using their sum or average)
  combined_r2 <- (mf_r2_value_opt + mpuf_r2_value_opt) / 2
  
  # Since optim minimizes, we return the negative of combined R-squared
  return(-combined_r2)
}

# Initial guesses for ka and kd
start_values <- c(ka = 10, kd = 0.00000004)

# Optimize ka and kd
opt_results <- optim(par = start_values, fn = objective_function, method = "L-BFGS-B",
                     lower = c(0.001, 0.00000000000001), upper = c(100, 0.01))

# Print optimized values
opt_ka <- opt_results$par[1]
opt_kd <- opt_results$par[2]
cat("Optimized ka:", opt_ka, "\n")
cat("Optimized kd:", opt_kd, "\n")
