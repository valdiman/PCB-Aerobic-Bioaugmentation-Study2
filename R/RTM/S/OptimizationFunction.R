
install.packages("minpack.lm")
library(minpack.lm)


# Define the objective function
objective_function <- function(params) {
  # Extract parameters
  ro <- params[1]
  ko <- params[2]
  kdf <- params[3]
  kds <- params[4]
  ka <- params[5]
  
  # Update parameters in the model
  parms <- list(ro = ro, ko = ko, kdf = kdf, kds = kds, f = 0.6, ka = ka, kb = 0)
  
  # Solve ODE
  out <- ode(y = cinit, times = t.1, func = rtm.PCB19, parms = parms) # change the function name for other congeners.
  out <- as.data.frame(out)
  colnames(out) <- c("time", "Cs", "Cw", "Cf", "Ca", "Cpuf")
  
  # Transform Cf and Cpuf to mass/cm and mass/puf
  Vf <- 0.000000069 # L/cm SPME
  Vpuf <- 29 # cm3 volume of PUF
  out$mf <- out$Cf * Vf  # [ng/cm]
  out$mpuf <- out$Cpuf * Vpuf / 1000  # [ng/puf]
  
  # Merge model results with observed data
  model_results <- out %>%
    mutate(time = as.numeric(time)) %>%
    select(time, mf, mpuf)
  
  comparison_data <- model_results %>%
    left_join(observed_data, by = "time")
  
  # Calculate the residual sum of squares for `mf` and `mpuf`
  residual_mf <- comparison_data$mf - comparison_data$mf_Control
  residual_mpuf <- comparison_data$mpuf - comparison_data$mpuf_Control
  
  # Sum of squared errors (objective value)
  sum(residual_mf^2, na.rm = TRUE) + sum(residual_mpuf^2, na.rm = TRUE)
}

# Initial guesses for the parameters
initial_params <- c(ro = 700, ko = 10, kdf = 1.2, kds = 0.0005, ka = 65)

# Optimize the parameters
result <- optim(
  par = initial_params,
  fn = objective_function,
  method = "L-BFGS-B",  # Limited-memory BFGS with bounds
  lower = c(0, 0, 0, 0, 0),  # Set lower bounds
  upper = c(100, 50, 10, 1, 500)  # Set upper bounds
)

# Optimized parameters
optimized_params <- result$par
print(optimized_params)
