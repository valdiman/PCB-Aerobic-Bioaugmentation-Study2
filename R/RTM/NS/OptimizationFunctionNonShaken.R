# General optimization function
# From the individual codes in the folders AVL and NBH
# needs to run the Read data, Organize data, Reactive transport function sections,
# Estimating Cpw, Cinit, and t.1, and observed_data lines. 

# Install packages
install.packages("minpack.lm")
install.packages("deSolve")

# Load libraries
{
  library(minpack.lm)
  library(deSolve) # solving differential equations
}

# Optimization function
objective_function <- function(params) {
  # Extract parameters
  ro <- params[1]
  ko <- params[2]
  
  # Update parameters in the model
  parms <- list(ro = ro, ko = ko, kb = 0.3) # For PCB 4, need to add kb value
  
  # Solve ODE
  out <- ode(y = cinit, times = t.1, func = rtm.PCB19, parms = parms) # Need to change the func to the specific PCB (rtm.PCB4 or rtm.PCB19)
  out <- as.data.frame(out)
  colnames(out) <- c("time", "Cs", "Cpw", "Cw", "Cf", "Ca", "Cpuf")
  
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
  
  # Calculate the mean of the model and observed values for 'mf' and 'mpuf'
  mean_mf_model <- mean(comparison_data$mf, na.rm = TRUE)
  mean_mf_observed <- mean(comparison_data$mf_Control, na.rm = TRUE)
  mean_mpuf_model <- mean(comparison_data$mpuf, na.rm = TRUE)
  mean_mpuf_observed <- mean(comparison_data$mpuf_Control, na.rm = TRUE)
  
  # Calculate the squared differences of the means (comparison of averages)
  diff_mf <- (mean_mf_model - mean_mf_observed)^2
  diff_mpuf <- (mean_mpuf_model - mean_mpuf_observed)^2
  
  # Return the sum of squared differences between the averages
  return(diff_mf + diff_mpuf)
}

# Initial guesses for the parameters
initial_params <- c(ro = 500, ko = 5)
# Need to run code lines from code to be optimized:
# Initial concentration
# Cinitial
# t.1
# observed data

# Optimize the parameters
result <- optim(
  par = initial_params,
  fn = objective_function,
  method = "L-BFGS-B",  # Limited-memory BFGS with bounds
  lower = c(200, 0.1),  # Set lower bounds
  upper = c(700, 10)  # Set upper bounds
)

# Set scipen option to avoid scientific notation
options(scipen = 999)

# Optimized parameters
optimized_params <- result$par
# Print optimized parameters
print(optimized_params)

# Reset scipen option to default
options(scipen = 0)



