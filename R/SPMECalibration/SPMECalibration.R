
# Packages and libraries --------------------------------------------------
# Install packages
install.packages("readxl")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("deSolve")

# Load libraries
{
  library(readxl) # read excel data
  library(dplyr) # organize data
  library(ggplot2) # plotting
  library(deSolve) # solving differential equations
}

# Experimental data -------------------------------------------------------
# Read data from Excel sheets
data <- read_excel("Data/SPMECalibration.xlsx", sheet = "data")

# Calibration data --------------------------------------------------------
# Select individual congeners from datasets
pcb.ind <- "PCB4"

# Extract relevant columns from each dataset
pcbi <- data[, c("sample", "treatment", "replicate", "time", "length", pcb.ind)]


# Define the system of differential equations
system_equations <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    Cw <- state[1]
    mf <- state[2]
    
    Af <- 0.138
    Aw <- 4.9
    Vf <- 0.000000069
    L <- 30
    Kow <- 10^(4.65)
    Kfi <- 10^(1.06 * log10(Kow) - 1.16)
    
    kb <- parameters[1]
    kf <- parameters[2]
    
    dCwdt <- - kb * Cw
    dmfdt <- (kf * Aw / L / 1000) * (Cw - mf / (Kfi * Vf))
    
    list(c(dCwdt, dmfdt))
  })
}

# Initial conditions and parameters
Ct <- 630.2023 * 4.5
foc <- 0.03
Kow <- 10^(4.65)
logKoc <- 0.94 * log10(Kow) + 0.42
Kd <- foc * 10^(logKoc)
Cpw <- Ct / Kd * 1000
initial_conditions <- c(Cw = Cpw, mf = 0)
rate_b <- 0.05
mtcf.w <- 80
times <- seq(0, 100, by = 1)

# Define parameters list
parameters <- c(rate_b, mtcf.w)

# Cost function for optimization
cost_function <- function(params) {
  parameters <- params
  solution <- ode(y = initial_conditions, times = times, func = system_equations, parms = parameters)
  model_time <- solution[, "time"]
  model_mf <- solution[, "mf"]
  model_mf_interp <- approx(model_time, model_mf, xout = pcbi$time)$y
  sse <- sum((model_mf_interp - pcbi$PCB4 / pcbi$length)^2)
  return(sse)
}

# Track optimization steps
optimization_steps <- data.frame(iteration = integer(), parameters = I(list()), cost = numeric())

# Optimization function that records each step
optimization_function <- function(par) {
  cost <- cost_function(par)
  optimization_steps <<- rbind(optimization_steps, data.frame(iteration = nrow(optimization_steps) + 1, parameters = I(list(par)), cost = cost))
  return(cost)
}

# Perform optimization
initial_guess <- c(rate_b, mtcf.w)
opt_result <- optim(par = initial_guess, fn = optimization_function, method = "L-BFGS-B", lower = c(0, 0), upper = c(1, 200))
optimized_params <- opt_result$par
final_rate_b <- optimized_params[1]
final_mtcf_w <- optimized_params[2]

# Function to run the model with given parameters
run_model <- function(params) {
  solution <- ode(y = initial_conditions, times = times, func = system_equations, parms = params)
  return(solution[, "mf"])
}

# Bootstrap sampling
set.seed(123)
n_bootstrap <- 1000
bootstrap_results <- replicate(n_bootstrap, {
  sampled_params <- optimized_params + rnorm(2, 0, 0.01)
  run_model(sampled_params)
})

# Calculate mean and confidence intervals
mean_model <- apply(bootstrap_results, 1, mean)
ci_lower <- apply(bootstrap_results, 1, quantile, probs = 0.025)
ci_upper <- apply(bootstrap_results, 1, quantile, probs = 0.975)

# Calculate R-squared and RMSE
model_mf_interp <- approx(times, mean_model, xout = pcbi$time)$y
observed_values <- pcbi$PCB4 / pcbi$length
residuals <- observed_values - model_mf_interp

# R-squared
rss <- sum(residuals^2)
tss <- sum((observed_values - mean(observed_values))^2)
r_squared <- 1 - (rss / tss)

# RMSE
rmse <- sqrt(mean(residuals^2))

# Print R-squared and RMSE
cat("R-squared:", r_squared, "\n")
cat("RMSE:", rmse, "\n")

# Prepare data for ggplot
plot_data <- data.frame(time = times, mean_model = mean_model, ci_lower = ci_lower, ci_upper = ci_upper)
experimental_data <- data.frame(time = pcbi$time, observed = pcbi$PCB4 / pcbi$length)

# Plot using ggplot2
ggplot() +
  geom_line(data = plot_data, aes(x = time, y = mean_model), color = "blue") +
  geom_ribbon(data = plot_data, aes(x = time, ymin = ci_lower, ymax = ci_upper), fill = "lightblue", alpha = 0.5) +
  geom_point(data = experimental_data, aes(x = time, y = observed), color = "red") +
  labs(x = "Time", y = "mf") +
  theme_minimal()

