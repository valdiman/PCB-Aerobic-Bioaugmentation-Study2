# Packages and libraries --------------------------------------------------
# Install packages
install.packages("readxl")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("deSolve")
install.packages("gridExtra")

# Load libraries
{
  library(readxl) # read excel data
  library(dplyr) # organize data
  library(ggplot2) # plotting
  library(deSolve) # solving differential equations
  library(zoo)
  library(gridExtra)
}

# Experimental data -------------------------------------------------------
data <- read_excel("Data/SPMECalibration.xlsx", sheet = "data")

# Calibration data --------------------------------------------------------
pcb.ind <- "PCB4"
pcbi <- data[, c("sample", "treatment", "replicate", "time", "length", pcb.ind)]

# Define the objective function for optimization
objective_function <- function(params) {
  # Extract parameters
  cf <- params[1]
  kb <- params[2]
  
  # Set up experiment parameters
  Aw <- 4.9  # [cm2] water area
  Vw <- 30  # [ml or cm3]
  Af <- 0.138  # [cm2/cm] SPME area
  Vf <- 0.000000069  # [L/cm] SPME volume/area
  L <- 1  # [cm] SPME length normalization to 1 cm
  Kow <- 10^(4.65)
  Kfi <- 10^(1.06 * log10(Kow) - 1.16)  # [Lw/Lf]
  kf <- 10  # [cm/d] MTC from fiber/water
  
  # Recalculate Cpw with the new cf value
  Ct <- 630.2023 * cf  # [ng/g]
  foc <- 0.03
  logKoc <- 0.94 * log10(Kow) + 0.42
  Kd <- foc * 10^(logKoc)  # [L/kg]
  Cpw <- Ct / Kd * 1000  # [ng/L]
  
  # Define the system of differential equations
  system_equations <- function(t, state, parms) {
    mf <- state[1]
    Cpw <- state[2]
    
    dmfdt <- (kf * Aw / L / 1000) * (Cpw - mf / (Kfi * Vf * L))  # [ng/cmf/d]
    dCpwdt <- - kb * Cpw # [ng/L/d]
    
    list(c(dmfdt, dCpwdt))
  }
  
  # Initial conditions
  initial_conditions <- c(mf = 0, Cpw = Cpw)  # Initial conditions for mf and Cw
  times <- seq(0, 100, by = 1)
  
  # Solve the system of equations
  solution <- ode(y = initial_conditions, times = times, func = system_equations, parms = NULL, method = "lsoda")
  solution_df <- as.data.frame(solution)
  
  # Interpolate predicted values at the observation times
  observed <- pcbi$PCB4 / pcbi$length
  predicted <- approx(solution_df$time, solution_df$mf, xout = pcbi$time)$y
  
  # Calculate RSS, TSS, and R-squared
  rss <- sum((observed - predicted)^2)  # Residual sum of squares
  tss <- sum((observed - mean(observed))^2)  # Total sum of squares
  r_squared <- 1 - (rss / tss)  # R-squared
  
  return(-r_squared)  # Return negative R-squared (optim minimizes the objective function)
}

# Initial parameters for optimization
initial_params <- c(cf = 1, kb = 0.03)

# Run the optimization
opt_results <- optim(par = initial_params, fn = objective_function, method = "L-BFGS-B", 
                     lower = c(0, 0), upper = c(5, 1))

# Print the results
cat("Optimized parameters:\n")
cat("cf:", opt_results$par[1], "\n")
cat("kb:", opt_results$par[2], "\n")
cat("Optimized R-squared:", -opt_results$value, "\n")  # Note: value is negated R-squared

# Use the optimized parameters for final predictions
cf_optimized <- opt_results$par[1]
kb_optimized <- opt_results$par[2]

# Recalculate Cpw with the optimized cf value
# Set up experiment
Aw <- 4.9  # [cm2] water area
Vw <- 30  # [ml or cm3]

# SPME
Af <- 0.138  # [cm2/cm] SPME area
Vf <- 0.000000069  # [L/cm] SPME volume/area
L <- 1  # [cm] SPME length normalization to 1 cm
Kow <- 10^(4.65)
Kfi <- 10^(1.06 * log10(Kow) - 1.16)  # [Lw/Lf]
kf <- 10  # [cm/d] MTC from fiber/water

# Recalculate Cpw with the new cf value
Ct <- 630.2023 * cf_optimized  # [ng/g]
foc <- 0.03
logKoc <- 0.94 * log10(Kow) + 0.42
Kd <- foc * 10^(logKoc)  # [L/kg]
Cpw <- Ct / Kd * 1000  # [ng/L]

# Redefine the system of differential equations with optimized parameters
system_equations <- function(t, state, parms) {
  mf <- state[1]
  Cpw <- state[2]
  
  dmfdt <- (kf * Aw / L / 1000) * (Cpw - mf / (Kfi * Vf * L))  # [ng/cmf/d]
  dCpwdt <- - kb_optimized * Cpw # [ng/L/d]
  
  list(c(dmfdt, dCpwdt))
}

# Initial conditions and parameters
initial_conditions <- c(mf = 0, Cpw = Cpw)  # Initial condition for mf and Cw
times <- seq(0, 100, by = 1)  # Time steps from 0 to 100

# Solve the system of equations with optimized parameters
solution <- ode(y = initial_conditions, times = times, func = system_equations, parms = NULL, method = "lsoda")
solution_df <- as.data.frame(solution)

# Prepare data for ggplot
plot_data <- data.frame(time = solution_df$time, mf = approx(solution_df$time, solution_df$mf, xout = solution_df$time)$y)
observed <- pcbi$PCB4 / pcbi$length
experimental_data <- data.frame(time = pcbi$time, observed = observed)

# Plot the data using ggplot
ggplot() +
  geom_line(data = plot_data, aes(x = time, y = mf), color = "blue", size = 1.5) +
  geom_point(data = experimental_data, aes(x = time, y = observed), color = "red") +
  theme_bw() +
  labs(x = "Time [day]", y = "PCB 4 [ng/cm]") +
  theme(aspect.ratio = 1,
        panel.border = element_blank(),
        axis.line = element_line(size = 0.5),
        axis.line.x = element_line(color = "black"),
        axis.line.y.left = element_line(color = "black"),
        axis.line.y.right = element_blank(),
        axis.line.x.top = element_blank(),
        panel.ontop = FALSE) +
  ylim(-0.3, 1.2)
