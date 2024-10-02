
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
  library(zoo)
}

# Experimental data -------------------------------------------------------
# Read data from Excel sheets
data <- read_excel("Data/SPMECalibration.xlsx", sheet = "data")

# Calibration data --------------------------------------------------------
# Select individual congeners from datasets
pcb.ind <- "PCB10"

# Extract relevant columns from each dataset
pcbi <- data[, c("sample", "treatment", "replicate", "time", "length", pcb.ind)]

# 2nd version -------------------------------------------------------------
# Define constants outside the function
# Set up experiment
Aw <- 4.9  # [cm2] water area
Vw <- 30  # [ml or cm3]

# SPME
Af <- 0.138  # [cm2/cm] SPME area
Vf <- 0.000000069  # [L/cm] SPME volume/area
L <- 1  # [cm] SPME length normalization to 1 cm
Kow <- 10^(4.84)
Kfi <- 10^(1.06 * log10(Kow) - 1.16)  # [Lw/Lf]

# Calculate Cpw
Ct <- 75.33758507 # [ng/g]
foc <- 0.03
logKoc <- 0.94 * log10(Kow) + 0.42
Kd <- foc * 10^(logKoc)  # [L/kg] added a calibration factor to the Kd
Cpw <- Ct / Kd * 1000  # [ng/L]

# Define the system of differential equations with Cw = Cpw
system_equations <- function(t, state, parameters) {
  mf <- state[1]
  kf <- parameters[1]  # MTC from fiber/water
  
  dmfdt <- (kf * Aw / L / 1000) * (Cpw - mf / (Kfi * Vf * L))  # [ng/cmf/d]
  
  # Return the derivative for mf
  list(c(dmfdt))
}

# Initial conditions and parameters
initial_conditions <- c(mf = 0)  # Initial condition for mf
times <- seq(0, 100, by = 1)

# Define the cost function to calculate RMSE and R-squared for a given value of kf
cost_function <- function(kf) {
  parameters <- c(kf)  # Set the parameter kf
  
  # Solve the system of equations with given parameters
  solution <- ode(y = initial_conditions, 
                  times = times, 
                  func = system_equations, 
                  parms = parameters, 
                  method = "lsoda")
  
  # Convert solution to a data frame
  solution_df <- as.data.frame(solution)
  
  # Interpolate predicted values at the observation times
  predicted <- approx(solution_df$time, solution_df$mf, xout = pcbi$time)$y
  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  
  # Calculate R-squared
  rss <- sum((observed - predicted)^2)  # Residual sum of squares
  tss <- sum((observed - mean(observed))^2)  # Total sum of squares
  r_squared <- 1 - (rss / tss)
  
  return(c(rmse = rmse, r_squared = r_squared))
}

# Observed values
observed <- pcbi$PCB10 / pcbi$length

# Initial guess for kf
initial_kf <- 0.1  # Starting value for kf

# Perform the optimization using optim
optim_result <- optim(par = initial_kf, 
                      fn = function(kf) cost_function(kf)[1], 
                      method = "L-BFGS-B", 
                      lower = 0, upper = 1)  # Define the bounds for kf

# Extract the optimized kf value and the corresponding minimum RMSE and R-squared
optimized_kf <- optim_result$par
min_rmse <- cost_function(optimized_kf)["rmse"]
r_squared <- cost_function(optimized_kf)["r_squared"]

# Print the optimized kf, minimum RMSE, and R-squared
cat("Optimized kf:", optimized_kf, "\n")
cat("Minimum RMSE:", min_rmse, "\n")
cat("R-squared:", r_squared, "\n")

# Solve with optimized kf
parameters <- c(optimized_kf)
solution <- ode(y = initial_conditions, 
                times = times, 
                func = system_equations, 
                parms = parameters, 
                method = "lsoda")
solution_df <- as.data.frame(solution)

# Prepare data for ggplot
plot_data <- data.frame(time = solution_df$time)
plot_data$mf <- approx(solution_df$time, solution_df$mf, xout = plot_data$time)$y

# Interpolate predicted values for the error calculation
predicted <- approx(plot_data$time, plot_data$mf, xout = pcbi$time)$y

# Calculate errors
errors <- observed - predicted
# Create a rolling window of size, e.g., 10
rolling_sd <- rollapply(errors, width = 10, FUN = sd, fill = NA, align = "center")

# Handle missing values by interpolation
non_na_indices <- !is.na(rolling_sd)
rolling_sd_non_na <- rolling_sd[non_na_indices]
time_non_na <- seq_along(rolling_sd)[non_na_indices]

# Interpolate missing rolling_sd values
rolling_sd_interp <- approx(time_non_na, rolling_sd_non_na, xout = seq_along(plot_data$time), rule = 2)$y

# Smooth the interpolated rolling standard deviation using a spline
smooth_spline <- smooth.spline(seq_along(rolling_sd_interp), rolling_sd_interp, spar = 0.7)
rolling_sd_smooth <- predict(smooth_spline, seq_along(plot_data$time))$y

# Calculate error margins using smoothed standard deviation
error_margin <- 1.96 * rolling_sd_smooth

# Add error bounds to plot_data
plot_data$mf_upper <- plot_data$mf + error_margin
plot_data$mf_lower <- plot_data$mf - error_margin

# Prepare experimental data for plotting
experimental_data <- data.frame(time = pcbi$time, observed = observed)

# Plot using ggplot2
p.10 <- ggplot() +
  geom_line(data = plot_data, aes(x = time, y = mf), color = "blue", size = 1.5) +
  geom_ribbon(data = plot_data, aes(x = time, ymin = mf_lower, ymax = mf_upper), alpha = 0.2) +
  geom_point(data = experimental_data, aes(x = time, y = observed), color = "red") +
  theme_bw() +
  labs(x = "Time [day]", y = "PCB 10 [ng/cm]") +
  theme(aspect.ratio = 1,
        panel.border = element_blank(),
        axis.line = element_line(size = 0.5),
        axis.line.x = element_line(color = "black"),
        axis.line.y.left = element_line(color = "black"),
        axis.line.y.right = element_blank(),
        axis.line.x.top = element_blank(),
        panel.ontop = FALSE) +
  ylim(-0.05, 0.1)

# Print the plots
print(p.10)

# Save plot in folder
ggsave("Output/Plots/SPMECalibration/PCB10V01.png", plot = p.10, width = 6,
       height = 5, dpi = 500)
