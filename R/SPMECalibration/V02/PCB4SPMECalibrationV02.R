
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
# Read data from Excel sheets
data <- read_excel("Data/SPMECalibration.xlsx", sheet = "data")

# Calibration data --------------------------------------------------------
# Select individual congeners from datasets
pcb.ind <- "PCB4"

# Extract relevant columns from each dataset
pcbi <- data[, c("sample", "treatment", "replicate", "time", "length", pcb.ind)]

# Define constants outside the function
# Set up experiment
Aw <- 4.9  # [cm2] water area
Vw <- 30  # [ml or cm3]

# SPME
Af <- 0.138  # [cm2/cm] SPME area
Vf <- 0.000000069  # [L/cm] SPME volume/area
L <- 1  # [cm] SPME length normalization to 1 cm
Kow <- 10^(4.65)
Kfi <- 10^(1.06 * log10(Kow) - 1.16)  # [Lw/Lf]

# Calculate Cpw
ms <- 3  # [g]
Ct <- 630.2023  # [ng/g]
foc <- 0.03
logKoc <- 0.94 * log10(Kow) + 0.42
Kd <- foc * 10^(logKoc) * 0.2  # [L/kg] added a calibration factor to the Kd
Cpw <- Ct / Kd * 1000  # [ng/L]

# Calculate sediment-water MTC
MCO2 <- 44.0094 # g/mol CO2 molecular weight
MW.pcb <- 257.532 # g/mol PCB 19 molecular weight
D.co2.w <- 1.67606E-05 # cm2/s CO2's diffusion coefficient in water @ Tair = 25 C, patm = 1013.25 mbars 
D.pcb.water <- D.co2.w * (MW.pcb / MCO2)^(-0.5) # cm2/s PCB 19's diffusion coefficient in water @ Tair = 25 C, patm = 1013.25 mbars
bl <- 0.2 # cm boundary layer thickness
ks <- D.pcb.water / bl # [cm/s]
ks.m.d <- ks * 60 * 60 * 24 / 100 # [m/d]

# Define the system of differential equations with Cw = Cpw
system_equations <- function(t, state, parameters) {
  mf <- state[1]
  Cw <- state[2]
  Cpw <- state[3]
  
  kf <- parameters[1]  # MTC from fiber/water
  kb <- parameters[2] # biotransformation rate
  
  dmfdt <- (kf * Aw / L / 1000) * (Cw - mf / (Kfi * Vf * L))  # [ng/cmf/d]
  dCwdt <- ks * Aw / Vw * (60 * 60 * 24) * (Cpw - Cw) - kb * Cw # [ng/L/d]
  dCpwdt <- - kb * Cpw # [ng/L/d]
  
  # Return the derivative for mf and Cw
  list(c(dmfdt, dCwdt, dCpwdt))
}

# Initial conditions and parameters
initial_conditions <- c(mf = 0, Cw = 0, Cpw = Cpw)  # Initial condition for mf and Cw
times <- seq(0, 100, by = 1)

# Define the cost function to optimize kb and kf
cost_function <- function(params) {
  kf <- params[1]
  kb <- params[2]
  parameters <- c(kf, kb)
  
  # Solve the system of equations with given parameters
  solution <- ode(y = initial_conditions, 
                  times = times, 
                  func = system_equations, 
                  parms = parameters, 
                  method = "lsoda")
  
  # Convert solution to a data frame
  solution_df <- as.data.frame(solution)
  
  # Interpolate predicted values at the observation times
  predicted <- approx(solution_df$time, solution_df$mf, xout = pcbi$time, rule = 2)$y
  
  # Calculate RMSE
  rmse <- sqrt(mean((observed - predicted)^2))
  
  # Calculate R-squared
  rss <- sum((observed - predicted)^2)  # Residual sum of squares
  tss <- sum((observed - mean(observed))^2)  # Total sum of squares
  r_squared <- 1 - (rss / tss)
  
  # Return RMSE and R-squared
  return(c(rmse = rmse, r_squared = r_squared))
}

# Observed values
observed <- pcbi$PCB4 / pcbi$length

# Initial guesses for kb and kf
initial_params <- c(kf = 0.1, kb = 0.1)

# Perform the optimization using optim
optim_result <- optim(par = initial_params, 
                      fn = function(params) cost_function(params)[1], 
                      method = "L-BFGS-B", 
                      lower = c(0, 0), upper = c(1, 1))  # Bounds for kb and kf

# Extract the optimized parameters
optimized_kf <- optim_result$par[1]
optimized_kb <- optim_result$par[2]

# Calculate the minimum RMSE and R-squared using optimized parameters
min_rmse <- cost_function(c(optimized_kb, optimized_kf))["rmse"]
r_squared <- cost_function(c(optimized_kb, optimized_kf))["r_squared"]

# Print the optimized kb, kf, minimum RMSE, and R-squared
cat("Optimized kb:", optimized_kb, "\n")
cat("Optimized kf:", optimized_kf, "\n")
cat("Minimum RMSE:", min_rmse, "\n")
cat("R-squared:", r_squared, "\n")

# Solve with optimized parameters
parameters <- c(optimized_kb, optimized_kf)
solution <- ode(y = initial_conditions, 
                times = times, 
                func = system_equations, 
                parms = parameters, 
                method = "lsoda")
solution_df <- as.data.frame(solution)

# Prepare data for ggplot
plot_data <- data.frame(time = solution_df$time)
plot_data$mf <- approx(solution_df$time, solution_df$mf, xout = plot_data$time)$y
plot_data$Cw <- approx(solution_df$time, solution_df$Cw, xout = plot_data$time)$y

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
p_mf <- ggplot() +
  geom_line(data = plot_data, aes(x = time, y = mf), color = "blue", size = 1.5) +
  geom_ribbon(data = plot_data, aes(x = time, ymin = mf_lower, ymax = mf_upper), alpha = 0.2) +
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
  ylim(-0.2, 1.2)

p_Cw <- ggplot() +
  geom_line(data = plot_data, aes(x = time, y = Cw), color = "green", size = 1.5) +
  theme_bw() +
  labs(x = "Time [day]", y = "Cw [ng/L]") +
  theme(aspect.ratio = 1,
        panel.border = element_blank(),
        axis.line = element_line(size = 0.5),
        axis.line.x = element_line(color = "black"),
        axis.line.y.left = element_line(color = "black"),
        axis.line.y.right = element_blank(),
        axis.line.x.top = element_blank(),
        panel.ontop = FALSE)

# Combine the two plots into one figure
grid.arrange(p_mf, p_Cw, ncol = 2)

# Print the plots
print(p_mf)

# Save plot in folder
ggsave("Output/Plots/SPMECalibration/PCB4V01.png", plot = p_mf, width = 6,
       height = 5, dpi = 500)
