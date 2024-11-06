
# Packages and libraries --------------------------------------------------
# Install packages
install.packages("readxl")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("deSolve")
install.packages("tidyr")
install.packages("gridExtra")

# Load libraries
{
  library(readxl) # read excel data
  library(dplyr) # organize data
  library(tidyr)
  library(ggplot2) # plotting
  library(gridExtra)
  library(deSolve) # solving differential equations
}

# Experimental data -------------------------------------------------------
# Read data from Excel sheets
data <- read_excel("Data/SPMECalibration.xlsx", sheet = "data")

# Calibration data --------------------------------------------------------
# Select individual congeners from datasets
pcb.ind <- "PCB84"

# Extract relevant columns from each dataset
pcbi <- data[, c("sample", "treatment", "replicate", "time", "length", pcb.ind)]

# Define the model function
model_function <- function(time, state, parameters) {
  # Unpack parameters
  ko <- parameters["ko"]
  Af <- parameters["Af"]
  Vw <- parameters["Vw"]
  Vf <- parameters["Vf"]
  L <- parameters["L"]
  Kf <- parameters["Kf"]
  Cpw <- parameters["Cpw"]
  
  # Current value of mf (mass in the sampler)
  mf <- state["mf"]
  
  # Differential equation
  dmf_dt <- ko * Af * Vw / (Vf * L * 10^6) * (Cpw - mf / (Vf * Kf))
  
  # Return the rate of change
  list(c(dmf_dt))
}

# Initial values and time vector
initial_mf <- 0  # Starting mass in sampler, adjust if needed
time <- seq(0, max(pcbi$time), by = 1)  # Time points for simulation

# Define observed data
observed_data <- pcbi %>%
  mutate(PCB = PCB84 / length) %>%
  select(time, PCB) %>%
  rename(mf = PCB)

{
  Vw = 10 # cm3
  # Congener-specific constants
  Kow <- 10^(4.65) # PCB 4 octanol-water equilibrium partition coefficient
  # SPME fiber constants
  Af <- 0.138 # cm2/cm SPME area
  Vf <- 0.000000069 # L/cm SPME volume/area
  L <- 1 # cm SPME length normalization to 1 cm
  Kf <- 10^(1.06 * log10(Kow) - 1.16) # PCB 4-SPME equilibrium partition coefficient
  Ct <- 630.2023 # ng/g PCB 4 sediment concentration
  foc <- 0.03 # organic carbon % in sediment
  Kow <- 10^(4.65) # PCB 4 octanol-water equilibrium partition coefficient
  dUow <- -21338.96 # internal energy for the transfer of octanol-water for PCB 4 (J/mol)
  R <- 8.3144 # J/(mol K) molar gas constant
  Tst <- 25 #C air temperature
  Tst.1 <- 273.15 + Tst # air and standard temperature in K, 25 C
  Tw <- 20 # C water temperature
  Tw.1 <- 273.15 + Tw
  Kow.t <- Kow*exp(-dUow/R*(1/Tw.1-1/Tst.1))
  logKoc <- 0.94 * log10(Kow.t) + 0.42 # koc calculation
  K <- foc * 10^(logKoc) # L/kg sediment-water equilibrium partition coefficient
  Cpw <- Ct / K * 1000 # [ng/L]
}

# Wrapper function to calculate residuals for optimization
objective_function <- function(ko_value) {
  # Set parameters including current ko value
  parameters <- c(ko = ko_value, Af = Af, Vw = Vw, Vf = Vf, L = L, Kf = Kf, Cpw = Cpw)
  # Initial state
  state <- c(mf = initial_mf)
  
  # Solve the ODE
  ode_output <- ode(y = state, times = time, func = model_function, parms = parameters)
  
  # Debugging: Print ode_output to inspect the results
  print(head(ode_output))
  
  # Check if the mf column exists and is not empty
  if (!"mf" %in% colnames(ode_output)) {
    stop("mf column not found in ode_output")
  }
  
  modeled_mf <- ode_output[, "mf"]
  
  # Interpolate model output to observed time points and calculate residuals
  interpolated_mf <- approx(x = ode_output[, "time"], y = modeled_mf, xout = observed_data$time)$y
  residuals <- observed_data$mf - interpolated_mf
  
  # Return sum of squared residuals
  sum(residuals^2)
}

# Perform optimization
best_fit <- optim(par = c(ko = 0.01), fn = objective_function, method = "L-BFGS-B", lower = 0)

# Print the best fit result
print(best_fit)

# Run the model with best-fit ko and plot
best_parameters <- c(ko = best_ko, Af = Af, Vw = Vw, Vf = Vf, L = L, Kf = Kf, Cpw = Cpw)
model_output <- ode(y = c(mf = initial_mf), times = time, func = model_function, parms = best_parameters)

# Plot observed vs modeled data
plot(observed_data$time, observed_data$mf, col = "blue", pch = 19, xlab = "Time", ylab = "mf",
     main = "Observed vs Modeled mf")
lines(model_output[, "time"], model_output[, "mf"], col = "red", lwd = 2)
legend("topright", legend = c("Observed", "Modeled"), col = c("blue", "red"), pch = c(19, NA), lty = c(NA, 1))






# Define the system of differential equations
SPMEUptake <- function(t, state, parameters) {
    
    Vw = 10 # cm3
  
    # Congener-specific constants
    Kow <- 10^(4.65) # PCB 4 octanol-water equilibrium partition coefficient
  
    # SPME fiber constants
    Af <- 0.138 # cm2/cm SPME area
    Vf <- 0.000000069 # L/cm SPME volume/area
    L <- 1 # cm SPME length normalization to 1 cm
    Kf <- 10^(1.06 * log10(Kow) - 1.16) # PCB 4-SPME equilibrium partition coefficient
    
    # Passive sampler rates
    ko <- parms$ko # cm/d mass transfer coefficient to SPME
    
    # Bioremediation rate
    kb <- parms$kb
    
    Cpw <- state[1]
    mf <- state[2]
    
    #dCpwdt <-  - ko * Af / (Vf * L * 1000) * (Cpw - mf / (Vf * Kf)) # Cw = [ng/L], mf = [ng/]
    dCpwdt <- 0
    dmfdt <- ko * Af * Vw / (Vf * L * 1000 * 1000) * (Cpw - mf / (Vf * Kf)) # Cw = [ng/L], mf = [ng/cmf]
    
    return(list(c(dCpwdt, dmfdt)))
}

# Initial conditions and parameters
{
  Ct <- 630.2023 # ng/g PCB 4 sediment concentration
  foc <- 0.03 # organic carbon % in sediment
  Kow <- 10^(4.65) # PCB 4 octanol-water equilibrium partition coefficient
  dUow <- -21338.96 # internal energy for the transfer of octanol-water for PCB 4 (J/mol)
  R <- 8.3144 # J/(mol K) molar gas constant
  Tst <- 25 #C air temperature
  Tst.1 <- 273.15 + Tst # air and standard temperature in K, 25 C
  Tw <- 20 # C water temperature
  Tw.1 <- 273.15 + Tw
  Kow.t <- Kow*exp(-dUow/R*(1/Tw.1-1/Tst.1))
  logKoc <- 0.94 * log10(Kow.t) + 0.42 # koc calculation
  K <- foc * 10^(logKoc) # L/kg sediment-water equilibrium partition coefficient
  Cpw <- Ct / K * 1000 # [ng/L]
  M <- 0.1 # kg/L solid-water ratio
  Cwi <- Ct * M * 1000 / (1 + M * K)
}

cinit <- c(Cpw = Cpw, mf = 0)
parms <- list(ko = 0.000003, kb = 0) # Input 
t <- seq(0, 90, by = 1)  # Adjust according to your needs
# Run the ODE function without specifying parms
out <- ode(y = cinit, times = t, func = SPMEUptake, parms = parms)
head(out)

model <- as_tibble(out) %>%
  mutate(across(c(Cpw, mf), as.numeric),
         time = as.numeric(time)) %>%
  select(time, Cpw, mf)

model$tmass <- model$Cpw * 10 / 1000 + model$mf * 1

model$fmf_1 <- model$mf / model$tmass * 100
model$fmf_2 <- model$mf / (model$Cpw * 10 / 1000) *100

# Prepare model data for plotting
model_data_long <- model %>%
  pivot_longer(cols = c(Cpw, mf), 
               names_to = "variable", 
               values_to = "model_value") %>%
  mutate(type = "Model")

p_Cpw <- ggplot(model_data_long %>% filter(variable == "Cpw"), aes(x = time)) +
  geom_line(data = . %>% filter(type == "Model"),
            aes(y = model_value, color = "Model"), linewidth = 1) +
  labs(x = "Time", y = "Cpw [ng/L]") +
  scale_color_manual(values = c("Model" = "blue")) +
  theme_bw() +
  theme(legend.title = element_blank())

# Plot mpuf
p_mf <- ggplot(model_data_long %>% filter(variable == "mf"), aes(x = time)) +
  geom_line(data = . %>% filter(type == "Model"),
            aes(y = model_value, color = "Model"), linewidth = 1) +
  labs(x = "Time", y = "mf [ng/cm]") +
  scale_color_manual(values = c("Model" = "blue")) +
  theme_bw() +
  theme(legend.title = element_blank())

# Arrange plots side by side
grid.arrange(p_Cpw, p_mf, ncol = 2)

# Filter the data based on the condition PCB4 / length < 0.12
filtered_pcbi <- pcbi %>%
  dplyr::filter(PCB84 / length < 0.1)

# Plot the filtered values of PCB4 / length against time
plot(filtered_pcbi$time, filtered_pcbi$PCB60 / filtered_pcbi$length,
     xlab = "Time", ylab = "PCB60 / Length", 
     main = "Filtered PCB4 / Length Values over Time",
     pch = 19, col = "blue")

plot(pcbi$time, pcbi$PCB84 / pcbi$length,
     xlab = "Time", ylab = "PCB84 / Length", 
     main = "Filtered PCB84 / Length Values over Time",
     pch = 19, col = "blue")

