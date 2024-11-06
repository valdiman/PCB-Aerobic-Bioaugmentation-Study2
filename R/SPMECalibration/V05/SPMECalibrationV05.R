
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

