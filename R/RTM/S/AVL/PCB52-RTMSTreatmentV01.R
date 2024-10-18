# Code to model PCB 19 in laboratory experiments
# using sediment from Altavista, VI. Passive measurements
# of PCB 19 in the water and the air phases are predicted and
# linked to the water and air concentrations from the passive
# samplers.

# Packages and libraries --------------------------------------------------
# Install packages
install.packages("dplyr")
install.packages("reshape2")
install.packages("ggplot2")
install.packages("deSolve")
install.packages("tidyr")
install.packages("gridExtra")

# Load libraries
{
  library(dplyr) # organize data
  library(reshape2) # organize data
  library(ggplot2) # plotting
  library(deSolve) # solving differential equations
  library(tidyr)
  library(gridExtra)
}

# Read data ---------------------------------------------------------------
{
  exp.data.SPME <- read.csv("Data/SPME_data_long.csv")
  exp.data.PUF <- read.csv("Data/PUF_data_long.csv")
  # Add sampler
  exp.data.SPME$sampler <- "SPME"
  exp.data.PUF$sampler <- "PUF"
  # Combine
  exp.data <- rbind(exp.data.PUF, exp.data.SPME)
  # Select individual congener from datasets
  pcb.ind <- "PCB_52"
  # Extract relevant columns from each dataset
  pcbi <- exp.data[, c("ID", "Group", "time", "sampler", pcb.ind)]
}

# Organize data -----------------------------------------------------------
{
  # Remove lost sample(s), NA
  pcbi <- pcbi[!is.na(pcbi$PCB_52), ]
  # Pull congener-specific data from the dataset without averaging
  # Select SPME control samples
  pcbi.spme.control <- pcbi %>%
    filter(ID == "AVL_S", Group == "Control", sampler == "SPME") %>%
    rename("mf_Control" = PCB_52)
  # Select SPME treatment samples
  pcbi.spme.treatment <- pcbi %>%
    filter(ID == "AVL_S", Group == "Treatment", sampler == "SPME") %>%
    rename("mf_Treatment" = PCB_52)
  # Select PUF control samples
  pcbi.puf.control <- pcbi %>%
    filter(ID == "AVL_S", Group == "Control", sampler == "PUF") %>%
    rename("mpuf_Control" = PCB_52)
  # Select PUF treatment samples
  pcbi.puf.treatment <- pcbi %>%
    filter(ID == "AVL_S", Group == "Treatment", sampler == "PUF") %>%
    rename("mpuf_Treatment" = PCB_52)
  # Combine the mf and mpuf data for Control
  pcb_combined_control <- cbind(
    pcbi.spme.control %>%
      select(time, mf_Control),
    pcbi.puf.control %>%
      select(mpuf_Control)
  )
  # Add a row for time = 0 if needed
  pcb_combined_control <- rbind(
    data.frame(time = 0, mf_Control = 0, mpuf_Control = 0),
    pcb_combined_control
  )
  # Combine the mf and mpuf data for Treatment
  pcb_combined_treatment <- cbind(
    pcbi.spme.treatment %>%
      select(time, mf_Treatment),
    pcbi.puf.treatment %>%
      select(mpuf_Treatment)
  )
  # Add a row for time = 0 if needed
  pcb_combined_treatment <- rbind(
    data.frame(time = 0, mf_Treatment = 0, mpuf_Treatment = 0),
    pcb_combined_treatment
  )
}

# Reactive transport function ---------------------------------------------
rtm.PCB52 = function(t, state, parms){
  
  # Experimental conditions
  MH2O <- 18.0152 # g/mol water molecular weight
  MCO2 <- 44.0094 # g/mol CO2 molecular weight
  MW.pcb <- 291.976 # g/mol PCB 52 molecular weight
  
  # Bioreactor parameters
  Vw <- 100 # cm3 water volume
  Va <- 125 # cm3 headspace volumne
  Aaw <- 20 # cm2 
  Aws <- 30 # cm2
  
  # Congener-specific constants
  Kaw <- 0.0130452 # PCB 52 dimensionless Henry's law constant @ 25 C
  Kow <- 10^(5.84) # PCB 52 octanol-water equilibrium partition coefficient
  Koa <- 10^(7.724554861) # PCB 52 octanol-air equilibrium partition coefficient
  
  # PUF constants 
  Vpuf <- 0.000029 # m3 volume of PUF
  Kpuf <- 10^(0.6366 * log10(Koa) - 3.1774)# m3/g PCB 52-PUF equilibrium partition coefficient
  d <- 0.0213*100^3 # g/m3 density of PUF
  
  # SPME fiber constants
  Af <- 0.138 # cm2/cm SPME area
  Vf <- 0.000000069 # L/cm SPME volume/area
  L <- 1 # cm SPME length normalization to 1 cm
  Kf <- 10^(1.06 * log10(Kow) - 1.16) # PCB 52-SPME equilibrium partition coefficient
  
  # Sediment partitioning
  M <- 0.1 # kg/L solid-water ratio
  foc <- 0.03 # organic carbon % in particles
  K <- foc * (10^(0.94 * log10(Kow) + 0.42)) #L/kg sediment-water equilibrium partition coefficient
  
  # Air & water physical conditions
  D.water.air <- 0.2743615 # cm2/s water's diffusion coefficient in the gas phase @ Tair = 25 C, patm = 1013.25 mbars 
  D.co2.w <- 1.67606E-05 # cm2/s CO2's diffusion coefficient in water @ Tair = 25 C, patm = 1013.25 mbars 
  D.pcb.air <- D.water.air*(MW.pcb/MH2O)^(-0.5) # cm2/s PCB 31's diffusion coefficient in the gas phase (eq. 18-45)
  D.pcb.water <- D.co2.w*(MW.pcb/MCO2)^(-0.5) # cm2/s PCB 31's diffusion coefficient in water @ Tair = 25 C, patm = 1013.25 mbars
  v.H2O <- 0.010072884	# cm2/s kinematic viscosity of water @ Tair = 25
  V.water.air <- 0.003 # m/s water's velocity of air-side mass transfer without ventilation (eq. 20-15)
  V.co2.w <- 4.1*10^-2 # m/s mass transfer coefficient of CO2 in water side without ventilation
  SC.pcb.w <- v.H2O/D.pcb.water # Schmidt number PCB 52

  # kaw calculations (air-water mass transfer coefficient)
  # i) Kaw.a, air-side mass transfer coefficient
  Kaw.a <- V.water.air*(D.pcb.air/D.water.air)^(0.67) # [m/s]
  # ii) Kaw.w, water-side mass transfer coefficient for PCB 52. 600 is the Schmidt number of CO2 at 298 K
  Kaw.w <- V.co2.w*(SC.pcb.w/600)^(-0.5) # [m/s]
  # iii) kaw, overall air-water mass transfer coefficient for PCB 52
  kaw.o <- (1/(Kaw.a*Kaw) + (1/Kaw.w))^-1 # [m/s]
  # v) kaw, overall air-water mass transfer coefficient for PCB 52, units change
  kaw.o <- kaw.o*100*60*60*24 # [cm/d]
  
  # Bioavailability factor B
  B <- (Vw + M * Vw * K + Vf * L * 1000) / Vw
  
  # Passive sampler rates
  ro <- parms$ro # m3/d sampling rate for PUF
  ko <- parms$ko # cm/d mass transfer coefficient to SPME
  
  # Biotransformation, sortion and desorption rates
  kb <- parms$kb #1/d
  ka <- parms$ka #1/d
  kd <- parms$kd #1/d
  
  # derivatives dx/dt are computed below
  Cw <- state[1]
  mf <- state[2]
  Ca <- state[3]
  mpuf <- state[4]
  
  dCwdt <- (kaw.o * Aaw / Vw * (Ca / (Kaw) - Cw) + kd * Cw * K * M - ka * Cw - kb * Cw) / B # 864 to change second to days and um to m, Ca in [ng/L]
  dmfdt <- ko * Af /(L * 1000) * (Cw - mf / (Vf * L * Kf)) # Cw = [ng/L], mf = [ng/cmf]
  dCadt <- kaw.o * Aaw / Va * (Cw - Ca / Kaw)
  dpufdt <- ro * Ca * 1000 - ro * (mpuf / (Vpuf * d)) / (Kpuf) # Ca = [ng/L], mpuf = [ng]
  
  # The computed derivatives are returned as a list
  return(list(c(dCwdt, dmfdt, dCadt, dpufdt)))
}

# Initial conditions and run function
# Estimating Cpw (PCB 52 concentration in sediment porewater)
Ct <- 321.4900673 * 4.5 # ng/g PCB 52 sediment concentration
foc <- 0.03 # organic carbon % in sediment
Kow <- 10^(5.84) # PCB 52 octanol-water equilibrium partition coefficient
logKoc <- 0.94 * log10(Kow) + 0.42 # koc calculation
K <- foc * 10^(logKoc) # L/kg sediment-water equilibrium partition coefficient
M <- 0.1 # kg/L solid-water ratio
Cwi <- Ct * M * 1000 / (1 + M * K)
kb2 <- 0.12
Cwi <- Cwi * exp(-kb2 * 11) # n days?
cinit <- c(Cw = Cwi, mf = 0, Ca = 0, mpuf = 0)
parms <- list(ro = 0.0003, ko = 5, kb = 0.0, ka = 25, kd = 0.0001) # Input
t.1 <- unique(pcb_combined_control$time)
# Run the ODE function without specifying parms
out.1 <- ode(y = cinit, times = t.1, func = rtm.PCB52, parms = parms)
head(out.1)

# Ensure observed data is in a tibble
observed_data <- as_tibble(pcb_combined_control) %>%
  select(time, mf_Control, mpuf_Control)

# Convert model results to tibble, ensuring 'time' is numeric and rename columns if needed
model_results <- as_tibble(out.1) %>%
  mutate(time = as.numeric(time)) %>%
  rename(mf = "mf", mpuf = "mpuf") %>%
  select(time, mf, mpuf)

# Merge model results with observed data
comparison_data <- model_results %>%
  left_join(observed_data, by = "time")

# Calculate the averages of mf and mpuf within each group (e.g., per time)
grouped_comparison <- comparison_data %>%
  group_by(time) %>%  # Adjust the grouping variable if needed
  summarise(
    avg_mf_model = mean(mf, na.rm = TRUE),
    avg_mf_observed = mean(mf_Control, na.rm = TRUE),
    avg_mpuf_model = mean(mpuf, na.rm = TRUE),
    avg_mpuf_observed = mean(mpuf_Control, na.rm = TRUE)
  )

# Define function to calculate R-squared, handling NA values
mf_r2 <- function(predicted, observed) {
  # Remove NA values from both predicted and observed
  valid_indices <- complete.cases(predicted, observed)
  predicted <- predicted[valid_indices]
  observed <- observed[valid_indices]
  
  # Calculate R-squared
  ss_total <- sum((observed - mean(observed))^2)
  ss_residual <- sum((observed - predicted)^2)
  1 - (ss_residual / ss_total)
}

# Calculate R-squared values based on grouped average data
mf_r2_value <- mf_r2(grouped_comparison$avg_mf_model, grouped_comparison$avg_mf_observed)
mpuf_r2_value <- mf_r2(grouped_comparison$avg_mpuf_model, grouped_comparison$avg_mpuf_observed)

# Print R-squared values
print(paste("R-squared for mf (average): ", mf_r2_value))
print(paste("R-squared for mpuf (average): ", mpuf_r2_value))

# Plot
# Run the model with the new time sequence
cinit <- c(Cw = Cwi, mf = 0, Ca = 0, mpuf = 0)
t_daily <- seq(0, 40, by = 1)  # Adjust according to your needs
out_daily <- ode(y = cinit, times = t_daily, func = rtm.PCB52, parms = parms)

# Convert model results to tibble and ensure numeric values
model_results_daily_clean <- as_tibble(out_daily) %>%
  rename(mf = `mf`, mpuf = `mpuf`) %>%
  mutate(across(c(mf, mpuf, Cw), as.numeric),
         time = as.numeric(time)) %>%  # Ensure 'time' is numeric
  select(time, mf, mpuf)  # Select only the relevant columns for plotting

# Export data
write.csv(model_results_daily_clean, file = "Output/Data/RTM/S/AVL/PCB52AVLSControl.csv")

# Prepare model data for plotting
model_data_long <- model_results_daily_clean %>%
  pivot_longer(cols = c(mf, mpuf), 
               names_to = "variable", 
               values_to = "model_value") %>%
  mutate(type = "Model")

# Clean observed data and prepare for plotting
observed_data_clean <- observed_data %>%
  pivot_longer(cols = c(mf_Control, mpuf_Control), 
               names_to = "variable", 
               values_to = "observed_value") %>%
  mutate(variable = recode(variable, 
                           "mf_Control" = "mf", 
                           "mpuf_Control" = "mpuf"),
         type = "Observed")

plot_data_daily <- bind_rows(
  model_data_long %>%
    rename(value = model_value) %>%
    mutate(type = "Model"),
  observed_data_clean %>%
    rename(value = observed_value) %>%
    mutate(type = "Observed")
)

# Plot mf
p_mf <- ggplot(plot_data_daily %>% filter(variable == "mf"), aes(x = time)) +
  geom_line(data = . %>% filter(type == "Model"),
            aes(y = value, color = "Model"), linewidth = 1) +
  geom_point(data = . %>% filter(type == "Observed"),
             aes(y = value, color = "Observed"), size = 2) +
  labs(x = "Time", y = "mf [ng/cm]") +
  scale_color_manual(values = c("Model" = "blue", "Observed" = "red")) +
  theme_bw() +
  theme(legend.title = element_blank())

# Plot mpuf
p_mpuf <- ggplot(plot_data_daily %>% filter(variable == "mpuf"), aes(x = time)) +
  geom_line(data = . %>% filter(type == "Model"),
            aes(y = value, color = "Model"), linewidth = 1) +
  geom_point(data = . %>% filter(type == "Observed"),
             aes(y = value, color = "Observed"), size = 2) +
  labs(x = "Time", y = "mpuf [ng/PUF]") +
  scale_color_manual(values = c("Model" = "blue", "Observed" = "red")) +
  theme_bw() +
  theme(legend.title = element_blank())

# Arrange plots side by side
p.52 <- grid.arrange(p_mf, p_mpuf, ncol = 2)

# Save plot in folder
ggsave("Output/Plots/RTM/S/AVL/PCB52ALV_S_Control.png", plot = p.52, width = 15,
       height = 5, dpi = 500)


