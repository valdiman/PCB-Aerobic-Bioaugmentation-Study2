# Code to model PCB 4 in laboratory experiments
# using sediment from Altavista, VI. Passive measurements
# of PCB 4 in the water and the air phases are predicted and
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
  exp.data <- read.csv("Data/06_Dataset_final_PCBmass.csv")
  # Select individual congener from datasets
  pcb.ind <- "PCB_4"
  # Extract relevant columns
  pcbi <- exp.data[, c("ID", "Group", "time", "Sample_medium", pcb.ind)]
}

# Organize data -----------------------------------------------------------
{
  # Remove lost sample(s), NA
  pcbi <- pcbi[!is.na(pcbi$PCB_4), ]
  # Pull congener-specific data from the dataset without averaging
  # Select SPME control samples
  pcbi.spme.control <- pcbi %>%
    filter(ID == "AVL_S", Group == "Control", Sample_medium == "SPME") %>%
    rename("mf_Control" = PCB_4)
  # Select SPME treatment samples
  pcbi.spme.treatment <- pcbi %>%
    filter(ID == "AVL_S", Group == "Treatment", Sample_medium == "SPME") %>%
    rename("mf_Treatment" = PCB_4)
  # Select PUF control samples
  pcbi.puf.control <- pcbi %>%
    filter(ID == "AVL_S", Group == "Control", Sample_medium == "PUF") %>%
    rename("mpuf_Control" = PCB_4)
  # Select PUF treatment samples
  pcbi.puf.treatment <- pcbi %>%
    filter(ID == "AVL_S", Group == "Treatment", Sample_medium == "PUF") %>%
    rename("mpuf_Treatment" = PCB_4)
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
rtm.PCB4 = function(t, state, parms){
  
  # Experimental conditions
  MH2O <- 18.0152 # g/mol water molecular weight
  MCO2 <- 44.0094 # g/mol CO2 molecular weight
  MW.pcb <- 223.088 # g/mol PCB 4 molecular weight
  R <- 8.3144 # J/(mol K) molar gas constant
  Tst <- 25 #C air temperature
  Tst.1 <- 273.15 + Tst # air and standard temperature in K, 25 C
  Tw <- 20 # C water temperature
  Tw.1 <- 273.15 + Tw
  
  # Bioreactor parameters
  Vw <- 100 # cm3 water volume
  Va <- 125 # cm3 headspace volumne
  Aaw <- 20 # cm2 
  Aws <- 30 # cm2
  
  # Congener-specific constants
  Kaw <- 0.01344142 # PCB 4 dimensionless Henry's law constant @ 25 C
  dUaw <- 49662.48 # internal energy for the transfer of air-water for PCB 4 (J/mol)
  Kaw.t <- Kaw*exp(-dUaw / R * (1 / Tw.1 - 1 / Tst.1)) * Tw.1 / Tst.1
  Kow <- 10^(4.65) # PCB 4 octanol-water equilibrium partition coefficient
  dUow <-  -21338.96 # internal energy for the transfer of octanol-water for PCB 4 (J/mol)
  Kow.t <- Kow*exp(-dUow / R * (1 / Tw.1 -  1/ Tst.1))
  Koa <- 10^(6.789298102) # PCB 4 octanol-air equilibrium partition coefficient
  
  # PUF constants
  Apuf <- 7.07 # cm2
  Vpuf <- 29 # cm3 volume of PUF
  d <- 21300 # g/m3 density of PUF
  Kpuf <- 10^(0.6366 * log10(Koa) - 3.1774) # PCB 4-PUF equilibrium partition coefficient [m3/g]
  Kpuf <- Kpuf * d # [La/Lpuf]
  
  # SPME fiber constants
  Af <- 0.138 # cm2/cm SPME area
  Vf <- 0.000000069 * 1000 # cm3/cm SPME volume/area
  L <- 1 # cm SPME length normalization to 1 cm
  Kf <- 10^(1.06 * log10(Kow.t) - 1.16) # PCB 4-SPME equilibrium partition coefficient [Lf/Lw]
  
  # Air & water physical conditions
  D.water.air <- 0.2743615 # cm2/s water's diffusion coefficient in the gas phase @ Tair = 25 C, patm = 1013.25 mbars 
  D.co2.w <- 1.67606E-05 # cm2/s CO2's diffusion coefficient in water @ Tair = 25 C, patm = 1013.25 mbars 
  D.pcb.air <- D.water.air*(MW.pcb/MH2O)^(-0.5) # cm2/s PCB 4's diffusion coefficient in the gas phase (eq. 18-45)
  D.pcb.water <- D.co2.w*(MW.pcb/MCO2)^(-0.5) # cm2/s PCB 4's diffusion coefficient in water @ Tair = 25 C, patm = 1013.25 mbars
  v.H2O <- 0.010072884	# cm2/s kinematic viscosity of water @ Tair = 25
  V.water.air <- 0.003 # m/s water's velocity of air-side mass transfer without ventilation (eq. 20-15)
  V.co2.w <- 4.1*10^-2 # m/s mass transfer coefficient of CO2 in water side without ventilation
  SC.pcb.w <- v.H2O/D.pcb.water # Schmidt number PCB 4

  # kaw calculations (air-water mass transfer coefficient)
  # i) Kaw.a, air-side mass transfer coefficient
  Kaw.a <- V.water.air*(D.pcb.air/D.water.air)^(0.67) # [m/s]
  # ii) Kaw.w, water-side mass transfer coefficient for PCB 4. 600 is the Schmidt number of CO2 at 298 K
  Kaw.w <- V.co2.w*(SC.pcb.w/600)^(-0.5) # [m/s]
  # iii) kaw, overall air-water mass transfer coefficient for PCB 4
  kaw.o <- (1/(Kaw.a*Kaw.t) + (1/Kaw.w))^-1 # [m/s]
  # iv) kaw, overall air-water mass transfer coefficient for PCB 4, units change
  kaw.o <- kaw.o*100*60*60*24 # [cm/d]
  
  # Bioremediation rate
  kb <- parms$kb
  
  # Sorption and desorption rates
  kdf <- parms$kdf # 1/d
  kds <- parms$kds # 1/d
  f <- parms$f # fraction
  ka <- parms$ka # 1/d
  
  # Passive sampler rates
  ko <- parms$ko # cm/d mass transfer coefficient to SPME
  ro <- parms$ro # cm/d sampling rate for PUF
  
  # derivatives dx/dt are computed below
  Cs <- state[1]
  Cw <- state[2]
  Cf <- state[3]
  Ca <- state[4]
  Cpuf <- state[5]
  
  dCsdt <- - f * kdf * Cs - (1 - f) * kds * Cs + ka * Cw
  dCwdt <- - ka * Cw + f * kdf * Cs + (1 - f) * kds * Cs -
    kaw.o * Aaw / Vw * (Cw - Ca / Kaw.t) - 
    ko * Af * L / Vw * (Cw - Cf / Kf) -
    kb * Cw # [ng/L]
  dCfdt <- ko * Af / Vf * (Cw - Cf / Kf) # Cw = [ng/L], Cf = [ng/L]
  dCadt <- kaw.o * Aaw / Va * (Cw - Ca / Kaw.t) -
    ro * Apuf / Va * (Ca - Cpuf / Kpuf) # Ca = [ng/L]
  dCpufdt <- ro * Apuf / Vpuf * (Ca - Cpuf / Kpuf) # Ca = [ng/L], Cpuf = [ng/L]
  
  # The computed derivatives are returned as a list
  return(list(c(dCsdt, dCwdt, dCfdt, dCadt, dCpufdt)))
}

# Initial conditions and run function
{
  # Estimating Cs0 (PCB 4 concentration in particles)
  Ct <- 630.2023 # ng/g PCB 4 sediment concentration
  M <- 0.1 # kg/L solid-water ratio
  Cs0 <- Ct * M * 1000 # [ng/L]
}
cinit <- c(Cs = Cs0, Cw = 0, Cf = 0, Ca = 0, Cpuf = 0)
parms <- list(ro = 540.409, ko = 10, kdf = 3.6, kds = 0.001, f = 0.8,
              ka = 90, kb = 1) # Input
t.1 <- unique(pcb_combined_control$time)
# Run the ODE function without specifying parms
out.1 <- ode(y = cinit, times = t.1, func = rtm.PCB4, parms = parms)
head(out.1)

{
# Transform Cf and Cpuf to mass/cm and mass/puf
out.1 <- as.data.frame(out.1)
colnames(out.1) <- c("time", "Cs", "Cw", "Cf", "Ca", "Cpuf")

# Calculate masses based on volumes
{
  Vw <- 100 # cm3
  Va <- 125 # cm3
  Vf <- 0.000000069 # L/cm, SPME volume/area
  Vpuf <- 29 # cm3, volume of PUF
  out.1$ms <- out.1$Cs * Vw / 1000
  out.1$mw <- out.1$Cw * Vw / 1000
  out.1$ma <- out.1$Ca * Va / 1000
  out.1$mf <- out.1$Cf * Vf # [ng/cm]
  out.1$mpuf <- out.1$Cpuf * Vpuf / 1000  # [ng/puf]
  out.1$Mt <- out.1$ms + out.1$mw + out.1$ma + out.1$mf + out.1$mpuf
}

# Ensure observed data is in a tibble
observed_data <- as_tibble(pcb_combined_control) %>%
  select(time, mf_Control, mpuf_Control)

# Convert model results to tibble and select relevant columns
model_results <- as_tibble(out.1) %>%
  mutate(time = as.numeric(time)) %>%
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
cinit <- c(Cs = Cs0, Cw = 0, Cf = 0, Ca = 0, Cpuf = 0)
t_daily <- seq(0, 40, by = 1)  # Adjust according to your needs
out_daily <- ode(y = cinit, times = t_daily, func = rtm.PCB4, parms = parms)
head(out_daily)

# Transform Cf and Cpuf to mass/cm and mass/puf
out.daily <- as.data.frame(out_daily)
colnames(out.daily) <- c("time", "Cs", "Cw", "Cf", "Ca", "Cpuf")

# Calculate Mf and Mpuf based on volumes
out.daily$mf <- out.daily$Cf * Vf # [ng]
out.daily$mpuf <- out.daily$Cpuf * Vpuf / 1000  # [ng]

# Convert model results to tibble and ensure numeric values
model_results_daily_clean <- as_tibble(out.daily) %>%
  mutate(across(c(mf, mpuf, Cw, time), as.numeric)) %>%  # Ensure all relevant columns are numeric
  select(time, mf, mpuf)

# Export data
write.csv(model_results_daily_clean, file = "Output/Data/RTM/S/AVL/PCB4AVLControlFV.csv")

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
}

# Arrange plots side by side
p.4 <- grid.arrange(p_mf, p_mpuf, ncol = 2)

# Save plot in folder
ggsave("Output/Plots/RTM/S/AVL/PCB4ALV_S_ControlFV.png", plot = p.4, width = 15,
       height = 5, dpi = 500)


