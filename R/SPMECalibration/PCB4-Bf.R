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

# Reactive transport function ---------------------------------------------
rtm.PCB4 = function(t, state, parms){
  
  # Experimental conditions
  MH2O <- 18.0152 # g/mol water molecular weight
  MCO2 <- 44.0094 # g/mol CO2 molecular weight
  MW.pcb <- 223.088 # g/mol PCB 4 molecular weight
  
  # Bioreactor parameters
  Vpw <- 25 #cm3 porewater volume 
  Vw <- 100 # cm3 water volume
  Va <- 125 # cm3 headspace volumne
  Aaw <- 20 # cm2 
  Aws <- 30 # cm2
  
  # Congener-specific constants
  Kaw <- 0.01344142 # PCB 4 dimensionless Henry's law constant @ 25 C
  Kow <- 10^(4.65) # PCB 4 octanol-water equilibrium partition coefficient
  Koa <- 10^(6.521554861) # PCB 4 octanol-air equilibrium partition coefficient
  
  # PUF constants 
  Vpuf <- 0.000029 # m3 volume of PUF
  Kpuf <- 10^(0.6366 * log10(Koa) - 3.1774)# m3/g PCB 4-PUF equilibrium partition coefficient
  d <- 0.0213*100^3 # g/m3 density of PUF
  ro <- 0.00025 # m3/d sampling rate
  
  # SPME fiber constants
  Af <- 0.138 # cm2/cm SPME area
  Vf <- 0.000000069 # L/cm SPME volume/area
  L <- 1 # cm SPME length normalization to 1 cm
  Kf <- 10^(1.06 * log10(Kow) - 1.16) # PCB 4-SPME equilibrium partition coefficient
  ko <- 100 # cm/d PCB 4 mass transfer coefficient to SPME
  
  # Sediment partitioning
  M <- 0.1 # kg/L solid-water ratio
  foc <- 0.03 # organic carbon % in particles
  K <- foc * (10^(0.94 * log10(Kow) + 0.42)) #L/kg sediment-water equilibrium partition coefficient
  
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
  kaw.o <- (1/(Kaw.a*Kaw) + (1/Kaw.w))^-1 # [m/s]
  # v) kaw, overall air-water mass transfer coefficient for PCB 4, units change
  kaw.o <- kaw.o*100*60*60*24 # [cm/d]
  
  # Bioavailability factor B
  B <- 1 / (1 + M * K )#+ Vf * L * 1000 / Vw)
  B <- 1 / (1 + M * K * (1 + kb / (kd * M * K)))
  O <- (kb / (kd * M * K))^0.5
  #B <- 
  
  # Biotransformation rate
  kb <- 1 #23 # 1/d, value changes depending on experiment 0.023 from SPME calibration
  
  # Sortion and desorption constants
  ka <- parms$ka #1/d
  kd <- parms$kd #1/d
  
  # derivatives dx/dt are computed below
  Cw <- state[1]
  mf <- state[2]
  Ca <- state[3]
  mpuf <- state[4]
  
  dCwdt <- (kaw.o * Aaw / Vw * (Ca / (Kaw) - Cw) - kb * Cw) * B # 864 to change second to days and um to m, Ca in [ng/L]
  #dCwdt <- (kaw.o * Aaw / Vw * (Ca / (Kaw) - Cw) + kd * Cw * K * M - ka * Cw - kb * Cw) / B # 864 to change second to days and um to m, Ca in [ng/L]
  dmfdt <- ko * Af /(L * 1000) * (Cw - mf / (Vf * L * Kf)) # Cw = [ng/L], mf = [ng/cmf]
  dCadt <- kaw.o * Aaw / Va * (Cw - Ca / Kaw)
  dpufdt <- ro * Ca * 1000 - ro * (mpuf / (Vpuf * d)) / (Kpuf) # Ca = [ng/L], mpuf = [ng]
  
  # The computed derivatives are returned as a list
  return(list(c(dCwdt, dmfdt, dCadt, dpufdt)))
}

# Initial conditions and run function
# Estimating Cpw (PCB 4 concentration in sediment porewater)
Ct <- 630.2023 * 5  # ng/g PCB 4 sediment concentration
foc <- 0.03 # organic carbon % in sediment
Kow <- 10^(4.65) # PCB 4 octanol-water equilibrium partition coefficient
logKoc <- 0.94 * log10(Kow) + 0.42 # koc calculation
K <- foc * 10^(logKoc) # L/kg sediment-water equilibrium partition coefficient
ds <- 900 # g/L sediment density
M <- 0.1 # kg/L solid-water ratio
Cwi <- Ct * M * 1000 / (1 + M * K)
cinit <- c(Cw = Cwi, mf = 0, Ca = 0, mpuf = 0)
parms <- list(ka = 3, kd = 0.015) # Input 
t.1 <- seq(0, 40, by = 1)  # Adjust according to your needs
# Run the ODE function without specifying parms
out.1 <- ode(y = cinit, times = t.1, func = rtm.PCB4, parms = parms)
head(out.1)

model_results <- as_tibble(out.1) %>%
  rename(mf = `mf`, mpuf = `mpuf`) %>%
  mutate(across(c(mf, mpuf, Cw), as.numeric),
         time = as.numeric(time)) %>%  # Ensure 'time' is numeric
  select(time, mf, mpuf)  # Select only the relevant columns for plotting

# Prepare model data for plotting
model_data_long <- model_results %>%
  pivot_longer(cols = c(mf, mpuf), 
               names_to = "variable", 
               values_to = "model_value") %>%
  mutate(type = "Model")

plot_data_daily <- bind_rows(
  model_data_long %>%
    rename(value = model_value) %>%
    mutate(type = "Model")
)

# Plot mf
ggplot(plot_data_daily %>% filter(variable == "mf"), aes(x = time)) +
  geom_line(data = . %>% filter(type == "Model"),
            aes(y = value, color = "Model"), linewidth = 1) +
  labs(x = "Time", y = "mf [ng/cm]") +
  scale_color_manual(values = c("Model" = "blue")) +
  theme_bw() +
  theme(legend.title = element_blank())

