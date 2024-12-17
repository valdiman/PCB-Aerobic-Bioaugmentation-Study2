
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

# Sediment, Water, SPME, Air & Puf Model ----------------------------------
SedWatAirV01 = function(t, state, parms){
    
  # Experimental conditions
  MH2O <- 18.0152 # g/mol water molecular weight
  MCO2 <- 44.0094 # g/mol CO2 molecular weight
  MW.pcb <- 257.532 # g/mol PCB 17 molecular weight
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
  Kaw <- 0.015256157 # PCB 17 dimensionless Henry's law constant @ 25 C
  dUaw <- 52590.22 # internal energy for the transfer of air-water for PCB 17 (J/mol)
  Kaw.t <- Kaw*exp(-dUaw / R * (1 / Tw.1 - 1 / Tst.1)) * Tw.1 / Tst.1
  Kow <- 10^(5.25) # PCB 17 octanol-water equilibrium partition coefficient
  dUow <-  -22888.94 # internal energy for the transfer of octanol-water for PCB 17 (J/mol)
  Kow.t <- Kow*exp(-dUow / R * (1 / Tw.1 -  1/ Tst.1))
  Koa <- 10^(7.628540115) # PCB 17 octanol-air equilibrium partition coefficient
  
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
  # ii) Kaw.w, water-side mass transfer coefficient for PCB 17. 600 is the Schmidt number of CO2 at 298 K
  Kaw.w <- V.co2.w*(SC.pcb.w/600)^(-0.5) # [m/s]
  # iii) kaw, overall air-water mass transfer coefficient for PCB 17
  kaw.o <- (1/(Kaw.a*Kaw.t) + (1/Kaw.w))^-1 # [m/s]
  # iv) kaw, overall air-water mass transfer coefficient for PCB 17, units change
  kaw.o <- kaw.o*100*60*60*24 # [cm/d]
  
  # Bioavailability factor B. Not necessary here.
  B <- 1
  
  # Bioremediation rate
  kb <- parms$kb
  
  # Sorption and desorption rates
  kdf <- parms$kdf # 1/d
  kds <- parms$kds # 1/d
  f <- parms$f # fraction
  ka <- parms$ka # 1/d
  
  # derivatives dx/dt are computed below
  Cs <- state[1]
  Cw <- state[2]
  Ca <- state[3]
  
  Cw <- Cw / B
  
  dCsdt <- - f * kdf * Cs - (1 - f) * kds * Cs + ka * Cw
  dCwdt <- - ka * Cw + f * kdf * Cs + (1 - f) * kds * Cs -
    kaw.o * Aaw / Vw * (Cw - Ca / Kaw.t) -
    kb * Cw # [ng/L]
  dCadt <- kaw.o * Aaw / Va * (Cw - Ca / Kaw.t) # Ca = [ng/L]
  
  # The computed derivatives are returned as a list
  return(list(c(dCsdt, dCwdt, dCadt)))
}

# Initial conditions and run function
{
  Ct <- 307.3052312  # ng/g PCB 17 sediment concentration
  M <- 0.1 # kg/L solid-water ratio
  Cs0 <- Ct * M * 1000 # [ng/L]
}
cinit <- c(Cs = Cs0, Cw = 0, Ca = 0)
parms <- list(kdf = 2.17, kds = 0.0315, f = 0.8,
              ka = 179, kb = 0) # Input
t <- seq(from = 0, to = 40, by = 1)
# Run the ODE function without specifying parms
out.3 <- ode(y = cinit, times = t, func = SedWatAirV01, parms = parms)
head(out.3)
  
{
  df.3 <- as.data.frame(out.3)
  colnames(df.3) <- c("time", "Cs", "Cw", "Ca")
  ms <- 10 # [g]
  M <- 0.1 # kg/L solid-water ratio
  Vw <- 100 # [cm3]
  Va <- 125 # [cm3]
  df.3$Mp <- df.3$Cs * ms / (M * 1000)
  df.3$Mw <- df.3$Cw * Vw / 1000 # [ng]
  df.3$Ma <- df.3$Ca * Va / 1000 # [ng]
  df.3$Mt <- df.3$Mp + df.3$Mw + df.3$Ma # [ng]
  df.3$fp <- df.3$Mp / df.3$Mt * 100
  df.3$fw <- df.3$Mw / df.3$Mt * 100 # < 0.5%
  df.3$fa <- df.3$Ma / df.3$Mt * 100
}

# Create the plot with sediment, water
ggplot(data = df.3, aes(x = time)) +
  geom_line(aes(y = Cs, color = "Sediment"), linewidth = 1) +
  geom_line(aes(y = Cw, color = "Water"), linewidth = 1) +
  labs(title = "Concentration vs Time", 
       x = "Time", 
       y = "Concentration (ng/L)") +
  scale_color_manual(values = c("Sediment" = "brown", "Water" = "red"),
                     name = "Phase") +
  theme_minimal()

# Create the plot with Sediment
ggplot(data = df.3, aes(x = time)) +
  geom_line(aes(y = Cs, color = "Sediment"), linewidth = 1) +
  labs(title = "Concentration vs Time", 
       x = "Time", 
       y = "Concentration (ng/L)") +
  scale_color_manual(values = c("Sediment" = "brown"),
                     name = "Phase") +
  theme_minimal()

# Create the plot with water
ggplot(data = df.3, aes(x = time)) +
  geom_line(aes(y = Cw, color = "Water"), linewidth = 1) +
  labs(title = "Concentration vs Time", 
       x = "Time", 
       y = "Concentration (ng/L)") +
  scale_color_manual(values = c("Water" = "red"),
                     name = "Phase") +
  theme_minimal()

  