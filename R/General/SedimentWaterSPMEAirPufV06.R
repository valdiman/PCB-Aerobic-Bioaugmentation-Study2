
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
SedWatSPMEAirPufV01 = function(t, state, parms){
    
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
  
  # PUF constants 
  Apuf <- 7.07 # cm2
  Vpuf <- 29 # cm3 volume of PUF
  d <- 0.0213*100^3 # g/m3 density of PUF
  Kpuf <- 10^(0.6366 * log10(Koa) - 3.1774) # PCB 17-PUF equilibrium partition coefficient [m3/g]
  Kpuf <- Kpuf * d
  
  # SPME fiber constants
  Af <- 0.138 # cm2/cm SPME area
  Vf <- 0.000000069 * 1000 # cm3/cm SPME volume/area
  L <- 1 # cm SPME length normalization to 1 cm
  Kf <- 10^(1.06 * log10(Kow.t) - 1.16) # PCB 17-SPME equilibrium partition coefficient
  
  # Sediment partitioning
  M <- 0.1 # kg/L solid-water ratio
  foc <- 0.03 # organic carbon % in particles
  logKoc <- 0.94 * log10(Kow.t) + 0.42 # koc calculation
  K <- foc * 10^(logKoc) # L/kg sediment-water equilibrium partition coefficient
  
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
  
  # Bioavailability factor B
  B <- (Vw + M * Vw * K + Vf * Kf * L) / Vw
  
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
  
  Cw <- Cw / B
  
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
  Ct <- 307.3052312  # ng/g PCB 17 sediment concentration
  M <- 0.1 # kg/L solid-water ratio
  Cs0 <- Ct * M * 1000 # [ng/L]
}
cinit <- c(Cs = Cs0, Cw = 0, Cf = 0, Ca = 0, Cpuf = 0)
parms <- list(ro = 100, ko = 5, kdf = 2.17, kds = 0.0315, f = 0.8,
              ka = 179, kb = 0) # Input
t <- seq(from = 0, to = 40, by = 1)
# Run the ODE function without specifying parms
out.3 <- ode(y = cinit, times = t, func = SedWatSPMEAirPufV01, parms = parms)
head(out.3)
  
{
  df.3 <- as.data.frame(out.3)
  colnames(df.3) <- c("time", "Cs", "Cw", "Cf", "Ca", "Cpuf")
  ms <- 10 # [g]
  M <- 0.1 # kg/L solid-water ratio
  Vw <- 100 # [cm3]
  Va <- 125 # [cm3]
  Vf <- 0.000000069 * 1000 # cm3/cm SPME
  Vpuf <- 29 # cm3 volume of PUF
  
  df.3$Mp <- df.3$Cs * ms / (M * 1000)
  df.3$Mw <- df.3$Cw * Vw / 1000 # [ng]
  df.3$Mf <- df.3$Cf * Vf / 1000 # [ng]
  df.3$Ma <- df.3$Ca * Va / 1000 # [ng]
  df.3$Mpuf <- df.3$Cpuf * Vpuf / 1000 # [ng]
  df.3$Mt <- df.3$Mp + df.3$Mw + df.3$Mf + df.3$Ma + df.3$Mpuf # [ng]
  df.3$fp <- df.3$Mp / df.3$Mt * 100
  df.3$fw <- df.3$Mw / df.3$Mt * 100 # < 0.5%
  df.3$ff <- df.3$Mf / df.3$Mt * 100
  df.3$fa <- df.3$Ma / df.3$Mt * 100
  df.3$fpuf <- df.3$Mpuf / df.3$Mt * 100
}

# See Mf and Mpuf
df.3[1:10, c("Mf", "Mpuf")]

# Create the plot with puf
ggplot(data = df.3, aes(x = time)) +
  geom_line(aes(y = Mpuf, color = "mpuf"), linewidth = 1) +
  labs(title = "mass vs Time", 
       x = "Time", 
       y = "mass (ng/puf)") +
  scale_color_manual(values = c("mpuf" = "blue"),
                     name = "Phase") +
  theme_minimal()

# Create the plot with SPME
ggplot(data = df.3, aes(x = time)) +
  geom_line(aes(y = Mf, color = "SPME"), linewidth = 1) +
  labs(title = "Concentration vs Time", 
       x = "Time", 
       y = "Concentration (ng/cm)") +
  scale_color_manual(values = c("SPME" = "green"),
                     name = "Phase") +
  theme_minimal()

# Create the plot with Sediment
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

  