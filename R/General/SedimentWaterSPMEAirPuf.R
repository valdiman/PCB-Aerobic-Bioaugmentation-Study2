
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
  
  # Congener-specific constants
  Kaw <- 0.01344142 # PCB 4 dimensionless Henry's law constant @ 25 C
  dUaw <- 49662.48 # internal energy for the transfer of air-water for PCB 4 (J/mol)
  Kow <- 10^(4.65) # PCB 4 octanol-water equilibrium partition coefficient
  dUow <-  -21338.96 # internal energy for the transfer of octanol-water for PCB 4 (J/mol)
  Koa <- 10^(6.521554861) # PCB 4 octanol-air equilibrium partition coefficient
  
  # PUF constants 
  Apuf <- 7.07 # cm2
  Vpuf <- 0.000029 # m3 volume of PUF
  Kpuf <- 10^(0.6366 * log10(Koa) - 3.1774)# m3/g PCB 4-PUF equilibrium partition coefficient
  d <- 0.0213*100^3 # g/m3 density of PUF
  
  # SPME fiber constants
  Af <- 0.138 # cm2/cm SPME area
  Vf <- 0.000000069 # L/cm SPME volume/cm
  L <- 1 # cm SPME length normalization to 1 cm
  Kf <- 10^(1.06 * log10(Kow) - 1.16) # PCB 4-SPME equilibrium partition coefficient
  
  # Sediment partitioning
  M <- 0.1 # kg/L solid-water ratio
  foc <- 0.03 # organic carbon % in particles
  Kow.t <- Kow*exp(-dUow / R * (1 / Tw.1 -  1/ Tst.1))
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
  # i) Ka.w.t, ka.w corrected by water and air temps during experiment
  Kaw.t <- Kaw*exp(-dUaw/R*(1/Tw.1-1/Tst.1))*Tw.1/Tst.1
  # ii) Kaw.a, air-side mass transfer coefficient
  Kaw.a <- V.water.air*(D.pcb.air/D.water.air)^(0.67) # [m/s]
  # iii) Kaw.w, water-side mass transfer coefficient for PCB 4. 600 is the Schmidt number of CO2 at 298 K
  Kaw.w <- V.co2.w*(SC.pcb.w/600)^(-0.5) # [m/s]
  # iv) kaw, overall air-water mass transfer coefficient for PCB 4
  kaw.o <- (1/(Kaw.a*Kaw.t) + (1/Kaw.w))^-1 # [m/s]
  # v) kaw, overall air-water mass transfer coefficient for PCB 4, units change
  kaw.o <- kaw.o*100*60*60*24 # [cm/d]
  
  # Bioavailability factor B
  B <- (Vw + M * Vw * K + Vf * L * 1000 + Va * Kaw.t) / Vw
  
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
  mf <- state[3]
  Ca <- state[4]
  mpuf <- state[5]
  
  dCsdt <- (- f * kdf * Cs - (1 - f) * kds * Cs + ka * Cw) / B # [ng/L]
  dCwdt <- (- ka * Cw + f * kdf * Cs + (1 - f) * kds * Cs -
              kaw.o * Aaw / Vw * (Cw - Ca / Kaw.t) - 
              ko * Af / (Vf * 1000) * (Cw - mf / (Vf * Kf)) -
              kb * Cw) / B # [ng/L]
  dmfdt <- (ko * Af * Vw / (Vf * 1000 * 1000) * (Cw - mf / (Vf * Kf))) / B # Cw = [ng/L], mf = [ng/cmf]
  dCadt <- (kaw.o * Aaw / Va * (Cw - Ca / Kaw.t) -
              ro * Apuf / (Vpuf * 10^6) * (Ca - mpuf / (Vpuf * d * Kpuf * 1000)))/ B # Ca = [ng/L]
  dmpufdt <- (ro * Apuf * Va/ (Vpuf * 10^9) * (Ca - mpuf / (Vpuf * d * Kpuf * 1000))) / B # Ca = [ng/L], mpuf = [ng]
  
  # The computed derivatives are returned as a list
  return(list(c(dCsdt, dCwdt, dmfdt, dCadt, dmpufdt)))
}

# Initial conditions and run function
  {
    # Estimating Cs0 (PCB 4 concentration in particles)
    Ct <- 630.2023 # ng/g PCB 4 sediment concentration
    M <- 0.1 # kg/L solid-water ratio
    Cs0 <- Ct * M * 1000 # [ng/L]
  }
cinit <- c(Cs = Cs0, Cw = 0, mf = 0, Ca = 0, mpuf = 0)
parms <- list(ro = 50000, ko = 1, kdf = 4.2, kds = 0.01, f = 0.8,
              ka = 450, kb = 2200) # Input
t <- seq(from = 0, to = 40, by = 1)
# Run the ODE function without specifying parms
out.3 <- ode(y = cinit, times = t, func = SedWatSPMEAirPufV01, parms = parms)
head(out.3)
  
{
  df.3 <- as.data.frame(out.3)
  colnames(df.3) <- c("time", "Cs", "Cw", "mf", "Ca", "mpuf")
  Vw <- 100 # [cm3]
  Va <- 125 # [cm3]
  l <- 1 # [cm]
  df.3$Mp <- df.3$Cs * Vw / 1000
  df.3$Mw <- df.3$Cw * Vw / 1000
  df.3$Mf <- df.3$mf * l
  df.3$Ma <- df.3$Ca * Va / 1000
  df.3$Mpuf <- df.3$mpuf
  df.3$Mt <- (df.3$Cs + df.3$Cw) * Vw / 1000 + df.3$Ca * Va / 1000 + df.3$mf * l + df.3$mpuf # [ng]
  df.3$fp <- df.3$Mp / df.3$Mt * 100
  df.3$fw <- df.3$Mw / df.3$Mt * 100 # ~0.5%
  df.3$ff <- df.3$Mf / df.3$Mt * 100
  df.3$fa <- df.3$Ma / df.3$Mt * 100
  df.3$fpuf <- df.3$Mpuf / df.3$Mt * 100
}

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

# Create the plot with puf
ggplot(data = df.3, aes(x = time)) +
  geom_line(aes(y = mpuf, color = "mpuf"), linewidth = 1) +
  labs(title = "mass vs Time", 
       x = "Time", 
       y = "mass (ng/puf)") +
  scale_color_manual(values = c("mpuf" = "blue"),
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

# Create the plot with SPME
ggplot(data = df.3, aes(x = time)) +
  geom_line(aes(y = mf, color = "SPME"), linewidth = 1) +
  labs(title = "Concentration vs Time", 
       x = "Time", 
       y = "Concentration (ng/cm)") +
  scale_color_manual(values = c("SPME" = "green"),
                     name = "Phase") +
  theme_minimal()



  