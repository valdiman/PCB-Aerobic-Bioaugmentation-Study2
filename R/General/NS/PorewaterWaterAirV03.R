
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

# Porewater, Water & Air Model ---------------------------------------------
PwWaAirV03 = function(t, state, parms){
  
  # Bioreactor parameters
  Vw <- 100 # cm3 water volume
  Vpw <- 4 # cm3 porewater volume
  Va <- 125 # cm3 headspace volumne
  Aaw <- 20 # cm2 
  Aws <- 30 # cm2
  Apw <- 1166000 # [cm2]
  
  # PUF constants 
  Apuf <- 7.07 # cm2
  Vpuf <- 29 # cm3 volume of PUF
  d <- 0.0213 * 100^3 # g/m3 density of PUF
  Koa <- 10^(7.765196745) # PCB 19 octanol-air equilibrium partition coefficient
  Kpuf <- 10^(0.6366 * log10(Koa) - 3.1774) # PCB 19-PUF equilibrium partition coefficient [m3/g]
  Kpuf <- Kpuf * d # [La/Lpuf]
  
  # SPME fiber constants
  Af <- 0.138 # cm2/cm SPME area
  Vf <- 0.000000069 * 1000 # cm3/cm SPME volume/area
  L <- 1 # cm SPME length normalization to 1 cm
  Kow <- 10^(5.02) # PCB 19 octanol-water equilibrium partition coefficient
  dUow <- -20988.94 # internal energy for the transfer of octanol-water for PCB 19 (J/mol)
  R <- 8.3144 # J/(mol K) molar gas constant
  Tst <- 25 #C air temperature
  Tst.1 <- 273.15 + Tst # air and standard temperature in K, 25 C
  Tw <- 20 # C water temperature
  Tw.1 <- 273.15 + Tw
  Kow.t <- Kow * exp(-dUow / R * (1 / Tw.1 -  1/ Tst.1))
  Kf <- 10^(1.06 * log10(Kow.t) - 1.16) # PCB 19-SPME equilibrium partition coefficient
  
  # Sediment paritioning
  foc <- 0.03 # organic carbon % in sediment
  logKoc <- 0.94 * log10(Kow.t) + 0.42 # koc calculation
  Kd <- foc * 10^(logKoc) # L/kg sediment-water equilibrium partition coefficient
  
  # Pore water MTC
  bl <- 0.21 # cm boundary layer thickness
  ks <- 6.928 * 10^-6 * 60 * 60 * 24 / bl # [cm/d]
  ks.m.d <- ks / 100 # [m/d]
  
  kaw.o <- 132.23 # [cm/d]
  Kaw.t <- 0.012
  
  # sediment desorption
  ksed <- parms$ksed # [1/d]
  
  # Bioremediation rate
  kb <- parms$kb
  
  # Passive sampler rates
  ko <- parms$ko # cm/d mass transfer coefficient to SPME
  ro <- parms$ro # cm/d sampling rate for PUF
  
  # derivatives dx/dt are computed below
  Cs <- state[1]
  Cpw <- state[2]
  Cw <- state[3]
  Cf <- state[4]
  Ca <- state[5]
  Cpuf <- state[6]
  
  dCsdt <- - ksed * (Cs - Cpw * Kd / 1000)
  dCpwdt <- ksed * Apw / Vpw * (Cs / Kd * 1000 - Cpw) -
    kb * Cpw - ks * Aws / Vpw * (Cpw - Cw)
  dCwdt <- ks * Aws / Vw * (Cpw - Cw) -
    ko * Af * L / Vw * (Cw - Cf / Kf) -
    kaw.o * Aaw / Vw * (Cw - Ca / Kaw.t) -
    kb * Cw # [ng/L]
  dCfdt <- ko * Af / Vf * (Cw - Cf / Kf) # Cw = [ng/L], Cf = [ng/L]
  dCadt <- kaw.o * Aaw / Va * (Cw - Ca / Kaw.t) -
    ro * Apuf / Va * (Ca - Cpuf / Kpuf) # Ca = [ng/L]
  dCpufdt <- ro * Apuf / Vpuf * (Ca - Cpuf / Kpuf) # Ca = [ng/L], Cpuf = [ng/L]
  
  # The computed derivatives are returned as a list
  return(list(c(dCsdt, dCpwdt, dCwdt, dCfdt, dCadt, dCpufdt)))
}

# Initial conditions and run function
Ct <- 259.8342356 # ng/g PCB 19 sediment concentration
cinit <- c(Cs = Ct, Cpw = 0, Cw = 0, Cf = 0, Ca = 0, Cpuf = 0) # [ng/L]
parms <- list(ro = 1, ko = 0.1, ksed = 0.001, kb = 0.0) # Input
t <- seq(from = 0, to = 80, by = 1)
# Run the ODE function without specifying parms
out.1 <- ode(y = cinit, times = t, func = PwWaAirV03, parms = parms)
head(out.1)

{
  df.1 <- as.data.frame(out.1)
  colnames(df.1) <- c("time", "Cs", "Cpw", "Cw", "Cf", "Ca", "Cpuf")
  ms <- 10 # [g]
  Vpw <- 4 # [cm3]
  Vw <- 100 # [cm3]
  Va <- 125 # cm3
  Vf <- 0.000000069 # L/cm SPME
  Vpuf <- 29 # cm3 volume of PUF
  df.1$Ms <- df.1$Cs * ms
  df.1$Mpw <- df.1$Cpw * Vpw / 1000
  df.1$Mw <- df.1$Cw * Vw / 1000
  df.1$Mf <- df.1$Cf * Vf
  df.1$Ma <- df.1$Ca * Va / 1000
  df.1$Mpuf <- df.1$Cpuf * Vpuf / 1000
  df.1$Mt <- df.1$Cs * ms + df.1$Cpw * Vpw / 1000 +
    df.1$Cw * Vw / 1000 + df.1$Cf * Vf +
    df.1$Ca * Va / 1000 + df.1$Cpuf * Vpuf / 1000 # [ng]
  df.1$fs <- df.1$Cs * ms / df.1$Mt * 100 # [ng]
  df.1$fpw <- (df.1$Cpw * Vpw / 1000) / df.1$Mt * 100 # [ng]
  df.1$fw <- (df.1$Cw * Vw / 1000) / df.1$Mt * 100 # [ng]
  df.1$ff <- df.1$Cf * Vf / df.1$Mt * 100
  df.1$fa <- (df.1$Ca * Va / 1000) / df.1$Mt * 100# [ng]
  df.1$fpuf <- df.1$Cpuf * Vpuf / 1000 / df.1$Mt * 100 # [ng]
}

# Create the plot with all
ggplot(data = df.1, aes(x = time)) +
  geom_line(aes(y = fpw, color = "Porewater"), linewidth = 1) +       # Line for fpw
  geom_line(aes(y = fw, color = "Water"), linewidth = 1) +          # Line for fw
  geom_line(aes(y = fa, color = "Air"), linewidth = 1) +  # Line for Fa
  labs(title = "Fraction vs Time", 
       x = "Time", 
       y = "Fraction") +
  scale_color_manual(values = c("Porewater" = "blue", "Water" = "red",
                                "Air" = "purple"),
                     name = "Phase") +
  theme_minimal()

# Create the plot with sediment
ggplot(data = df.1, aes(x = time)) +
  geom_line(aes(y = Cpw, color = "Porewater"), linewidth = 1) +
  labs(title = "Concentration vs Time", 
       x = "Time", 
       y = "Concentration (ng/L)") +
  scale_color_manual(values = c("Porewater" = "blue"),
                     name = "Phase") +
  theme_minimal()

# Create the plot with water
ggplot(data = df.1, aes(x = time)) +
  geom_line(aes(y = Cw, color = "Water"), linewidth = 1) +
  labs(title = "Concentration vs Time", 
       x = "Time", 
       y = "Concentration (ng/L)") +
  scale_color_manual(values = c("Water" = "red"),
                     name = "Phase") +
  theme_minimal()

# Create the plot with air
ggplot(data = df.1, aes(x = time)) +
  geom_line(aes(y = Ca, color = "Air"), linewidth = 1) +  # Line for Fa
  labs(title = "Concentration vs Time", 
       x = "Time", 
       y = "Concentration (ng/L)") +
  scale_color_manual(values = c("Air" = "purple"),
                     name = "Phase") +
  theme_minimal()

