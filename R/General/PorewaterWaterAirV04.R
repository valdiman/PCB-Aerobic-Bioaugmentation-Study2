
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
  
  # Sediment paritioning
  Kow <- 10^(5.02) # PCB 19 octanol-water equilibrium partition coefficient
  dUow <- -20988.94 # internal energy for the transfer of octanol-water for PCB 19 (J/mol)
  R <- 8.3144 # J/(mol K) molar gas constant
  Tst <- 25 #C air temperature
  Tst.1 <- 273.15 + Tst # air and standard temperature in K, 25 C
  Tw <- 20 # C water temperature
  Tw.1 <- 273.15 + Tw
  Kow.t <- Kow * exp(-dUow / R * (1 / Tw.1 -  1/ Tst.1))
  foc <- 0.03 # organic carbon % in sediment
  logKoc <- 0.94 * log10(Kow.t) + 0.42 # koc calculation
  Kd <- foc * 10^(logKoc) # L/kg sediment-water equilibrium partition coefficient
  
  # Pore water MTC
  bl <- 0.21 # cm boundary layer thickness
  ks <- 6.928 * 10^-6 * 60 * 60 * 24 / bl / 50# [cm/d]
  ks.m.d <- ks / 100 # [m/d]
  
  kaw.o <- 132.23 /50 # [cm/d]
  Kaw.t <- 0.012
  
  # sediment desorption
  ksed <- parms$ksed # [1/d]
  
  # derivatives dx/dt are computed below
  Cs <- state[1]
  Cpw <- state[2]
  Cw <- state[3]
  Ca <- state[4]
  
  dCsdt <- - ksed * (Cs - Cpw * Kd / 1000)
  dCpwdt <- ksed * Apw / Vpw * (Cs / Kd * 1000 - Cpw) -
    ks * Aws / Vpw * (Cpw - Cw)
  dCwdt <- ks * Aws / Vw * (Cpw - Cw) -
    kaw.o * Aaw / Vw * (Cw - Ca / Kaw.t)
  dCadt <- kaw.o * Aaw / Va * (Cw - Ca / Kaw.t)
  
  # The computed derivatives are returned as a list
  return(list(c(dCsdt, dCpwdt, dCwdt, dCadt)))
}

# Initial conditions and run function
Ct <- 259.8342356 # ng/g PCB 19 sediment concentration
cinit <- c(Cs = Ct, Cpw = 0, Cw = 0, Ca = 0) # [ng/L]
parms <- list(ksed = 0.001/50) # Input
t <- seq(from = 0, to = 30, by = 1)
# Run the ODE function without specifying parms
out.1 <- ode(y = cinit, times = t, func = PwWaAirV03, parms = parms)
head(out.1)

# Process output data
df.1 <- as.data.frame(out.1)
colnames(df.1) <- c("time", "Cs", "Cpw", "Cw", "Ca")

# Mass calculations (in ng)
ms <- 10 # Mass of sediment in g
Vpw <- 4 # Porewater volume in cm3
Vw <- 100 # Water volume in cm3
Va <- 125 # Air volume in cm3

df.1$Ms <- df.1$Cs * ms # Mass in sediment (ng)
df.1$Mpw <- df.1$Cpw * Vpw / 1000 # Mass in porewater (ng)
df.1$Mw <- df.1$Cw * Vw / 1000 # Mass in water (ng)
df.1$Ma <- df.1$Ca * Va / 1000 # Mass in air (ng)

# Total mass in the system (ng)
df.1$Mt <- df.1$Ms + df.1$Mpw + df.1$Mw + df.1$Ma
