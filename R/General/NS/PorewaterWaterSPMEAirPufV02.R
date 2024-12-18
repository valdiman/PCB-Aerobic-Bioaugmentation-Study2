
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
PwWaSPMEAirPufV02 = function(t, state, parms){
  
  # Experimental conditions
  R <- 8.3144 # J/(mol K) molar gas constant
  Tst <- 25 #C air temperature
  Tst.1 <- 273.15 + Tst # air and standard temperature in K, 25 C
  Tw <- 20 # C water temperature
  Tw.1 <- 273.15 + Tw
  
  # Bioreactor parameters
  Vw <- 100 # cm3 water volume
  Vpw <- 4 # cm3 porewater volume
  Va <- 125 # cm3 headspace volumne
  Aaw <- 20 # cm2 
  Aws <- 30 # cm2
  Apw <- 1166000 # [cm2]
  
  # Congener-specific constants
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
  # Pore water diffusion rate
  ks <- 0.001 # [1/d]
  
  kaw.o <- 132.23 # [cm/d]
  Kaw.t <- 0.012
  
  # Passive sampler rates
  ko <- parms$ko # cm/d mass transfer coefficient to SPME
  ro <- parms$ro # cm/d sampling rate for PUF
  
  # Bioremediation rate
  kb <- parms$kb
  
  # derivatives dx/dt are computed below
  Cpw <- state[1]
  Cw <- state[2]
  Cf <- state[3]
  Ca <- state[4]
  Cpuf <- state[5]
  
  dCpwdt <- - kb * Cpw - ks * Vw / Vpw * (Cpw - Cw)
  dCwdt <- ks * Vw / Vw * (Cpw - Cw) -
    kaw.o * Aaw / Vw * (Cw - Ca / Kaw.t) -
    ko * Af * L / Vw * (Cw - Cf / Kf) -
    kb * Cw # [ng/L]
  dCfdt <- ko * Af / Vf * (Cw - Cf / Kf) # Cw = [ng/L], Cf = [ng/L]
  dCadt <- kaw.o * Aaw / Va * (Cw - Ca / Kaw.t) -
    ro * Apuf / Va * (Ca - Cpuf / Kpuf) # Ca = [ng/L]
  dCpufdt <- ro * Apuf / Vpuf * (Ca - Cpuf / Kpuf) # Ca = [ng/L], Cpuf = [ng/L]
    
  # The computed derivatives are returned as a list
  return(list(c(dCpwdt, dCwdt, dCfdt, dCadt, dCpufdt)))
}

# Initial conditions and run function
{
  Cs <- 259.8342356 # ng/g PCB 19 sediment concentration
  foc <- 0.03 # organic carbon % in sediment
  Kow <- 10^(5.02) # PCB 19 octanol-water equilibrium partition coefficient
  dUow <- -20988.94 # internal energy for the transfer of octanol-water for PCB 19 (J/mol)
  R <- 8.3144 # J/(mol K) molar gas constant
  Tst <- 25 #C air temperature
  Tst.1 <- 273.15 + Tst # air and standard temperature in K, 25 C
  Tw <- 20 # C water temperature
  Tw.1 <- 273.15 + Tw
  Kow.t <- Kow*exp(-dUow / R * (1 / Tw.1 -  1/ Tst.1))
  logKoc <- 0.94 * log10(Kow.t) + 0.42 # koc calculation
  Kd <- foc * 10^(logKoc) # L/kg sediment-water equilibrium partition coefficient
  Cpw <- Cs / Kd * 1000 # [ng/L]
}
cinit <- c(Cpw = Cpw, Cw = 0, Cf = 0, Ca = 0, Cpuf = 0) # [ng/L]
parms <- list(ro = 1, ko = 1, kb = 0) # Input
t <- seq(from = 0, to = 30, by = 1)
# Run the ODE function without specifying parms
out.2 <- ode(y = cinit, times = t, func = PwWaSPMEAirPufV02, parms = parms)
head(out.2)

{
  df.2 <- as.data.frame(out.2)
  colnames(df.2) <- c("time", "Cpw", "Cw", "Cf", "Ca", "Cpuf")
  Vpw <- 4 # [cm3]
  Vw <- 100 # [cm3]
  Va <- 125 # cm3
  Vf <- 0.000000069 * 1000 # [cm3/cm SPME]
  Vpuf <- 29 # [cm3 volume of PUF]
  df.2$Mpw <- df.2$Cpw * Vpw / 1000
  df.2$Mw <- df.2$Cw * Vw / 1000
  df.2$Mf <- df.2$Cf * Vf / 1000 # [ng]
  df.2$Ma <- df.2$Ca * Va / 1000 # [ng]
  df.2$Mpuf <- df.2$Cpuf * Vpuf / 1000 # [ng]
  df.2$Mt <- df.2$Mpw + df.2$Mw + df.2$Mf + df.2$Ma + df.2$Mpuf # [ng]
  df.2$fpw <- df.2$Mpw / df.2$Mt # [ng]
  df.2$fw <- df.2$Mw / df.2$Mt # [ng]
  df.2$ff <- df.2$Mf / df.2$Mt # [ng]
  df.2$fa <- df.2$Ma / df.2$Mt # [ng]
  df.2$fpuf <- df.2$Mpuf / df.2$Mt # [ng]
}

# Create the plot with all
ggplot(data = df.2, aes(x = time)) +
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
ggplot(data = df.2, aes(x = time)) +
  geom_line(aes(y = Cpw, color = "Porewater"), linewidth = 1) +
  labs(title = "Concentration vs Time", 
       x = "Time", 
       y = "Concentration (ng/L)") +
  scale_color_manual(values = c("Porewater" = "blue"),
                     name = "Phase") +
  theme_minimal()

# Create the plot with water
ggplot(data = df.2, aes(x = time)) +
  geom_line(aes(y = Cw, color = "Water"), linewidth = 1) +
  labs(title = "Concentration vs Time", 
       x = "Time", 
       y = "Concentration (ng/L)") +
  scale_color_manual(values = c("Water" = "red"),
                     name = "Phase") +
  theme_minimal()

# Create the plot with air
ggplot(data = df.2, aes(x = time)) +
  geom_line(aes(y = Ca, color = "Air"), linewidth = 1) +  # Line for Fa
  labs(title = "Concentration vs Time", 
       x = "Time", 
       y = "Concentration (ng/L)") +
  scale_color_manual(values = c("Air" = "purple"),
                     name = "Phase") +
  theme_minimal()

