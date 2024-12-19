
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
PwWaAirV06 = function(t, state, parms){
  
  # Bioreactor parameters
  Vw <- 100 # cm3 water volume
  Vpw <- 4 # cm3 porewater volume
  Va <- 125 # cm3 headspace volumne
  Aaw <- 20 # cm2 
  Aws <- 30 # cm2
  Apw <- 1166000 # [cm2]
  ms <- 10 # [g]
  n <- 0.42 # [%] porosity
  ds <- 1540 # [g/L] sediment density
  M <- ds * (1 - n) / n # [g/L]
  Vs <- ms / M * 1000 # [cm3]
  
  # Pore water MTC
  bl <- 0.21 # cm boundary layer thickness
  ks <- 6.928 * 10^-6 * 60 * 60 * 24 / bl # [cm/d]
  ks.m.d <- ks / 100 # [m/d]
  
  kaw.o <- 132.23 # [cm/d]
  Kaw.t <- 0.012
  
  # Bioremediation rate
  kb <- parms$kb # [1/d]
  
  # sediment desorption
  ksed <- parms$ksed # [1/d]
  
  # derivatives dx/dt are computed below
  Cs <- state[1]
  Cpw <- state[2]
  Cw <- state[3]
  Ca <- state[4]
  
  dCsdt <- - ksed * (Cs - Cpw) # Desorption from sediment to porewater
  dCpwdt <- ksed * Vs / Vpw * (Cs - Cpw) -
    ks * Aws / Vpw * (Cpw - Cw) -
    kb * Cpw
  dCwdt <- ks * Aws / Vw * (Cpw - Cw) -
    kaw.o * Aaw / Vw * (Cw - Ca / Kaw.t) -
    kb * Cw # [ng/L]
  dCadt <- kaw.o * Aaw / Va * (Cw - Ca / Kaw.t) # Ca = [ng/L]
    
  # The computed derivatives are returned as a list
  return(list(c(dCsdt, dCpwdt, dCwdt, dCadt)))
}

# Initial conditions and run function
{
  Ct <- 259.8342356 # ng/g PCB 19 sediment concentration
  n <- 0.42 # [%] porosity
  ds <- 1540 # [g/L] sediment density
  M <- ds * (1 - n) / n # [g/L]
  Cs0 <- Ct * M # [ng/L]
}
cinit <- c(Cs = Cs0, Cpw = 0, Cw = 0, Ca = 0) # [ng/L]
parms <- list(ksed = 0.001, kb = 0.0) # Input ksed from Koelmas et al 2010
t <- seq(from = 0, to = 80, by = 1)
# Run the ODE function without specifying parms
out.1 <- ode(y = cinit, times = t, func = PwWaAirV06, parms = parms)
head(out.1)

{
  df.1 <- as.data.frame(out.1)
  colnames(df.1) <- c("time", "Cs", "Cpw", "Cw", "Ca")
  ms <- 10 # [g]
  Vpw <- 4 # [cm3]
  Vw <- 100 # [cm3]
  Va <- 125 # cm3
  df.1$Ms <- df.1$Cs * ms / M
  df.1$Mpw <- df.1$Cpw * Vpw / 1000
  df.1$Mw <- df.1$Cw * Vw / 1000
  df.1$Ma <- df.1$Ca * Va / 1000
  df.1$Mt <- df.1$Ms + df.1$Mpw + df.1$Mw + df.1$Ma # [ng]
  df.1$fs <- df.1$Ms / df.1$Mt # [ng]
  df.1$fpw <- df.1$Mpw / df.1$Mt # [ng]
  df.1$fw <- df.1$Mw  / df.1$Mt # [ng]
  df.1$fa <- df.1$Ma / df.1$Mt # [ng]
}

# Concentration plot
ggplot(data = df.1, aes(x = time)) +
  geom_line(aes(y = Cpw, color = "Porewater"), linewidth = 1) +       # Line for fpw
  geom_line(aes(y = Cw, color = "Water"), linewidth = 1) +          # Line for fw
  geom_line(aes(y = Ca, color = "Air"), linewidth = 1) +  # Line for Fa
  labs(title = "Fraction vs Time", 
       x = "Time", 
       y = "Concentration") +
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

