
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
PwWaAirV02 = function(t, state, parms){
  
  # Bioreactor parameters
  Vw <- 100 # cm3 water volume
  Vpw <- 4 # cm3 porewater volume
  Va <- 125 # cm3 headspace volumne
  Aaw <- 20 # cm2 
  Aws <- 30 # cm2
  Apw <- 1166000 # [cm2]
  
  # Pore water diffision rate
  ks <- 0.001 # [1/d]
  
  kaw.o <- 132.23 # [cm/d]
  Kaw.t <- 0.012
  
  # Bioremediation rate
  kb <- parms$kb
  
  # derivatives dx/dt are computed below
  Cpw <- state[1]
  Cw <- state[2]
  Ca <- state[3]
  
  dCpwdt <- - kb * Cpw - ks * Vw / Vpw * (Cpw - Cw)
  dCwdt <- ks * Vw / Vw * (Cpw - Cw) -
    kaw.o * Aaw / Vw * (Cw - Ca / Kaw.t) -
    kb * Cw # [ng/L]
  dCadt <- kaw.o * Aaw / Va * (Cw - Ca / Kaw.t) # Ca = [ng/L]
    
  # The computed derivatives are returned as a list
  return(list(c(dCpwdt, dCwdt, dCadt)))
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
cinit <- c(Cpw = Cpw, Cw = 0, Ca = 0) # [ng/L]
parms <- list(kb = 0.0) # Input
t <- seq(from = 0, to = 30, by = 1)
# Run the ODE function without specifying parms
out.2 <- ode(y = cinit, times = t, func = PwWaAirV02, parms = parms)
head(out.2)

{
  df.2 <- as.data.frame(out.2)
  colnames(df.2) <- c("time", "Cpw", "Cw", "Ca")
  Vpw <- 4 # [cm3]
  Vw <- 100 # [cm3]
  Va <- 125 # cm3
  df.2$Mpw <- df.2$Cpw * Vpw / 1000
  df.2$Mw <- df.2$Cw * Vw / 1000
  df.2$Ma <- df.2$Ca * Va / 1000
  df.2$Mt <- df.2$Mpw + df.2$Mw + df.2$Ma # [ng]
  df.2$fpw <- df.2$Mpw / df.2$Mt # [ng]
  df.2$fw <- df.2$Mw / df.2$Mt # [ng]
  df.2$fa <- df.2$Ma / df.2$Mt # [ng]
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

