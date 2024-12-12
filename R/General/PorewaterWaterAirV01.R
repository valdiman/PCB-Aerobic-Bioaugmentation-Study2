
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
PwWaAirV01 = function(t, state, parms){
  
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
  Vpw <- 2.5 # cm3 porewater volume
  Va <- 125 # cm3 headspace volumne
  Aaw <- 20 # cm2 
  Aws <- 30 # cm2
  
  # Congener-specific constants
  Kaw <- 0.018048667 # PCB 19 dimensionless Henry's law constant @ 25 C
  dUaw <- 51590.22 # internal energy for the transfer of air-water for PCB 19 (J/mol)
  Kaw.t <- Kaw*exp(-dUaw / R * (1 / Tw.1 - 1 / Tst.1)) * Tw.1 / Tst.1
  Kow <- 10^(5.02) # PCB 19 octanol-water equilibrium partition coefficient (low value!!)
  dUow <-  -20988.94 # internal energy for the transfer of octanol-water for PCB 19 (J/mol)
  Kow.t <- Kow*exp(-dUow / R * (1 / Tw.1 -  1/ Tst.1))
  Koa <- 10^(7.765196745) # PCB 19 octanol-air equilibrium partition coefficient
  
  # Sediment partitioning
  M <- 0.1 # kg/L solid-water ratio
  foc <- 0.03 # organic carbon % in particles
  logKoc <- 0.94 * log10(Kow.t) + 0.42 # koc calculation
  Kd <- foc * 10^(logKoc) # L/kg sediment-water equilibrium partition coefficient
  
  # Air & water physical conditions
  D.water.air <- 0.2743615 # cm2/s water's diffusion coefficient in the gas phase @ Tair = 25 C, patm = 1013.25 mbars 
  D.co2.w <- 1.67606E-05 # cm2/s CO2's diffusion coefficient in water @ Tair = 25 C, patm = 1013.25 mbars 
  D.pcb.air <- D.water.air*(MW.pcb/MH2O)^(-0.5) # cm2/s PCB 19's diffusion coefficient in the gas phase (eq. 18-45)
  D.pcb.water <- D.co2.w*(MW.pcb/MCO2)^(-0.5) # cm2/s PCB 19's diffusion coefficient in water @ Tair = 25 C, patm = 1013.25 mbars
  v.H2O <- 0.010072884	# cm2/s kinematic viscosity of water @ Tair = 25
  V.water.air <- 0.003 # m/s water's velocity of air-side mass transfer without ventilation (eq. 20-15)
  V.co2.w <- 4.1*10^-2 # m/s mass transfer coefficient of CO2 in water side without ventilation
  SC.pcb.w <- v.H2O/D.pcb.water # Schmidt number PCB 19
  
  # Pore water MTC
  bl <- 0.21 # cm boundary layer thickness
  ks <- D.pcb.water * 60 * 60 * 24 / bl # [cm/d]
  ks.m.d <- ks / 100 # [m/d]
  
  # kaw calculations (air-water mass transfer coefficient)
  # i) Kaw.a, air-side mass transfer coefficient
  Kaw.a <- V.water.air*(D.pcb.air/D.water.air)^(0.67) # [m/s]
  # ii) Kaw.w, water-side mass transfer coefficient for PCB 19. 600 is the Schmidt number of CO2 at 298 K
  Kaw.w <- V.co2.w*(SC.pcb.w/600)^(-0.5) # [m/s]
  # iii) kaw, overall air-water mass transfer coefficient for PCB 19
  kaw.o <- (1/(Kaw.a*Kaw.t) + (1/Kaw.w))^-1 # [m/s]
  # iv) kaw, overall air-water mass transfer coefficient for PCB 19, units change
  kaw.o <- kaw.o * 100 * 60 * 60 * 24 # [cm/d]
  
  # Bioremediation rate
  kb <- parms$kb
  
  # derivatives dx/dt are computed below
  Cpw <- state[1]
  Cw <- state[2]
  Ca <- state[3]
  
  dCpwdt <- - kb * Cpw - ks * Aws / Vpw * (Cpw - Cw)
  dCwdt <- ks * Aws / Vw * (Cpw - Cw) -
    kaw.o * Aaw / Vw * (Cw - Ca / Kaw.t) -
    kb * Cw # [ng/L]
  dCadt <- kaw.o * Aaw / Va * (Cw - Ca / Kaw.t) # Ca = [ng/L]
    
  # The computed derivatives are returned as a list
  return(list(c(dCpwdt, dCwdt, dCadt)))
}

# Initial conditions and run function
{
  Ct <- 259.8342356 # ng/g PCB 19 sediment concentration
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
  Cpw <- Ct / Kd * 1000 # [ng/L]
}
cinit <- c(Cpw = Cpw, Cw = 0, Ca = 0) # [ng/L]
parms <- list(kb = 0.0) # Input
t <- seq(from = 0, to = 80, by = 1)
# Run the ODE function without specifying parms
out.1 <- ode(y = cinit, times = t, func = PwWaAirV01, parms = parms)
head(out.1)

{
  df.1 <- as.data.frame(out.1)
  colnames(df.1) <- c("time", "Cpw", "Cw", "Ca")
  Vpw <- 2.5 # [cm3]
  Vw <- 100 # [cm3]
  Va <- 125 # cm3
  df.1$Mpw <- df.1$Cpw * Vpw / 1000
  df.1$Mw <- df.1$Cw * Vw / 1000
  df.1$Ma <- df.1$Ca * Va / 1000
  df.1$Mt <- df.1$Cpw * Vpw / 1000 + df.1$Cw * Vw / 1000 + df.1$Ca * Va / 1000 # [ng]
  df.1$fpw <- (df.1$Cpw * Vpw / 1000) / df.1$Mt # [ng]
  df.1$fw <- (df.1$Cw * Vw / 1000) / df.1$Mt # [ng]
  df.1$fa <- (df.1$Ca * Va / 1000) / df.1$Mt # [ng]
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

