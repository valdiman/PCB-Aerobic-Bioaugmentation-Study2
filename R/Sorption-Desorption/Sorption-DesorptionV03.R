
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
  Kaw <- 0.01344142 # PCB 4 dimensionless Henry's law constant @ 25 C
  dUaw <- 49662.48 # internal energy for the transfer of air-water for PCB 4 (J/mol)
  Kow <- 10^(4.65) # PCB 4 octanol-water equilibrium partition coefficient
  dUow <-  -21338.96 # internal energy for the transfer of octanol-water for PCB 4 (J/mol)
  Koa <- 10^(6.521554861) # PCB 4 octanol-air equilibrium partition coefficient
  
  # SPME fiber constants
  Af <- 0.138 # cm2/cm SPME area
  Vf <- 0.000000069 # L/cm SPME volume/area
  L <- 1 # cm SPME length normalization to 1 cm
  Kf <- 10^(1.06 * log10(Kow) - 1.16) # PCB 4-SPME equilibrium partition coefficient
  
  # Congener-specific constants
  Kow <- 10^(4.65) # PCB 4 octanol-water equilibrium partition coefficient
  dUow <-  -21338.96 # internal energy for the transfer of octanol-water for PCB 4 (J/mol)
  
  # Sediment partitioning
  foc <- 0.03 # organic carbon % in particles
  Kow.t <- Kow*exp(-dUow / R * (1 / Tw.1 -  1/ Tst.1))
  logKoc <- 0.94 * log10(Kow.t) + 0.42 # koc calculation
  K <- foc * 10^(logKoc) # L/kg sediment-water equilibrium partition coefficient
  
  Kaw.t <- 0.00939
  kaw.o <- 104
  
  # Sorption and desorption rates
  kdf <- parms$kdf # 1/d
  kds <- parms$kds # 1/d
  f <- parms$f # fraction
  kaf <- parms$kaf # 1/d
  kas <- parms$kas # 1/d
  ko <- parms$ko
  
  # Bioremediation rate
  kb <- parms$kb
  
  # derivatives dx/dt are computed below
  Cs <- state[1]
  Cw <- state[2]
  mf <- state[3]
  Ca <- state[4]
  
  # Define the maximum mf based on 5% of Cw mass
  mf_max <- 0.05 * Vf * L * Cw
  
  # Derivative for water concentration
  
  dCsdt <- - (f * kdf * Cs) - ((1 - f) * kds * Cs) + kaf * Cw + kas * Cw
  dCwdt <- - kaf * Cw - kas * Cw  + f * kdf * Cs + (1 - f) * kds * Cs - kb * Cw + kaw.o * Aaw / Vw * (Ca / (Kaw.t) - Cw)
  
  # Differential equation for dmfdt with saturation condition
  if (mf < mf_max) {
    dmfdt <- ko * Af * Vw / (Vf * L * 1e6) * (Cw - mf / (Vf * Kf))
  } else {
    dmfdt <- ko * Af * Vw / (Vf * L * 1e6) * (Cw - mf_max / (Vf * Kf))
  }
  
  dCadt <- kaw.o * Aaw / Va * (Cw - Ca / Kaw.t)
  
  # The computed derivatives are returned as a list
  return(list(c(dCsdt, dCwdt, dmfdt, dCadt)))
}

# Parameters and initial state
Ct <- 630.2023 # ng/g PCB 4 sediment concentration
M <- 0.1 # kg/L solid-water ratio
Cs0 <- Ct * M * 1000 # [ng/L]
cinit <- c(Cs = Cs0, Cw = 0, mf = 0, Ca = 0)
parms <- list(kdf = 0.1, kds = 0.01, f = 0.6, kaf = 0.1, kas = 0.01, kb = 0, ko = 1)
t.1 <- seq(from = 0, to = 40, by = 1)
# Run the ODE function without specifying parms
out.1 <- ode(y = cinit, times = t.1, func = rtm.PCB4, parms = parms,  maxsteps = 10000)
head(out.1)

Vw <- 100 # cm3 water volume
Va <- 125 # cm3 headspace volumne

# total mass calculations
{
  df <- as.data.frame(out.1)
  colnames(df) <- c("time", "Cs", "Cw", "Ca")
  df$Ct <- (df$Cs * Vw / 1000) + (df$Cw * Vw / 1000) + (df$Ca * Va / 1000)
}

# Fraction calculations
{
  df$Sf <- (df$Cs * Vw / 1000) / df$Ct
  df$Wf <- (df$Cw * Vw / 1000) / df$Ct
  df$Af <- (df$Ca * Va / 1000) / df$Ct
}

# Create the plot with all three lines
ggplot(data = df, aes(x = time)) +
  geom_line(aes(y = Cs * Vw / 1000, color = "Sediment (Cs)"), linewidth = 1) +       # Line for Cs
  geom_line(aes(y = Cw * Vw / 1000, color = "Water (Cw)"), linewidth = 1) +          # Line for Cw
  geom_line(aes(y = Ca * Va / 1000, color = "Air (Ca)"), linewidth = 1) +          # Line for Ca
  geom_line(aes(y = Ct, color = "Total (Ctotal)"), linewidth = 1) +  # Line for total concentration Cs + Cw
  labs(title = "Concentrations vs Time", 
       x = "Time", 
       y = "Mass (ng)") +
  scale_color_manual(values = c("Sediment (Cs)" = "blue", "Water (Cw)" = "red",
                                "Air (Ca)" = "orange", "Total (Ctotal)" = "purple"),
                     name = "Concentrations") +
  theme_minimal()

