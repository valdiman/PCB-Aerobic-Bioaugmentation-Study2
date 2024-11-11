
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

# Sorption Desorption V01 function ---------------------------------------------
SorptionDesorptionV01 = function(t, state, parms){
  
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
  
  # Congener-specific constants
  Kow <- 10^(4.65) # PCB 4 octanol-water equilibrium partition coefficient
  dUow <-  -21338.96 # internal energy for the transfer of octanol-water for PCB 4 (J/mol)
  
  # Sediment partitioning
  M <- 0.1 # kg/L solid-water ratio
  foc <- 0.03 # organic carbon % in particles
  Kow.t <- Kow*exp(-dUow / R * (1 / Tw.1 -  1/ Tst.1))
  logKoc <- 0.94 * log10(Kow) + 0.42 # koc calculation
  K <- foc * 10^(logKoc) # L/kg sediment-water equilibrium partition coefficient
  
  # Bioavailability factor B
  B <- (Vw + M * Vw * K) / Vw
  
  # Sorption and desorption rates
  kdf <- parms$kdf # 1/d
  kds <- parms$kds # 1/d
  f <- parms$f # fraction
  ka <- parms$ka # 1/d
  # Bioremediation rate
  kb <- parms$kb
  
  # derivatives dx/dt are computed below
  Cs <- state[1]
  Cw <- state[2]
  
  # Derivative for water concentration
  #dCsdt <- (- f * kdf * Cs * (t<= 1) - (1 - f) * kds * Cs * (t>1) + ka * Cw) / B
  #dCwdt <- (- ka * Cw + f * kdf * Cs * (t<= 1) + (1 - f) * kds * Cs * (t>1) - kb * Cw) / B
  
  dCsdt <- (- f * kdf * Cs - (1 - f) * kds * Cs + ka * Cw) / B
  dCwdt <- (- ka * Cw + f * kdf * Cs + (1 - f) * kds * Cs - kb * Cw) / B
  
  # The computed derivatives are returned as a list
  return(list(c(dCsdt, dCwdt)))
}

# Parameters and initial state
{
  # Estimating Cs0 (PCB 4 concentration in particles)
  Ct <- 630.2023 # ng/g PCB 4 sediment concentration
  M <- 0.1 # kg/L solid-water ratio
  Cs0 <- Ct * M * 1000 # [ng/L]
}
cinit <- c(Cs = 63020.23, Cw = 0) # [ng/L]
parms <- list(kdf = 5, kds = 0.1, f = 0.6, ka = 450, kb = 0) # Input
t <- seq(from = 0, to = 40, by = 1)
# Run the ODE function without specifying parms
out.1 <- ode(y = cinit, times = t, func = SorptionDesorptionV01, parms = parms)
head(out.1)

{
  df.1 <- as.data.frame(out.1)
  colnames(df.1) <- c("time", "Cs", "Cw")
  df.1$Ct <- df.1$Cs + df.1$Cw
  Vw <- 100 #[cm3]
  df.1$Mt <- (df.1$Cs + df.1$Cw) * Vw / 1000 # [ng]
  df.1$fP <- (df.1$Cs * Vw / 1000) / df.1$Mt # [ng]
  df.1$fd <- (df.1$Cw * Vw / 1000) / df.1$Mt # [ng]
}

# Create the plot with all three lines
ggplot(data = df.1, aes(x = time)) +
  geom_line(aes(y = Cs, color = "Sediment (Cs)"), linewidth = 1) +       # Line for Cs
  geom_line(aes(y = Cw, color = "Water (Cw)"), linewidth = 1) +          # Line for Cw
  geom_line(aes(y = Ct, color = "Total (Ctotal)"), linewidth = 1) +  # Line for total concentration Cs + Cw
  labs(title = "Concentrations vs Time", 
       x = "Time", 
       y = "Concentration (ng/L)") +
  scale_color_manual(values = c("Sediment (Cs)" = "blue", "Water (Cw)" = "red",
                                "Total (Ctotal)" = "purple"),
                     name = "Concentrations") +
  theme_minimal()
  ylim(0, 650)

# Sorption Desorption V02 function ----------------------------------------
# Water concentration and SPME
SorptionDesorptionV02 = function(t, state, parms){
  
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
  Kf <- 10^(1.06 * log10(Kow) - 1.16) # PCB 4-SPME equilibrium partition coefficient [Lw/Lf]
  
  # Sediment partitioning
  M <- 0.1 # kg/L solid-water ratio
  foc <- 0.03 # organic carbon % in particles
  Kow.t <- Kow*exp(-dUow / R * (1 / Tw.1 -  1/ Tst.1))
  logKoc <- 0.94 * log10(Kow) + 0.42 # koc calculation
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
  ro <- parms$ro # cm/d sampling rate for PUF
  ko <- parms$ko # cm/d mass transfer coefficient to SPME
  
  # derivatives dx/dt are computed below
  Cs <- state[1]
  Cw <- state[2]
  mf <- state[3]

  dCsdt <- (- f * kdf * Cs - (1 - f) * kds * Cs + ka * Cw) / B # [ng/L]
  dCwdt <- (- ka * Cw + f * kdf * Cs + (1 - f) * kds * Cs + kb * Cw -
              ko * Af / Vf /1000 * (Cw - mf / (Vf * Kf))) / B # [ng/L]
  dmfdt <- (ko * Af * Vw / (Vf * 1000 * 1000) * (Cw - mf / (Vf * Kf))) / B # Cw = [ng/L], mf = [ng/cmf]
  
  # The computed derivatives are returned as a list
  return(list(c(dCsdt, dCwdt, dmfdt)))
}

# Initial conditions and run function
{
  # Estimating Cs0 (PCB 4 concentration in particles)
  Ct <- 630.2023 # ng/g PCB 4 sediment concentration
  M <- 0.1 # kg/L solid-water ratio
  Cs0 <- Ct * M * 1000 # [ng/L]
}
cinit <- c(Cs = Cs0, Cw = 0, mf = 0)
parms <- list(ko = 1, kdf = 5, kds = 0.1, f = 0.6,
              ka = 450, kb = 0) # Input 
t <- seq(from = 0, to = 40, by = 1)
# Run the ODE function without specifying parms
out.2 <- ode(y = cinit, times = t, func = SorptionDesorptionV02, parms = parms)
head(out.2)

{
  df.2 <- as.data.frame(out.2)
  colnames(df.2) <- c("time", "Cs", "Cw", "mf")
  Vw <- 100 # [cm3]
  l <- 1 # [cm]
  df.2$Mp <- df.2$Cs * Vw / 1000
  df.2$Mw <- df.2$Cw * Vw / 1000
  df.2$Mf <- df.2$mf * l
  df.2$Mt <- (df.2$Cs + df.2$Cw) * Vw / 1000 + df.2$mf * l # [ng]
}

# Create the plot with all three lines
ggplot(data = df.2, aes(x = time)) +
  geom_line(aes(y = Cs, color = "Sediment (Cs)"), linewidth = 1) +       # Line for Cs
  geom_line(aes(y = Cw, color = "Water (Cw)"), linewidth = 1) +          # Line for Cw
  geom_line(aes(y = mf, color = "SPME (ng/cm)"), linewidth = 1) +  # Line for total concentration Cs + Cw
  labs(title = "Concentrations vs Time", 
       x = "Time", 
       y = "Concentration (ng/L)") +
  scale_color_manual(values = c("Sediment (Cs)" = "blue", "Water (Cw)" = "red",
                                "Total (Ctotal)" = "purple"),
                     name = "Concentrations") +
  theme_minimal() +
  ylim(0, 1)

# Sorption Desorption V05 function ----------------------------------------
# Passive samplers where included
  
SorptionDesorptionV05 = function(t, state, parms){
    
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
  logKoc <- 0.94 * log10(Kow) + 0.42 # koc calculation
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
  ro <- parms$ro # cm/d sampling rate for PUF
  ko <- parms$ko # cm/d mass transfer coefficient to SPME
  
  # derivatives dx/dt are computed below
  Cs <- state[1]
  Cw <- state[2]
  mf <- state[3]
  Ca <- state[4]
  mpuf <- state[5]
  
  dCsdt <- (- f * kdf * Cs - (1 - f) * kds * Cs + ka * Cw) / B # [ng/L]
  dCwdt <- (- ka * Cw + f * kdf * Cs + (1 - f) * kds * Cs + kb * Cw -
              ko * Af / Vf /1000 * (Cw - mf / (Vf * Kf))) / B # [ng/L]
  dmfdt <- (ko * Af * Vw / (Vf * 1000 * 1000) * (Cw - mf / (Vf * Kf))) / B # Cw = [ng/L], mf = [ng/cmf]
  
  dCadt <- (kaw.o * Aaw / Va * (Cw - Ca / Kaw.t) -
              ro * Apuf / 1000 * (Ca - mpuf / (Vpuf * d * Kpuf * 1000)))/ B # Ca = [ng/L]
  
  dmpufdt <- (ro * Apuf / 1000 * (Ca - mpuf / (Vpuf * d * Kpuf * 1000))) / B # Ca = [ng/L], mpuf = [ng]
  
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
  parms <- list(ro = 1, ko = 500, kdf = 5, kds = 0.1, f = 0.6,
                ka = 450, kb = 0) # Input 
  t <- seq(from = 0, to = 40, by = 1)
  # Run the ODE function without specifying parms
  out.5 <- ode(y = cinit, times = t, func = SorptionDesorptionV05, parms = parms)
  head(out.5)
  
  {
    df.5 <- as.data.frame(out.2)
    colnames(df.2) <- c("time", "Cs", "Cw", "mf", "Ca", "mpuf")
    Vw <- 100 # [cm3]
    Va <- 125 # [cm3]
    l <- 1 # [cm]
    df.2$Mp <- df.2$Cs * Vw / 1000
    df.2$Mw <- df.2$Cw * Vw / 1000
    df.2$Ma <- df.2$Ca * Va / 1000
    df.2$Mf <- df.2$mf * l
    df.2$Mpuf <- df.2$mpuf
    df.2$Mwf <- df.2$Cw * Vw / 1000 + df.2$mf * l
    df.2$Mapuf <- df.2$Ca * Va / 1000 + df.2$mpuf
    df.2$Mt <- (df.2$Cs + df.2$Cw) * Vw / 1000 + df.2$Ca * Va / 1000 + df.2$mf * l + df.2$mpuf # [ng]
  }

  1000/24/60
  
  