
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
  
  # Pore water MTC
  bl <- 0.21 # cm boundary layer thickness
  ks <- 6.928 * 10^-6 * 60 * 60 * 24 / bl # [cm/d]
  ks.m.d <- ks / 100 # [m/d]
  
  kaw.o <- 132.23 # [cm/d]
  Kaw.t <- 0.012
  
  # Sediment desorption
  ksed <- parms$ksed # [1/d]
  
  # Diffusion term (new term for dCsdt)
  D <- parms$D # diffusion coefficient [cm^2/d], should be specified
  
  # Initial state variables
  Cs <- state[1]
  Cpw <- state[2]
  Cw <- state[3]
  Ca <- state[4]
  
  # Apply radial diffusion equation for Cs
  r <- 1  # Assuming a uniform radial distance, you can expand this if you have more spatial info
  dCsdr <-  (Cs  - 0) / r # derivative with respect to r (simplified as an example)
  d2Csdr2 <- (Cs  - 0) / r^2 # second derivative with respect to r (simplified as an example)
  dCsdt <- D * (1/r) * (r^2 * d2Csdr2)  # The simplified diffusion term
  
  # Porewater and other equations (no change in the existing equations)
  dCpwdt <- ksed * Apw / Vpw * (Cs - Cpw) -
    ks * Aws / Vpw * (Cpw - Cw)
  dCwdt <- ks * Aws / Vw * (Cpw - Cw) -
    kaw.o * Aaw / Vw * (Cw - Ca / Kaw.t)
  dCadt <- kaw.o * Aaw / Va * (Cw - Ca / Kaw.t)
  
  # Return the system of ODEs
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
parms <- list(ksed = 0.001, D = 0.1) # Diffusion coefficient 'D' needs to be defined
t <- seq(from = 0, to = 80, by = 1)
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

df.1$Ms <- df.1$Cs * ms / (ds * (1-n)) # Mass in sediment (ng)
df.1$Mpw <- df.1$Cpw * Vpw / 1000 # Mass in porewater (ng)
df.1$Mw <- df.1$Cw * Vw / 1000 # Mass in water (ng)
df.1$Ma <- df.1$Ca * Va / 1000 # Mass in air (ng)
# Total mass in the system (ng)
df.1$Mt <- df.1$Ms + df.1$Mpw + df.1$Mw + df.1$Ma

