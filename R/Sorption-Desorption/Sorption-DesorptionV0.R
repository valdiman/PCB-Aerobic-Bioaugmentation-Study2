
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
  
  # Congener-specific constants
  Kow <- 10^(4.65) # PCB 4 octanol-water equilibrium partition coefficient
  dUow <-  -21338.96 # internal energy for the transfer of octanol-water for PCB 4 (J/mol)
  
  # Sediment partitioning
  foc <- 0.03 # organic carbon % in particles
  Kow.t <- Kow*exp(-dUow / R * (1 / Tw.1 -  1/ Tst.1))
  logKoc <- 0.94 * log10(Kow.t) + 0.42 # koc calculation
  K <- foc * 10^(logKoc) # L/kg sediment-water equilibrium partition coefficient
  
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
  
  # Determine the desorption rate based on time
  # Sorption start after day 2
  if (t <= 1) {  # If t is less than or equal to 1 day
    # Use fast desorption
    dCsdt <- - (f * kdf * Cs)
  } else {
    # Use slow desorption
    dCsdt <- - ((1 - f) * kds * Cs) + ka * Cw
  }
  
  # Derivative for water concentration
  dCwdt <- - ka * Cw * (t > 1) + (f * kdf * Cs * (t <= 1)) + ((1 - f) * kds * Cs * (t > 1))
  
  # The computed derivatives are returned as a list
  return(list(c(dCsdt, dCwdt)))
}

# Parameters and initial state
parms <- list(kdf = 1, kds = 0.001, f = 0.6, ka = 0.1, kb = 0.0)
cinit <- c(Cs = 63020.23, Cw = 0)
t.1 <- seq(from = 0, to = 75, by = 1)
# Run the ODE function without specifying parms
out.1 <- ode(y = cinit, times = t.1, func = rtm.PCB4, parms = parms)
head(out.1)

# Assuming out.1 is your data matrix
{
  df <- as.data.frame(out.1)
  colnames(df) <- c("time", "Cs", "Cw")
  df$Ct <- df$Cs + df$Cw
}

# Create the plot with all three lines
ggplot(data = df, aes(x = time)) +
  geom_line(aes(y = Cs, color = "Sediment (Cs)"), size = 1) +       # Line for Cs
  geom_line(aes(y = Cw, color = "Water (Cw)"), size = 1) +          # Line for Cw
  geom_line(aes(y = Ct, color = "Total (Ctotal)"), size = 1) +  # Line for total concentration Cs + Cw
  labs(title = "Concentrations vs Time", 
       x = "Time", 
       y = "Mass (ng)") +
  scale_color_manual(values = c("Sediment (Cs)" = "blue", "Water (Cw)" = "red",
                                "Total (Ctotal)" = "purple"),
                     name = "Concentrations") +
  theme_minimal()

