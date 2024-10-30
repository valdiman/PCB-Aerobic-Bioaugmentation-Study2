
# Reactive transport function ---------------------------------------------
rtm.PCB4 = function(t, state, parms){
  
  # Bioreactor parameters
  Vw <- 100 # cm3 water volume
  
  # Congener-specific constants
  Kow <- 10^(4.65) # PCB 4 octanol-water equilibrium partition coefficient
  
  # Sediment partitioning
  Cs <- 630.2023 * 1 # ng/g PCB 4 sediment concentration
  M <- 0.1 # kg/L solid-water ratio
  foc <- 0.03 # organic carbon % in particles
  logKoc <- 0.94 * log10(Kow) + 0.42 # koc calculation
  K <- foc * 10^(logKoc) # L/kg sediment-water equilibrium partition coefficient
  
  # SPME fiber constants
  Af <- 0.138 # cm2/cm SPME area
  Vf <- 0.000000069 # L/cm SPME volume/area
  L <- 1 # cm SPME length normalization to 1 cm
  Kf <- 10^(1.06 * log10(Kow) - 1.16) # PCB 4-SPME equilibrium partition coefficient [Lw/Lf]
  
  B <- (Vw + M * Vw * K + Vf * L * 1000) / Vw
  
  # Sorption and desorption rates
  kr <- parms$kr # 1/d
  f <- parms$f # fraction
  ka <- parms$ka # 1/d
  
  ko <- parms$ko # cm/d
  
  # derivatives dx/dt are computed below
  Cs <- state[1]
  Cw <- state[2]
  mf <- state[3]
  
  dCsdt <- (- Cs * kr * f + ka * Cw)
  dCwdt <- (- ka * Cw +  Cs * kr * f)
  dmfdt <- ko * Af * Vw / (Vf * 10000) * (Cw - mf / (Vf * Kf)) # Cw = [ng/L], mf = [ng/cmf]
  
  # The computed derivatives are returned as a list
  return(list(c(dCsdt, dCwdt, dmfdt)))
}
{
  Ct <- 630.2023 * 5# ng/g PCB 4 sediment concentration
  foc <- 0.03 # organic carbon % in sediment
  Kow <- 10^(4.65) # PCB 4 octanol-water equilibrium partition coefficient
  dUow <- -21338.96 # internal energy for the transfer of octanol-water for PCB 4 (J/mol)
  R <- 8.3144 # J/(mol K) molar gas constant
  Tst <- 25 #C air temperature
  Tst.1 <- 273.15 + Tst # air and standard temperature in K, 25 C
  Tw <- 20 # C water temperature
  Tw.1 <- 273.15 + Tw
  Kow.t <- Kow*exp(-dUow/R*(1/Tw.1-1/Tst.1))
  logKoc <- 0.94 * log10(Kow.t) + 0.42 # koc calculation
  K <- foc * 10^(logKoc) # L/kg sediment-water equilibrium partition coefficient
  M <- 0.1 # kg/L solid-water ratio
  Cw0 <- Ct * M * 1000 / (1 + M * K)
}

cinit <- c(Cs = 0, Cw = Cw0, mf = 0)
parms <- list(kr = 0.01, f = 0.6, ka = 0.03, ko = 10) # Input 
t.1 <- seq(from = 0, to = 75, by = 1)
# Run the ODE function without specifying parms
out.1 <- ode(y = cinit, times = t.1, func = rtm.PCB4, parms = parms)
head(out.1)

# Assuming out.1 is your data matrix
df <- as.data.frame(out.1)
colnames(df) <- c("time", "Cs", "Cw", "mf")

# Add a new column for the total concentration
df$Ctotal <- df$Cs + df$Cw

# Create the plot with all three lines
ggplot(data = df, aes(x = time)) +
  geom_line(aes(y = Cs, color = "Sediment (Cs)"), size = 1) +       # Line for Cs
  geom_line(aes(y = Cw, color = "Water (Cw)"), size = 1) +          # Line for Cw
  geom_line(aes(y = Ctotal, color = "Total (Ctotal)"), size = 1) +  # Line for total concentration Cs + Cw
  labs(title = "Concentrations vs Time", 
       x = "Time", 
       y = "Concentration (ng/L)") +
  scale_color_manual(values = c("Sediment (Cs)" = "blue", "Water (Cw)" = "red",
                                "Total (Ctotal)" = "purple"),
                     name = "Concentrations") +
  theme_minimal()

# Create the plot with all three lines
ggplot(data = df, aes(x = time, y = mf)) +
  geom_line(size = 1) +
  labs(title = "Concentrations vs Time", 
       x = "Time", 
       y = "SPME (ng/cm)") +
  theme_minimal()

