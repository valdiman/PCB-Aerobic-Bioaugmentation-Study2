
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
  
  K <- 2128 # L/kg
  M <- 0.1 # kg/L
  
  # Bioavailability factor B
  Vw <- 10 # cm3
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
parms <- list(kdf = 500, kds = 0.1, f = 0.6, ka = 200, kb = 0) # Input 
cinit <- c(Cs = 63020.23, Cw = 0)
t.1 <- seq(from = 0, to = 50, by = 1)
# Run the ODE function without specifying parms
out.1 <- ode(y = cinit, times = t.1, func = rtm.PCB4, parms = parms)
head(out.1)

{
  df <- as.data.frame(out.1)
  colnames(df) <- c("time", "Cs", "Cw")
  df$Ct <- df$Cs + df$Cw
  Vw <- 10 #[cm3]
  df$Mt <- (df$Cs + df$Cw) * Vw / 1000 # [ng]
  df$fP <- (df$Cs * Vw / 1000) / df$Mt # [ng]
  df$fd <- (df$Cw * Vw / 1000) / df$Mt # [ng]
}

# Create the plot with all three lines
ggplot(data = df, aes(x = time)) +
  geom_line(aes(y = Cs, color = "Sediment (Cs)"), linewidth = 1) +       # Line for Cs
  geom_line(aes(y = Cw, color = "Water (Cw)"), linewidth = 1) +          # Line for Cw
  geom_line(aes(y = Ct, color = "Total (Ctotal)"), linewidth = 1) +  # Line for total concentration Cs + Cw
  labs(title = "Concentrations vs Time", 
       x = "Time", 
       y = "Mass (ng)") +
  scale_color_manual(values = c("Sediment (Cs)" = "blue", "Water (Cw)" = "red",
                                "Total (Ctotal)" = "purple"),
                     name = "Concentrations") +
  theme_minimal()
  ylim(10000, 65000)

