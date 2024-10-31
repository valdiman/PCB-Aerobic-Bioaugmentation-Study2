
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
  
  # Sorption and desorption rates
  kdf <- parms$kdf # 1/d
  kds <- parms$kds # 1/d
  f <- parms$f # fraction
  ka <- parms$ka # 1/d
  
  # derivatives dx/dt are computed below
  Cs <- state[1]
  Cw <- state[2]
  
  # Determine the desorption rate based on time
  if (t <= 1) {  # If t is less than or equal to 5 day
    # Use fast desorption
    dCsdt <- - (f * kdf * Cs) + ka * Cw
  } else {
    # Use slow desorption
    dCsdt <- - ((1 - f) * kds * Cs) + ka * Cw
  }
  
  # Derivative for water concentration
  dCwdt <- - ka * Cw + (f * kdf * Cs * (t <= 1)) + ((1 - f) * kds * Cs * (t > 1))
  
  # The computed derivatives are returned as a list
  return(list(c(dCsdt, dCwdt)))
}

cinit <- c(Cs = 63020.23, Cw = 0)
parms <- list(kdf = 0.5, kds = 0.000001, f = 0.6, ka = 0.03) # Input 
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

