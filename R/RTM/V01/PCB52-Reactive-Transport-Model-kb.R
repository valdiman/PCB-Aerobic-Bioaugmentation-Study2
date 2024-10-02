# Code to model PCB 52 in laboratory experiments
# using sediment from Altavista, VI. Passive measurements
# of PCB 52 in the water and the air phases are predicted and
# linked to the water and air concentrations from the passive
# samplers.

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

# Read data ---------------------------------------------------------------
exp.data.SPME <- read.csv("Data/SPME_data_long.csv")
exp.data.PUF <- read.csv("Data/PUF_data_long.csv")

# Add sampler
exp.data.SPME$sampler <- "SPME"
exp.data.PUF$sampler <- "PUF"

# Combine
exp.data <- rbind(exp.data.PUF, exp.data.SPME)

# Select individual congener from datasets
pcb.ind <- "PCB_52"

# Extract relevant columns from each dataset
pcbi <- exp.data[, c("ID", "Group", "time", "sampler", pcb.ind)]

# Organize data -----------------------------------------------------------
# Remove lost sample(s), NA
pcbi <- pcbi[!is.na(pcbi$PCB_52), ]

# SPME = SPME fiber sampler
# PUF = PUF sampler

# Pull congener-specific data from the dataset without averaging
pcbi.spme.control <- pcbi %>%
  filter(ID == "AVL_NS", Group == "Control", sampler == "SPME") %>%
  rename("mf_Control" = PCB_4)

# Need to add an extra point to 16 days with value of NA
pcbi.spme.treatment <- pcbi %>%
  filter(ID == "AVL_NS", Group == "Treatment", sampler == "SPME") %>%
  rename("mf_Treatment" = PCB_4) %>%
  add_row(time = 16, mf_Treatment = NA, ID = "AVL_NS",
          Group = "Treatment", sampler = "SPME")

pcbi.puf.control <- pcbi %>%
  filter(ID == "AVL_NS", Group == "Control", sampler == "PUF") %>%
  rename("mpuf_Control" = PCB_4)

pcbi.puf.treatment <- pcbi %>%
  filter(ID == "AVL_NS", Group == "Treatment", sampler == "PUF") %>%
  rename("mpuf_Treatment" = PCB_4)

# Combine the mf and mpuf data for Control
pcb_combined_control <- cbind(
  pcbi.spme.control %>%
    select(time, mf_Control),
  pcbi.puf.control %>%
    select(mpuf_Control)
)

# Add a row for time = 0 if needed
pcb_combined_control <- rbind(
  data.frame(time = 0, mf_Control = 0, mpuf_Control = 0),
  pcb_combined_control
)

# Combine the mf and mpuf data for Treatment
pcb_combined_treatment <- cbind(
  pcbi.spme.treatment %>%
    select(time, mf_Treatment),
  pcbi.puf.treatment %>%
    select(mpuf_Treatment)
)

# Add a row for time = 0 if needed
pcb_combined_treatment <- rbind(
  data.frame(time = 0, mf_Treatment = 0, mpuf_Treatment = 0),
  pcb_combined_treatment
)

# Reactive transport function ---------------------------------------------

rtm.PCB52 = function(t, c, parms){
  
  # Experimental conditions
  MH2O <- 18.0152 # g/mol water molecular weight
  Mco2 <- 44.0094 # g/mol CO2 molecular weight
  MW.pcb <- 291.976 # g/mol PCB 52 molecular weight
  
  # Bioreactor parameters
  Vw <- 100 # cm3 water volume
  Va <- 125 # cm3 headspace volumne
  Aaw <- 20 # cm2 
  Aws <- 30 # cm2
  
  # Congener-specific constants
  Kaw <-  0.0130452 # PCB52 dimensionless Henry's law constant @ 25 C
  Kow <- 10^(5.84) # PCB52 octanol-water equilibrium partition coefficient
  Koa <-  10^(8.351339075) # PCB52 octanol-air equilibrium partition coefficient
  
  # PUF constants 
  Vpuf <- 0.000029 # m3 volume of PUF
  Kpuf <- 10^(0.6366 * log10(Koa) - 3.1774)# m3/g PCB 52-PUF equilibrium partition coefficient
  d <- 0.0213*100^3 # g/m3 density of PUF
  ro <- 0.0005 # m3/d sampling rate
  
  # SPME fiber constants
  Af <- 0.138 # cm2/cm SPME area
  Vf <- 0.000000069 # L/cm SPME volume/area
  L <- 1 # cm SPME length average
  Kf <- 10^(1.06 * log10(Kow) - 1.16) # PCB 52-SPME equilibrium partition coefficient
  ko <- 70 # cm/d PCB 52 mass transfer coefficient to SPME
  
  # Air & water physical conditions
  D.water.air <- 0.2743615 # cm2/s water's diffusion coefficient in the gas phase @ Tair = 25 C, patm = 1013.25 mbars 
  D.co2.w <- 1.67606E-05 # cm2/s CO2's diffusion coefficient in water @ Tair = 25 C, patm = 1013.25 mbars 
  D.pcb.air <- D.water.air*(MW.pcb/MH2O)^(-0.5) # cm2/s PCB 4's diffusion coefficient in the gas phase (eq. 18-45)
  D.pcb.water <- D.co2.w*(MW.pcb/MCO2)^(-0.5) # cm2/s PCB 4's diffusion coefficient in water @ Tair = 25 C, patm = 1013.25 mbars
  v.H2O <- 0.010072884	# cm2/s kinematic viscosity of water @ Tair = 25
  V.water.air <- 0.003 # m/s water's velocity of air-side mass transfer without ventilation (eq. 20-15)
  V.co2.w <- 4.1*10^-2 # m/s mass transfer coefficient of CO2 in water side without ventilation
  SC.pcb.w <- v.H2O/D.pcb.water # Schmidt number PCB 4
  bl <- 0.21 # cm boundary layer thickness
  ks <- ks <- D.pcb.water / bl # [cm/s]
  ks.m.d <- ks * 60 * 60 * 24 / 100 # [m/d]
  
  # kaw calculations (air-water mass transfer coefficient)
  # i) Kaw.a, air-side mass transfer coefficient
  Kaw.a <- V.water.air*(D.pcb.air/D.water.air)^(0.67) # [m/s]
  # ii) Kaw.w, water-side mass transfer coefficient for PCB 4. 600 is the Schmidt number of CO2 at 298 K
  Kaw.w <- V.co2.w*(SC.pcb.w/600)^(-0.5) # [m/s]
  # iii) kaw, overall air-water mass transfer coefficient for PCB 4
  kaw.o <- (1/(Kaw.a*Kaw) + (1/Kaw.w))^-1 # [m/s]
  # v) kaw, overall air-water mass transfer coefficient for PCB 4, units change
  kaw.o <- kaw.o*100*60*60*24 # [cm/d]
  
  # Biotransformation rate
  kb <- 0.06 # 1/d, value changes depending on experiment, i.e., control = 0, treatments LB400 = 0.130728499

  # derivatives dx/dt are computed below
  Cpw <- state[1]
  Cw <- state[2]
  mf <- state[3]
  Ca <- state[4]
  mpuf <- state[5]
  
  dCpwdt <- -kb * Cpw
  dCwdt <- kaw.o * Aaw / Vw * (Ca / (Kaw) - Cw) + ks * Aws * 60 * 60 * 24 / Vw * (Cpw - Cw) - kb * Cw # 864 to change second to days and um to m, Ca in [ng/L]
  dmfdt <- ko * Af /(L * 1000) * (Cw - mf / (Vf * L * Kf)) # Cw = [ng/L], mf = [ng/cmf]
  dCadt <- kaw.o * Aaw / Va * (Cw - Ca / Kaw)
  dpufdt <- ro * Ca * 1000 - ro * (mpuf / (Vpuf * d)) / (Kpuf) # Ca = [ng/L], mpuf = [ng]
  
  # The computed derivatives are returned as a list
  return(list(c(dCpwdt, dCwdt, dmfdt, dCadt, dpufdt)))
}

# Initial conditions and run function
# Estimating Cpw (PCB 52 concentration in sediment porewater)
Ct <- 321.4900673 # ng/g PCB 52 sediment concentration
foc <- 0.03 # organic carbon % in sediment
Kow <- 10^(5.84) # PCB52 octanol-water equilibrium partition coefficient
logKoc <- 0.94 * log10(Kow) + 0.42 # koc calculation
Kd <- foc * 10^(logKoc) # L/kg sediment-water equilibrium partition coefficient
Cpw <- Ct / Kd * 1000 # [ng/L]

cinit <- c(Cpw = Cpw, Cw = 0, mf = 0, Ca = 0, mpuf = 0)
t.1 <- unique(pcb_combined_control$time)
# Run the ODE function without specifying parms
out.1 <- ode(y = cinit, times = t.1, func = rtm.PCB4)
head(out.1)



# Plotting ----------------------------------------------------------------
# Plot of predicted vs experimental LB400 data
names(out.plot.kb) <- c("time", "PCB 52 predicted Cw", "PCB 52 predicted mass fiber",
                     "PCB 52 predicted Ca", "PCB 52 predicted mass PUF")
# Predicted
# spme
out.plot.mspme <- subset(out.plot.kb,
                      select = c(time, `PCB 52 predicted mass fiber`))
# puf
out.plot.mpuf <- subset(out.plot.kb,
                        select = c(time, `PCB 52 predicted mass PUF`))
pred.mspme.kb <- melt(out.plot.mspme, id.var = c("time"),
                variable.name = "compartment", value.name = "mass")
pred.mpuf.kb <- melt(out.plot.mpuf, id.var = c("time"),
                  variable.name = "compartment", value.name = "mass")
# Experimental
# spme
pcb.52.mspme <- subset(pcb.52,
                    select = c(time,`PCB 52 mass in SPME (ng/cm); LB400`))
# puf
pcb.52.mpuf <- subset(pcb.52,
                      select = c(time, `PCB 52 mass in PUF (ng); LB400`))
exp.mspme <- melt(pcb.52.mspme, id.var = c("time"), variable.name = "compartment",
               value.name = "mass")
exp.mpuf <- melt(pcb.52.mpuf, id.var = c("time"), variable.name = "compartment",
                 value.name = "mass")  
names(pcb.52.mspme) <- c("time", "PCB 52 measured mass fiber")
names(pcb.52.mpuf) <- c("time", "PCB 52 measured mass PUF")

mspme.kb <- ggplot(data = pred.mspme.kb, aes(x = time, y = mass)) +
  geom_line(colour = 2) +
  geom_point(data = exp.mspme, aes(x = time, y = mass), color = 2) +
  theme_bw() +
  theme(aspect.ratio = 6/10) + #Change the plot title for the SPME fiber plot here!
  ggtitle("Treatment with LB400 - SPME Fiber Results for PCB 52") +
  xlab(expression(bold("Time (day)"))) +
  ylab(expression(bold("PCB 52 fiber mass accumulated (ng/cm)"))) +
  ylim(0, 0.25) +
  xlim(0, 80)

print(mspme.kb)

mpuf.kb <- ggplot(data = pred.mpuf.kb, aes(x = time, y = mass)) +
  geom_line(colour = 4) +
  geom_point(data = exp.mpuf, aes(x = time, y = mass), color = 4) +
  theme_bw() +
  theme(aspect.ratio = 6/10) + #And the PUF plot here!
  ggtitle("Treatment with LB400 - PUF Results for PCB 52") +
  xlab(expression(bold("Time (day)"))) +
  ylab(expression(bold("PCB 52 PUF mass accumulated (ng/PUF)"))) +
  ylim(0, 40) +
  xlim(0, 80)

print(mpuf.kb)
