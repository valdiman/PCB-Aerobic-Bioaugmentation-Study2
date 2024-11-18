
# Packages and libraries --------------------------------------------------
# Install packages
install.packages("readxl")
install.packages("dplyr")
install.packages("ggplot2")

# Load libraries
{
  library(readxl) # read excel data
  library(dplyr) # organize data
  library(ggplot2) # plotting
}

# Experimental data -------------------------------------------------------
# Read data from Excel sheets
data <- read_excel("Data/SPMECalibration.xlsx", sheet = "data")

# Select individual congeners from datasets
pcb.ind <- "PCB52"

# Extract relevant columns from each dataset
pcbi <- data[, c("sample", "treatment", "replicate", "time", "length", pcb.ind)]

# SPME
Af <- 0.138 # [cm2/cm] SPME area
Vf <- 0.000000069 # [L/cm] SPME volume/area
L <- 1 # [cm] SPME length normalization to 1 cm
Kow <- 10^(5.84)
Kfi <- 10^(1.06 * log10(Kow) - 1.16) # [Lw/Lf]

# Calculate Cpw
Vw <- 30 # [ml or cm3]
ms <- 3 # [gr]
Ct <- 321.4900673 # [ng/g]
foc <- 0.03
logKoc <- 0.94 * log10(Kow) + 0.42
Kd <- foc * 10^(logKoc) # [L/kg]
Cpw <- Ct / Kd * 1000 # [ng/L]

# Calculate SPME volumes --------------------------------------------------
VSPME <- Vf * pcbi$length * 1000 # [mL]
# NOte: it should be < 0.01 mL. This should be fine.

# Depletion calculations --------------------------------------------------
Depl.1 <- Kfi * VSPME / (Kd * ms + Kfi * VSPME) * 100 # [%], no need to change units
# All experiments yield < 5%. This should be fine.

Depl.2 <- Kfi * (VSPME / Vw)
# Not sure about this
