# Biodegradation rate calaculations

# Packages and libraries --------------------------------------------------
# Install packages
install.packages("tidyverse")
install.packages("ggplot2")

# Load libraries
library(tidyr)
library(dplyr)
library(ggplot2)

# Read data ---------------------------------------------------------------

df <- read.csv("Data/PCBBehaviorSPME.csv")

# Reshape the data to a long format for easier handling
df_long <- df %>%
  pivot_longer(cols = starts_with("PCB_"), 
               names_to = "PCB", 
               values_to = "Concentration")

# Fit the exponential decay model for each PCB
model_results <- df_long %>%
  group_by(PCB) %>%
  do({
    model <- nls(Concentration ~ M_w0 * exp(-k * time), 
                 data = ., 
                 start = list(M_w0 = max(.$Concentration), k = 0.1))
    k_value <- coef(model)["k"]
    M_w0 <- coef(model)["M_w0"]
    
    tibble(k = k_value, M_w0 = M_w0)
  }) %>%
  ungroup()

print(model_results)
