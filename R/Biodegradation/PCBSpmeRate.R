# Biodegradation rate calaculations

# Packages and libraries --------------------------------------------------
# Install packages
install.packages("tidyverse")
install.packages("ggplot2")
install.packages("dplyr")

# Load libraries
library(tidyr)
library(dplyr)
library(ggplot2)

# Read data ---------------------------------------------------------------

df <- read.csv("Data/PCBBehaviorSPME.csv")

logKow <- as.numeric(df[13, 4:ncol(df)])
PCB_names <- colnames(df)[4:ncol(df)]

# Create logKow data frame
logKow_df <- data.frame(PCB = PCB_names, logKow = logKow)

df <- df[-13, ]


# Reshape the data to a long format for easier handling
df_long <- df %>%
  pivot_longer(cols = starts_with("PCB_"), 
               names_to = "PCB", 
               values_to = "Concentration")

# Join logKow values with df_long
df_long <- df_long %>%
  left_join(logKow_df, by = "PCB")

# Ensure time is numeric
df_long <- df_long %>%
  mutate(time = as.numeric(time))

df_long_filtered <- df_long %>%
  group_by(PCB) %>%
  filter(sum(Concentration == 0) < 3) %>%  # Keep PCBs with less than 3 zeros
  ungroup()  # Ungroup after filtering

# Fit the exponential decay model for PCB_1 to PCB_99 and extract k, M_w0, and R²
model_results <- df_long_filtered %>%
  filter(grepl("^PCB_[1-9][0-9]?$|^PCB_[1-9]$", PCB)) %>%  # Filter for PCB_1 to PCB_99
  group_by(PCB) %>%
  do({
    model <- nls(Concentration ~ M_w0 * exp(-k * time), 
                 data = ., 
                 start = list(M_w0 = max(.$Concentration, na.rm = TRUE), k = 0.1),
                 na.action = na.omit)  # Use na.omit to handle NAs
    
    # Extract coefficients
    k_value <- coef(model)["k"]
    M_w0 <- coef(model)["M_w0"]
    
    # Calculate R²
    residuals <- resid(model)
    total_variance <- sum((.$Concentration - mean(.$Concentration, na.rm = TRUE))^2)
    residual_variance <- sum(residuals^2)
    r_squared <- 1 - (residual_variance / total_variance)
    
    # Create a tibble with k, M_w0, and R²
    tibble(k = k_value, M_w0 = M_w0, R_squared = r_squared)
  }) %>%
  ungroup()





k_logKow_df <- model_results %>%
  left_join(df_long_filtered %>%
              select(PCB, logKow) %>%
              distinct(), by = "PCB")  # Ensure only unique PCBs are joined

k_model <- lm(k ~ logKow, data = k_logKow_df)

summary(k_model)

ggplot(k_logKow_df, aes(x = logKow, y = k)) +
  geom_point() +                      # Scatter plot of k vs logKow
  geom_smooth(method = "lm", se = TRUE, color = "blue") +  # Add regression line
  labs(title = "Regression of k against logKow",
       x = "logKow",
       y = "Decay Constant k") +
  theme_minimal()

