
# Read data
# Organize data
# PCB4-RTMSControlV05 needs to be run before

Ct <- 630.2023 * 1  # ng/g PCB 4 sediment concentration
foc <- 0.03 # organic carbon % in sediment
Kow <- 10^(4.65) # PCB 4 octanol-water equilibrium partition coefficient
logKoc <- 0.94 * log10(Kow) + 0.42 # koc calculation
K <- foc * 10^(logKoc) # L/kg sediment-water equilibrium partition coefficient
ds <- 900 # g/L sediment density
M <- 0.1 # kg/L solid-water ratio
Cwi <- Ct * ds / (1 + M * K)
cinit <- c(Cw = Cwi, mf = 0, Ca = 0, mpuf = 0)
t.1 <- unique(pcb_combined_control$time)

# kb ----------------------------------------------------------------------
# Parameters that don't change
parms <- list(ka = 1, kd = 0.00000001)

# Sensitivity analysis: range of kb values
kb_values <- c(0, 0.02, 0.3, 0.5, 0.7)  # Varying kb

# Store results for each kb value
results_kb <- list()

# Run simulation for each kb value
for (kb in kb_values) {
  parms$kb <- kb
  result <- ode(y = cinit, times = t.1, func = rtm.PCB4, parms = parms)
  
  # Check if result is valid
  if (is.null(result) || nrow(result) == 0) {
    print(paste("No valid result for kb =", kb))
  } else {
    # Store the valid result in the results list
    results_kb[[as.character(kb)]] <- result
  }
}

# Extract results into a single data frame for plotting
all_results_mf_kb <- do.call(rbind, lapply(names(results_kb), function(kb) {
  result <- results_kb[[kb]]
  data.frame(time = result[, "time"], mf = result[, "mf"], kb = as.factor(kb))
}))

# Create a ggplot object for mf
p_mf_kb <- ggplot(all_results_mf_kb, aes(x = time, y = mf, color = kb)) +
  geom_line(linewidth = 1) +
  labs(x = "Time", y = "mf (ng/cm)") +
  theme_minimal() +
  theme(legend.position = "top")

# Determine the y-axis limits based on all kb values
mf_ranges_kb <- sapply(all_results_mf_kb$mf, range, na.rm = TRUE)  # Get range for all mf values
ylim_range_kb <- range(mf_ranges_kb, na.rm = TRUE)  # Overall range for ylim

# Recreate the plot with y-axis limits
p_mf_kb <- p_mf_kb + coord_cartesian(ylim = ylim_range_kb)

# Add a second plot for mpuf values (optional)
all_results_mpuf_kb <- do.call(rbind, lapply(names(results_kb), function(kb) {
  result <- results_kb[[kb]]
  data.frame(time = result[, "time"], mpuf = result[, "mpuf"], kb = as.factor(kb))
}))

# Create a ggplot object for mpuf
p_mpuf_kb <- ggplot(all_results_mpuf_kb, aes(x = time, y = mpuf, color = kb)) +
  geom_line(linewidth = 1, linetype = "dashed") +
  labs(x = "Time", y = "mpuf (ng/cm)") +
  theme_minimal() +
  theme(legend.position = "top")

# Combine both plots using grid.arrange
grid.arrange(p_mf_kb, p_mpuf_kb, ncol = 2)

# ka ----------------------------------------------------------------------
# Parameters that don't change
parms <- list(kd = 0.000687, kb = 0.02)  # Keep other parameters constant

# Sensitivity analysis: range of ka values
ka_values <- c(0.2, 2, 5, 10)  # Varying ka

# Store results for ka sensitivity analysis
results <- list()

# Loop through ka values and store results
for (ka in ka_values) {
  parms$ka <- ka
  result <- ode(y = cinit, times = t.1, func = rtm.PCB4, parms = parms)
  # Check if result is valid
  if (is.null(result) || nrow(result) == 0) {
    print(paste("No valid result for ka =", ka))
  } else {
    # Store the valid result in the results list
    results[[as.character(ka)]] <- result
  }
}

# Extract results into a single data frame for plotting
all_results_mf <- do.call(rbind, lapply(names(results), function(ka) {
  result <- results[[ka]]
  data.frame(time = result[, "time"], mf = result[, "mf"], ka = as.factor(ka))
}))

all_results_mpuf <- do.call(rbind, lapply(names(results), function(ka) {
  result <- results[[ka]]
  data.frame(time = result[, "time"], mpuf = result[, "mpuf"], ka = as.factor(ka))
}))

# Create a ggplot object for mf
p_mf <- ggplot(all_results_mf, aes(x = time, y = mf, color = ka)) +
  geom_line(linewidth = 1) +
  labs(x = "Time", y = "mf (ng/cm)") +
  theme_minimal() +
  theme(legend.position = "top")

# Create a ggplot object for mpuf
p_mpuf <- ggplot(all_results_mpuf, aes(x = time, y = mpuf, color = ka)) +
  geom_line(linewidth = 1, linetype = "dashed") +
  labs(x = "Time", y = "mpuf (ng/cm)") +
  theme_minimal() +
  theme(legend.position = "top")

# Combine both plots using grid.arrange
grid.arrange(p_mf, p_mpuf, ncol = 2)

# kd ----------------------------------------------------------------------
# Parameters that don't change
parms <- list(ka = 5, kb = 0.02)

# Sensitivity analysis: range of kd values
kd_values <- c(0.00000001, 0.0001, 0.001, 0.01, 0.02)  # Varying kd

# Store results for each kd value
results_kd <- list()

# Run simulation for each kd value
for (kd in kd_values) {
  parms$kd <- kd
  result <- ode(y = cinit, times = t.1, func = rtm.PCB4, parms = parms)
  
  # Check if result is valid
  if (is.null(result) || nrow(result) == 0) {
    print(paste("No valid result for kd =", kd))
  } else {
    # Store the valid result in the results list
    results_kd[[as.character(kd)]] <- result
  }
}

# Extract results into a single data frame for plotting
all_results_mf_kd <- do.call(rbind, lapply(names(results_kd), function(kd) {
  result <- results_kd[[kd]]
  data.frame(time = result[, "time"], mf = result[, "mf"], kd = as.factor(kd))
}))

# Create a ggplot object for mf
p_mf_kd <- ggplot(all_results_mf_kd, aes(x = time, y = mf, color = kd)) +
  geom_line(linewidth = 1) +
  labs(x = "Time", y = "mf (ng/cm)") +
  theme_minimal() +
  theme(legend.position = "top")

# Determine the y-axis limits based on all kd values
mf_ranges_kd <- sapply(all_results_mf_kd$mf, range, na.rm = TRUE)  # Get range for all mf values
ylim_range_kd <- range(mf_ranges_kd, na.rm = TRUE)  # Overall range for ylim

# Recreate the plot with y-axis limits
p_mf_kd <- p_mf_kd + coord_cartesian(ylim = ylim_range_kd)

# Add a second plot for mpuf values (optional)
all_results_mpuf_kd <- do.call(rbind, lapply(names(results_kd), function(kd) {
  result <- results_kd[[kd]]
  data.frame(time = result[, "time"], mpuf = result[, "mpuf"], kd = as.factor(kd))
}))

# Create a ggplot object for mpuf
p_mpuf_kd <- ggplot(all_results_mpuf_kd, aes(x = time, y = mpuf, color = kd)) +
  geom_line(linewidth = 1, linetype = "dashed") +
  labs(x = "Time", y = "mpuf (ng/cm)") +
  theme_minimal() +
  theme(legend.position = "top")

# Combine both plots using grid.arrange
grid.arrange(p_mf_kd, p_mpuf_kd, ncol = 2)
