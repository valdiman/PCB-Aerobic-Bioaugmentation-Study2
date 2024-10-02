# AVL_S vs. AVL_NS

# Upload libraries
library(ggplot2)

# Read PUF data -----------------------------------------------------------
PCB_data <- read.csv("Data/PUF_data_long.csv")

# Include tPCB and PCB_ columns to be checked for zeros
PCB_columns <- c("tPCB", grep("^PCB_", names(PCB_data), value = TRUE))

# AVL_S vs. AVL_NS 16 days Analysis ---------------------------------------------------
# Filter the data for the specific Time and ID
subset_data.16 <- subset(PCB_data, Time == 16 & (ID == "AVL_S" | ID == "AVL_NS"))

# Select only ID = Control
subset_data.16.2 <- subset(subset_data.16, Group == "Control")

# Initialize a vector to store the p-values
p_values <- setNames(rep(NA, length(PCB_columns)), PCB_columns)

# Perform the t-test for each PCB variable, filtering out zeros in that specific variable, and store the p-values
for (pcb in PCB_columns) {
  # Check if there are any zeros in the specific PCB variable
  if (!any(subset_data.16.2[[pcb]] == 0)) {
    # Perform the t-test if no zeros are found and store the p-value
    t_test_result <- t.test(log10(subset_data.16.2[[pcb]]) ~ subset_data.16.2$ID)
    p_values[pcb] <- t_test_result$p.value
  }
}

# Convert the p-values vector to a data frame
p_values_df <- data.frame(PCB = names(p_values), p_value = p_values)
p_values_df <- na.omit(p_values_df)
p_values_df$PCB <- factor(p_values_df$PCB,
                          levels = unique(p_values_df$PCB))

# Create the plot using ggplot2
plot.16 <- ggplot(p_values_df, aes(x = PCB, y = -log10(p_value))) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") + # Adding a line at p = 0.05 for significance threshold
  labs(title = "Control: AVL_S vs. AVL_NS 16 d", x = "",
       y = "-log10 (p-value)") +
  theme_bw() +
  theme(aspect.ratio = 10/15,
        axis.text.x = element_text(face = "bold", size = 10, angle = 60, hjust = 1,
                                   color = "black"),
        axis.text.y = element_text(face = "bold", size = 10, color = "black"),
        axis.title.y = element_text(face = "bold", size = 12, color = "black"),
        axis.ticks = element_line(linewidth = 0.6, color = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position = "top")

# See plot
plot.16

# Save plot in folder
ggsave("Output/Plots/StatAnalysis/PUF/PlotPUFAVL_S_AVL_NS_Control_16d.png",
       plot = plot.16, width = 15, height = 10, dpi = 300)

# Filter the data for the specific Time and ID
subset_data.35 <- subset(PCB_data, Time == 35 & (ID == "AVL_S" | ID == "AVL_NS"))

# Select only ID = Control
subset_data.35.2 <- subset(subset_data.35, Group == "Control")

# Initialize a vector to store the p-values
p_values <- setNames(rep(NA, length(PCB_columns)), PCB_columns)

# Perform the t-test for each PCB variable, filtering out zeros in that specific variable, and store the p-values
for (pcb in PCB_columns) {
  # Check if there are any zeros in the specific PCB variable
  if (!any(subset_data.35.2[[pcb]] == 0)) {
    # Perform the t-test if no zeros are found and store the p-value
    t_test_result <- t.test(log10(subset_data.35.2[[pcb]]) ~ subset_data.35.2$ID)
    p_values[pcb] <- t_test_result$p.value
  }
}

# Convert the p-values vector to a data frame
p_values_df <- data.frame(PCB = names(p_values), p_value = p_values)
p_values_df <- na.omit(p_values_df)
p_values_df$PCB <- factor(p_values_df$PCB,
                          levels = unique(p_values_df$PCB))

# Create the plot using ggplot2
plot.35 <- ggplot(p_values_df, aes(x = PCB, y = -log10(p_value))) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") + # Adding a line at p = 0.05 for significance threshold
  labs(title = "Control: AVL_S vs. AVL_NS 35 d", x = "",
       y = "-log10 (p-value)") +
  theme_bw() +
  theme(aspect.ratio = 10/15,
        axis.text.x = element_text(face = "bold", size = 10, angle = 60, hjust = 1,
                                   color = "black"),
        axis.text.y = element_text(face = "bold", size = 10, color = "black"),
        axis.title.y = element_text(face = "bold", size = 12, color = "black"),
        axis.ticks = element_line(linewidth = 0.6, color = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position = "top")

# See plot
plot.35

# Save plot in folder
ggsave("Output/Plots/StatAnalysis/PUF/PlotPUFAVL_S_AVL_NS_Control_35d.png",
       plot = plot.35, width = 15, height = 10, dpi = 300)

# Select only ID = Treatment
subset_data.16.2 <- subset(subset_data.16, Group == "Treatment")

# Initialize a vector to store the p-values
p_values <- setNames(rep(NA, length(PCB_columns)), PCB_columns)

# Perform the t-test for each PCB variable, filtering out zeros in that specific variable, and store the p-values
for (pcb in PCB_columns) {
  # Check if there are any zeros in the specific PCB variable
  if (!any(subset_data.16.2[[pcb]] == 0)) {
    # Perform the t-test if no zeros are found and store the p-value
    t_test_result <- t.test(log10(subset_data.16.2[[pcb]]) ~ subset_data.16.2$ID)
    p_values[pcb] <- t_test_result$p.value
  }
}

# Convert the p-values vector to a data frame
p_values_df <- data.frame(PCB = names(p_values), p_value = p_values)
p_values_df <- na.omit(p_values_df)
p_values_df$PCB <- factor(p_values_df$PCB,
                          levels = unique(p_values_df$PCB))

# Create the plot using ggplot2
plot.16 <- ggplot(p_values_df, aes(x = PCB, y = -log10(p_value))) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") + # Adding a line at p = 0.05 for significance threshold
  labs(title = "Treatment: AVL_S vs. AVL_NS 16 d", x = "",
       y = "-log10 (p-value)") +
  theme_bw() +
  theme(aspect.ratio = 10/15,
        axis.text.x = element_text(face = "bold", size = 10, angle = 60, hjust = 1,
                                   color = "black"),
        axis.text.y = element_text(face = "bold", size = 10, color = "black"),
        axis.title.y = element_text(face = "bold", size = 12, color = "black"),
        axis.ticks = element_line(linewidth = 0.6, color = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position = "top")

# See plot
plot.16

# Save plot in folder
ggsave("Output/Plots/StatAnalysis/PUF/PlotPUFAVL_S_AVL_NS_Treatment_16d.png",
       plot = plot.16, width = 15, height = 10, dpi = 300)

# AVL_S vs. AVL_NS 35 days Analysis ---------------------------------------------------

# Filter the data for the specific Time and ID
subset_data.35 <- subset(PCB_data, Time == 35 & (ID == "NBH_NS" | ID == "AVL_NS"))

# Select only ID = Treatment
subset_data.35.2 <- subset(subset_data.35, Group == "Treatment")

# Initialize a vector to store the p-values
p_values <- setNames(rep(NA, length(PCB_columns)), PCB_columns)

# Perform the t-test for each PCB variable, filtering out zeros in that specific variable, and store the p-values
for (pcb in PCB_columns) {
  # Check if there are any zeros in the specific PCB variable
  if (!any(subset_data.35.2[[pcb]] == 0)) {
    # Perform the t-test if no zeros are found and store the p-value
    t_test_result <- t.test(log10(subset_data.35.2[[pcb]]) ~ subset_data.35.2$ID)
    p_values[pcb] <- t_test_result$p.value
  }
}

# Convert the p-values vector to a data frame
p_values_df <- data.frame(PCB = names(p_values), p_value = p_values)
p_values_df <- na.omit(p_values_df)
p_values_df$PCB <- factor(p_values_df$PCB,
                          levels = unique(p_values_df$PCB))

# Create the plot using ggplot2
plot.35 <- ggplot(p_values_df, aes(x = PCB, y = -log10(p_value))) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") + # Adding a line at p = 0.05 for significance threshold
  labs(title = "Treatment: AVL_S vs. AVL_NS 35 d", x = "",
       y = "-log10 (p-value)") +
  theme_bw() +
  theme(aspect.ratio = 10/15,
        axis.text.x = element_text(face = "bold", size = 10, angle = 60, hjust = 1,
                                   color = "black"),
        axis.text.y = element_text(face = "bold", size = 10, color = "black"),
        axis.title.y = element_text(face = "bold", size = 12, color = "black"),
        axis.ticks = element_line(linewidth = 0.6, color = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position = "top")

# See plot
plot.35

# Save plot in folder
ggsave("Output/Plots/StatAnalysis/PUF/PlotPUFAVL_S_AVL_NS_Treatment_35d.png",
       plot = plot.35, width = 15, height = 10, dpi = 300)

