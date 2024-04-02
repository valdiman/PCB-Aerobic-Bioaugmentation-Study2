# Script to analysis PUF data using ANOVA and Tukey's Honest test
# for lab data from bioaugmentation experiments
# Needs to change tPCB, LCPCB, PCB_4, PCB_19

# Install Packages
install.packages("stats")

# Load Libraries
library(stats)

# Read PUF data -----------------------------------------------------------
PCB_data <- read.csv("Data/SPME_data_long.csv")

# Display the new dataset
print(PCB_data)

# Prepare data for analysis -----------------------------------------------
# Subset the data per ID
AVL_S <- subset(PCB_data, ID == "AVL_S")
AVL_NS <- subset(PCB_data, ID == "AVL_NS")
NBH_NS <- subset(PCB_data, ID == "NBH_NS")

# (1) Perform Anova for AVL_S ---------------------------------------------
# Ensure 'Time' is a factor for ANOVA
AVL_S$Time <- as.factor(AVL_S$Time)

# Perform ANOVA
anova_result <- aov(PCB_4 ~ Time*Group, data = AVL_S)
summary(anova_result)

# Perform Tukey's Honest Significant Difference test
tukey_result <- TukeyHSD(anova_result, "Time:Group") # No difference found for tPCB!

# Access the relevant element of Tukey's HSD result
tukey_result_group <- tukey_result$"Time:Group"
# Remove rows with NA values
tukey_result_clean <- tukey_result_group[complete.cases(tukey_result_group), ]
# Subset to include only significant results
significant_results <- tukey_result_clean[tukey_result_clean[, 4] < 0.05, ]
print(significant_results)

# Export results (significant only)
write.csv(significant_results,
          file = "Output/Data/StatAnalysis/SPME/AVL_S_TukeyResultsPCB4.csv")

# (2) Perform Anova for AVL_NS --------------------------------------------
# Ensure 'Time' is a factor for ANOVA
AVL_NS$Time <- as.factor(AVL_NS$Time)

# Perform ANOVA
anova_result <- aov(PCB_4 ~ Time*Group, data = AVL_NS)
summary(anova_result) # No difference found for tPCB!

# Perform Tukey's Honest Significant Difference test
tukey_result <- TukeyHSD(anova_result, "Time:Group") # No difference found for tPCB!

# Access the relevant element of Tukey's HSD result
tukey_result_group <- tukey_result$"Time:Group"
# Remove rows with NA values
tukey_result_clean <- tukey_result_group[complete.cases(tukey_result_group), ]
# Subset to include only significant results
significant_results <- tukey_result_clean[tukey_result_clean[, 4] < 0.05, ]
print(significant_results)

# Export results (significant only)
write.csv(significant_results,
          file = "Output/Data/StatAnalysis/SPME/AVL_NS_TukeyResultsPCB4.csv")

# (3) Perform Anova for NBH_NS --------------------------------------------
# Ensure 'Time' is a factor for ANOVA
NBH_NS$Time <- as.factor(NBH_NS$Time)

# Perform ANOVA
anova_result <- aov(PCB_4 ~ Time*Group, data = NBH_NS)
summary(anova_result)

# Perform Tukey's Honest Significant Difference test
tukey_result <- TukeyHSD(anova_result, "Time:Group")

# Subset to include only significant results
significant_results <- tukey_result$`Time:Group`[tukey_result$`Time:Group`[, 4] < 0.05, ]
print(significant_results)

# Export results (significant only)
write.csv(significant_results, file = "Output/Data/StatAnalysis/SPME/NBH_NS_TukeyResultstPCB.csv")

# (4) Perform Anova for AVL_S and AVL_NS, only treatment ------------------
# (4.1) Analyze control
# Filter rows from AVL_S and AVL_SN where Group is "Control"
control_AVL_S <- AVL_S[AVL_S$Group == "Control", ]
control_AVL_NS <- AVL_NS[AVL_NS$Group == "Control", ]

# Filter rows with matching Time values
matching_times <- intersect(control_AVL_S$Time, control_AVL_NS$Time)

# Filter rows with matching Time values from AVL_S and AVL_SN
filtered_AVL_S <- control_AVL_S[control_AVL_S$Time %in% matching_times, ]
filtered_AVL_NS <- control_AVL_NS[control_AVL_NS$Time %in% matching_times, ]

# Concatenate the filtered data frames
AVL <- rbind(filtered_AVL_S, filtered_AVL_NS)

AVL$Time <- as.factor(AVL$Time)
AVL$Group <- factor(AVL$Group)

# Perform ANOVA
anova_result <- aov(tPCB ~ Time*ID, data = AVL)
summary(anova_result)

# Perform Tukey's Honest Significant Difference test
tukey_result <- TukeyHSD(anova_result, "Time:ID")

# Subset to include only significant results
significant_results <- tukey_result$`Time:ID`[tukey_result$`Time:ID`[, 4] < 0.05, ]
print(significant_results)

# Export results (significant only)
write.csv(significant_results, file = "Output/Data/StatAnalysis/SPME/AVL_Control_TukeyResultstPCB.csv")

# (4.2) Analyze treatment
# Filter rows from AVL_S and AVL_SN where Group is "Treatment"
treatment_AVL_S <- AVL_S[AVL_S$Group == "Treatment", ]
treatment_AVL_NS <- AVL_NS[AVL_NS$Group == "Treatment", ]

# Filter rows with matching Time values
matching_times <- intersect(treatment_AVL_S$Time, treatment_AVL_NS$Time)

# Filter rows with matching Time values from AVL_S and AVL_SN
filtered_AVL_S <- treatment_AVL_S[treatment_AVL_S$Time %in% matching_times, ]
filtered_AVL_NS <- treatment_AVL_NS[treatment_AVL_NS$Time %in% matching_times, ]

# Concatenate the filtered data frames
AVL <- rbind(filtered_AVL_S, filtered_AVL_NS)

AVL$Time <- as.factor(AVL$Time)
AVL$Group <- factor(AVL$Group)

# Perform t.test
ttest_result <- t.test(tPCB ~ ID, data = AVL)
print(ttest_result) # No difference found!

