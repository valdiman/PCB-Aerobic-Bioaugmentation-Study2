# Script to analysis SPME data using ANOVA and Tukey's Honest test
# for lab data from bioaugmentation experiments

# Install Packages
install.packages("lattice")

# Load Libraries
library(lattice)

# Read PUF data -----------------------------------------------------------
PCB_data <- read.csv("Data/SPME_data_long.csv")

# See data
print(PCB_data)

# tPCB --------------------------------------------------------------------
# Plot data
xyplot(tPCB ~ factor(Time)|ID, group = Group,
       data = PCB_data, auto.key = TRUE)

xyplot(log10(tPCB) ~ factor(Time)|ID, group = Group,
       data = PCB_data, auto.key = TRUE)

# Stat Analysis (t test) --------------------------------------------------
# two-sample t test
# (1) AVL_S
t.test(log10(tPCB) ~ Group, data = subset(PCB_data, Time == 3 & ID == "AVL_S"))
t.test(log10(tPCB) ~ Group, data = subset(PCB_data, Time == 11 & ID == "AVL_S"))
t.test(log10(tPCB) ~ Group, data = subset(PCB_data, Time == 35 & ID == "AVL_S"))
# (2) AVL_NS
t.test(log10(tPCB) ~ Group, data = subset(PCB_data, Time == 16 & ID == "AVL_NS"))
t.test(log10(tPCB) ~ Group, data = subset(PCB_data, Time == 35 & ID == "AVL_NS"))
t.test(log10(tPCB) ~ Group, data = subset(PCB_data, Time == 75 & ID == "AVL_NS"))
# (3) NBH_NS
t.test(log10(tPCB) ~ Group, data = subset(PCB_data, Time == 16 & ID == "NBH_NS"))
t.test(log10(tPCB) ~ Group, data = subset(PCB_data, Time == 35 & ID == "NBH_NS"))
t.test(log10(tPCB) ~ Group, data = subset(PCB_data, Time == 75 & ID == "NBH_NS"))
# (4) AVL_S vs. AVL_NS
# Time 16 days
tmp.data = subset(PCB_data, Time == 16 & (ID == "AVL_S" | ID == "AVL_NS"))
# Plot data
xyplot(tPCB ~ factor(Time)|ID, group = Group, data = tmp.data,
       auto.key = TRUE)
xyplot(log10(tPCB) ~ factor(Time)|ID, group = Group, data = tmp.data,
       auto.key = TRUE)

# two-sample t test
t.test(log10(tPCB) ~ ID, data = subset(tmp.data, Group == "Control"))

# Time 35 days
tmp.data = subset(PCB_data, Time == 35 & (ID == "AVL_S" | ID == "AVL_NS"))
# Plot data
xyplot(tPCB ~ factor(Time)|ID, group = Group, data = tmp.data,
       auto.key = TRUE)
xyplot(log10(tPCB) ~ factor(Time)|ID, group = Group, data = tmp.data,
       auto.key = TRUE)

# two-sample t test
t.test(log10(tPCB) ~ ID, data = subset(tmp.data, Group == "Control"))
t.test(log10(tPCB) ~ ID, data = subset(tmp.data, Group == "Treatment"))

# PCB 4 -------------------------------------------------------------------
# Plot data
xyplot(PCB_4 ~ factor(Time)|ID, group = Group,
       data = PCB_data, auto.key = TRUE)

xyplot(log10(PCB_4) ~ factor(Time)|ID, group = Group,
       data = PCB_data, auto.key = TRUE)

# Stat Analysis (t test) --------------------------------------------------
# two-sample t test
# (1) AVL_S
t.test(log10(PCB_4) ~ Group, data = subset(PCB_data, Time == 3 & ID == "AVL_S"))
t.test(log10(PCB_4) ~ Group, data = subset(PCB_data, Time == 11 & ID == "AVL_S"))
t.test(log10(PCB_4) ~ Group, data = subset(PCB_data, Time == 35 & ID == "AVL_S"))
# (2) AVL_NS
t.test(log10(PCB_4) ~ Group, data = subset(PCB_data, Time == 16 & ID == "AVL_NS"))
t.test(log10(PCB_4) ~ Group, data = subset(PCB_data, Time == 35 & ID == "AVL_NS"))
t.test(log10(PCB_4) ~ Group, data = subset(PCB_data, Time == 75 & ID == "AVL_NS"))
# (3) NBH_NS
t.test(log10(PCB_4) ~ Group, data = subset(PCB_data, Time == 16 & ID == "NBH_NS"))
t.test(log10(PCB_4) ~ Group, data = subset(PCB_data, Time == 35 & ID == "NBH_NS"))
t.test(log10(PCB_4) ~ Group, data = subset(PCB_data, Time == 75 & ID == "NBH_NS"))

# (4) AVL_S vs. AVL_NS
# Time 16 days
tmp.data = subset(PCB_data, Time == 16 & (ID == "AVL_S" | ID == "AVL_NS"))
# Plot data
xyplot(PCB_4 ~ factor(Time)|ID, group = Group, data = tmp.data,
       auto.key = TRUE)
xyplot(log10(PCB_4) ~ factor(Time)|ID, group = Group, data = tmp.data,
       auto.key = TRUE)

# two-sample t test
t.test(log10(PCB_4) ~ ID, data = subset(tmp.data, Group == "Control"))

# Time 35 days
tmp.data = subset(PCB_data, Time == 35 & (ID == "AVL_S" | ID == "AVL_NS"))
# Plot data
xyplot(PCB_4 ~ factor(Time)|ID, group = Group, data = tmp.data,
       auto.key = TRUE)
xyplot(log10(PCB_4) ~ factor(Time)|ID, group = Group, data = tmp.data,
       auto.key = TRUE)

# two-sample t test
t.test(log10(PCB_4) ~ ID, data = subset(tmp.data, Group == "Control"))
t.test(log10(PCB_4) ~ ID, data = subset(tmp.data, Group == "Treatment"))

# PCB 19 -------------------------------------------------------------------
# Plot data
xyplot(PCB_19 ~ factor(Time)|ID, group = Group,
       data = PCB_data, auto.key = TRUE)

xyplot(log10(PCB_19) ~ factor(Time)|ID, group = Group,
       data = PCB_data, auto.key = TRUE)

# Stat Analysis (t test) --------------------------------------------------
# two-sample t test
# (1) AVL_S
t.test(log10(PCB_19) ~ Group, data = subset(PCB_data, Time == 3 & ID == "AVL_S"))
t.test(log10(PCB_19) ~ Group, data = subset(PCB_data, Time == 11 & ID == "AVL_S"))
t.test(log10(PCB_19) ~ Group, data = subset(PCB_data, Time == 35 & ID == "AVL_S"))
# (2) AVL_NS
t.test(log10(PCB_19) ~ Group, data = subset(PCB_data, Time == 16 & ID == "AVL_NS"))
t.test(log10(PCB_19) ~ Group, data = subset(PCB_data, Time == 35 & ID == "AVL_NS"))
t.test(log10(PCB_19) ~ Group, data = subset(PCB_data, Time == 75 & ID == "AVL_NS"))
# (3) NBH_NS
t.test(log10(PCB_19) ~ Group, data = subset(PCB_data, Time == 16 & ID == "NBH_NS"))
t.test(log10(PCB_19) ~ Group, data = subset(PCB_data, Time == 35 & ID == "NBH_NS"))
t.test(log10(PCB_19) ~ Group, data = subset(PCB_data, Time == 75 & ID == "NBH_NS"))

# (4) AVL_S vs. AVL_NS
# Time 16 days
tmp.data = subset(PCB_data, Time == 16 & (ID == "AVL_S" | ID == "AVL_NS"))
# Plot data
xyplot(PCB_19 ~ factor(Time)|ID, group = Group, data = tmp.data,
       auto.key = TRUE)
xyplot(log10(PCB_19) ~ factor(Time)|ID, group = Group, data = tmp.data,
       auto.key = TRUE)

# two-sample t test
t.test(log10(PCB_19) ~ ID, data = subset(tmp.data, Group == "Control"))

# Time 35 days
tmp.data = subset(PCB_data, Time == 35 & (ID == "AVL_S" | ID == "AVL_NS"))
# Plot data
xyplot(PCB_19 ~ factor(Time)|ID, group = Group, data = tmp.data,
       auto.key = TRUE)
xyplot(log10(PCB_19) ~ factor(Time)|ID, group = Group, data = tmp.data,
       auto.key = TRUE)

# two-sample t test
t.test(log10(PCB_19) ~ ID, data = subset(tmp.data, Group == "Control"))
t.test(log10(PCB_19) ~ ID, data = subset(tmp.data, Group == "Treatment"))





