
# Calculate SPME volume per length ----------------------------------------
# Parameters
D_inner_um <- 110  # Inner diameter in um
thickness_um <- 7  # Thickness in um

# Convert diameters and thickness from micrometers to centimeters
D_inner_cm <- D_inner_um / 10000  # 1 cm = 10000 um
thickness_cm <- thickness_um / 10000  # 1 cm = 10000 um

# Calculate the outer diameter
D_outer_cm <- D_inner_cm + 2 * thickness_cm

# Calculate the volume of the outer cylinder per cm of height (in cm^3)
V_outer_per_cm <- pi * (D_outer_cm / 2)^2  # Volume per cm of height

# Calculate the volume of the inner cylinder per cm of height (in cm^3)
V_inner_per_cm <- pi * (D_inner_cm / 2)^2  # Volume per cm of height

# Calculate the volume of the wall per cm of height (in cm^3)
V_wall_per_cm <- V_outer_per_cm - V_inner_per_cm

# Convert to liters (1 cm^3 = 0.001 L)
V_wall_per_cm_liters <- V_wall_per_cm * 0.001

# Display results
V_wall_per_cm_liters # [L/cm]

# Calculate SPME are per length -------------------------------------------
# Parameters
D_inner_um <- 110  # Inner diameter in micrometers
thickness_um <- 7  # Thickness of the wall in micrometers

# Convert to centimeters
D_inner_cm <- D_inner_um * 1e-4  # Convert inner diameter to cm
D_outer_cm <- (D_inner_um + 2 * thickness_um) * 1e-4  # Convert outer diameter to cm

# Calculate the surface area of the wall per cm of length
A_wall <- pi * (D_outer_cm + D_inner_cm)  # Area in cm^2

# Display result
A_wall # [cm2/cm]
