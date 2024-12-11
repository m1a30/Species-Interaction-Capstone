
## Visualizations ####
library(qgraph)

# Blanton Ridge #######

# Reading in the data ####
mean_interactions_br_df <- read.csv("mean_interaction_matrix_br.csv")

rownames(mean_interactions_br_df) <- mean_interactions_br_df$RowNames
mean_interactions_br_df <- mean_interactions_br_df[, -which(colnames(mean_interactions_br_df) == "RowNames")]

all_species_br <- unique(c(rownames(mean_interactions_br_df), colnames(mean_interactions_br_df)))

# based on Bimler figure code!! #######
# Convert to matrix
mean_interactions_br_matrix <- as.matrix(mean_interactions_br_df)

# Ensure matrix is correctly oriented
mean_interactions_br_matrix <- t(mean_interactions_br_matrix)  # Arrows point towards focal species
# rownames gives the focal names and colnames gives all species names
all_species <- unique(c(rownames(mean_interactions_br_df), colnames(mean_interactions_br_df)))

# Create a square matrix
square_matrix <- matrix(0, nrow = length(all_species), ncol = length(all_species))
rownames(square_matrix) <- all_species
colnames(square_matrix) <- all_species

# Fill the square matrix with values from mean_interactions_br_matrix
for (i in seq_len(nrow(mean_interactions_br_df))) {
  for (j in seq_len(ncol(mean_interactions_br_df))) {
    square_matrix[rownames(mean_interactions_br_df)[i], colnames(mean_interactions_br_df)[j]] <- mean_interactions_br_df[i, j]
  }
}

mean_interactions_br_matrix <- square_matrix



# Define species groups dynamically
foundation <- c('CLAPUR', 'COLLOM', 'COLLIN', "GILCAP", "NAVSQU", 'EPIDEN', "PLAFIG", "PLECON")  # Replace with foundation species at Blanton Ridge
all_species <- rownames(mean_interactions_br_matrix)  # All species in the dataset
invasives_br <- setdiff(all_species, foundation)         # Invasives = all species except foundation

# Set up colors for nodes
names(all.sp) <- all_species

all.sp[invasives_br] <- 'forestgreen'  # Color for invasive species
all.sp[foundation] <- 'purple'    # Color for foundation species

# Plot competition and facilitation separately #####
# Competition only
qgraph(mean_interactions_br_matrix,
       layout = 'circle',
       posCol = rgb(red = 0, green = 0, blue = 0, alpha = 0),  # Facilitation = transparent
       negCol = 'orange',                                     # Competition = red
       color = all.sp,
       labels = rownames(mean_interactions_br_matrix), label.color = "white", label.cex = 1.25,
       fade = TRUE, directed = TRUE, curve = 2,
       title = 'Competition at Blanton Ridge', title.cex = 1)

# Facilitation only
qgraph(mean_interactions_br_matrix,
       layout = 'circle',
       posCol = 'royalblue4',                                 # Facilitation = blue
       negCol = rgb(red = 0, green = 0, blue = 0, alpha = 0), # Competition = transparent
       color = all.sp,
       labels = rownames(mean_interactions_br_matrix),  label.color = "white", label.cex = 1.25,
       fade = TRUE, directed = TRUE,
       title = 'Facilitation at Blanton Ridge', title.cex = 1)



# both facil and comp together
qgraph(mean_interactions_br_matrix,
       layout = 'circle',           # Circular layout
       posCol = 'royalblue4',       # Facilitation = blue
       negCol = 'orange',           # Competition = orange
       color = all.sp,              # Node colors
       labels = rownames(mean_interactions_br_matrix),  # Node labels
       fade = TRUE, curve = 2,   label.color = "white", label.cex = 1.25,
       title = 'Competition and Facilitation at Blanton Ridge', title.cex = 1)


# just looking at our focals from larger plot #####
# Define foundation species
foundation <- c('CLAPUR', 'COLLOM', 'COLLIN', "GILCAP", "NAVSQU", 'EPIDEN', "PLAFIG", "PLECON")


# Extract the submatrix: interactions toward foundation species
# Extract submatrix for foundation-to-foundation interactions
native_interactions <- mean_interactions_br_matrix[foundation, foundation]

# Plot foundation-to-foundation interactions
qgraph(native_interactions,
       layout = 'circle',
       negCol = 'orange',  # competition = orange
       posCol = 'royalblue4',      # facilitation = green
       color = 'purple',       # All nodes are foundations
       labels = foundation,
       fade = TRUE, directed = TRUE,   label.color = "white", label.cex = 1.25,
       title = 'Interactions Among Natives at BR', title.cex = 1)



# South Eugene Meadows #######
## Reading in the data ####
mean_interactions_sem_df <- read.csv("mean_interaction_matrix_sem.csv")

rownames(mean_interactions_sem_df) <- mean_interactions_sem_df$RowNames
mean_interactions_sem_df <- mean_interactions_sem_df[, -which(colnames(mean_interactions_sem_df) == "RowNames")]

all_species_sem <- unique(c(rownames(mean_interactions_sem_df), colnames(mean_interactions_sem_df)))

# based on Bimler figure code!! #######
# Convert to matrix
mean_interactions_sem_matrix <- as.matrix(mean_interactions_sem_df)

# Ensure matrix is correctly oriented
mean_interactions_sem_matrix <- t(mean_interactions_sem_matrix)  # Arrows point towards focal species

# Create a square matrix
square_matrix_sem <- matrix(0, nrow = length(all_species_sem), ncol = length(all_species_sem))
rownames(square_matrix_sem) <- all_species_sem
colnames(square_matrix_sem) <- all_species_sem

# Fill the square matrix with values from mean_interactions_sem_matrix
for (i in seq_len(nrow(mean_interactions_sem_df))) {
  for (j in seq_len(ncol(mean_interactions_sem_df))) {
    square_matrix_sem[rownames(mean_interactions_sem_df)[i], colnames(mean_interactions_sem_df)[j]] <- mean_interactions_sem_df[i, j]
  }
}

mean_interactions_sem_matrix <- square_matrix_sem



# Define species groups dynamically
natives_sem <- c('CLAPUR', 'COLLOM', 'COLLIN', "GILCAP", "NAVSQU", 'EPIDEN', "PLAFIG", "PLECON")  # Replace with foundation species at Blanton Ridge
all_species_sem <- rownames(mean_interactions_sem_matrix)  # All species in the dataset
invasives_sem <- setdiff(all_species_sem, natives)         # Invasives = all species except foundation

# Set up colors for nodes
names(all.sp) <- all_species_sem

all.sp[invasives_sem] <- 'forestgreen'  # Color for invasive species
all.sp[natives_sem] <- 'purple'    # Color for natives

# Plot competition and facilitation separately #####
# Competition only
qgraph(mean_interactions_sem_matrix,
       layout = 'circle',
       posCol = rgb(red = 0, green = 0, blue = 0, alpha = 0),  # Facilitation = transparent
       negCol = 'orange',                                     # Competition = orange
       color = all.sp,
       labels = rownames(mean_interactions_sem_matrix),
       fade = TRUE, directed = TRUE, curveAll = 0.5, label.color = "white", label.cex = 1.25,
       title = 'Competition at South Eugene Meadows', title.cex = 1)

# Facilitation only
qgraph(mean_interactions_sem_matrix,
       layout = 'circle',
       posCol = 'royalblue4',                                 # Facilitation = blue
       negCol = rgb(red = 0, green = 0, blue = 0, alpha = 0), # Competition = transparent
       color = all.sp,
       labels = rownames(mean_interactions_sem_matrix),
       fade = TRUE, directed = TRUE, curveAll = 0.5, label.color = "white", label.cex = 1.25,
       title = 'Facilitation at South Eugene Meadows', title.cex = 1)



# both facil and comp together
qgraph(mean_interactions_sem_matrix,
       layout = 'circle',           # Circular layout
       posCol = 'royalblue4',       # Facilitation = blue
       negCol = 'orange',           # Competition = orange
       color = all.sp,              # Node colors
       labels = rownames(mean_interactions_sem_matrix),  # Node labels
       fade = TRUE, curveAll = 0.5,
       title = 'Competition and Facilitation at South Eugene Meadows', title.cex = 2)


# just looking at our focals from larger plot #####
# Define foundation species
natives <- c('CLAPUR', 'COLLOM', 'COLLIN', "GILCAP", "NAVSQU", 'EPIDEN', "PLAFIG", "PLECON")


# Extract the submatrix: interactions toward foundation species
# Extract submatrix for foundation-to-foundation interactions
native_interactions_sem <- mean_interactions_sem_matrix[natives, natives]

# Plot foundation-to-foundation interactions
qgraph(native_interactions_sem,
       layout = 'circle',
       negCol = 'orange',  # competition = orange
       posCol = 'royalblue4',      # facilitation = green
       color = 'purple',       # All nodes are foundations
       labels = foundation,
       fade = TRUE, directed = TRUE, label.color = "white", label.cex = 1.25,
       title = 'Interactions Among Natives at SEM', title.cex = 1)



# Riverfront ########
## Reading in the data ####
mean_interactions_rf_df <- read.csv("mean_interaction_matrix_rf.csv")

rownames(mean_interactions_rf_df) <- mean_interactions_rf_df$RowNames
mean_interactions_rf_df <- mean_interactions_rf_df[, -which(colnames(mean_interactions_rf_df) == "RowNames")]

all_species_rf <- unique(c(rownames(mean_interactions_rf_df), colnames(mean_interactions_rf_df)))

# based on Bimler figure code!! #######
# Convert to matrix
mean_interactions_rf_matrix <- as.matrix(mean_interactions_rf_df)

# Ensure matrix is correctly oriented
mean_interactions_rf_matrix <- t(mean_interactions_rf_matrix)  # Arrows point towards focal species

# Create a square matrix
square_matrix_rf <- matrix(0, nrow = length(all_species_rf), ncol = length(all_species_rf))
rownames(square_matrix_rf) <- all_species_rf
colnames(square_matrix_rf) <- all_species_rf

# Fill the square matrix with values from mean_interactions_rf_matrix
for (i in seq_len(nrow(mean_interactions_rf_df))) {
  for (j in seq_len(ncol(mean_interactions_rf_df))) {
    square_matrix_rf[rownames(mean_interactions_rf_df)[i], colnames(mean_interactions_rf_df)[j]] <- mean_interactions_rf_df[i, j]
  }
}

mean_interactions_rf_matrix <- square_matrix_rf



# Define species groups dynamically
natives_rf <- c('CLAPUR', 'COLLOM', 'COLLIN', "GILCAP", "NAVSQU", 'EPIDEN', "PLAFIG", "PLECON")  # Replace with foundation species at Blanton Ridge
all_species_rf <- rownames(mean_interactions_rf_matrix)  # All species in the dataset
invasives_rf <- setdiff(all_species_rf, natives)         # Invasives = all species except foundation

# Set up colors for nodes
names(all.sp) <- all_species_rf

all.sp[invasives_rf] <- 'forestgreen'  # Color for invasive species
all.sp[natives_rf] <- 'purple'    # Color for natives

# Plot competition and facilitation separately #####
# Competition only
qgraph(mean_interactions_rf_matrix,
       layout = 'circle',
       posCol = rgb(red = 0, green = 0, blue = 0, alpha = 0),  # Facilitation = transparent
       negCol = 'orange',                                     # Competition = orange
       color = all.sp,
       labels = rownames(mean_interactions_rf_matrix),
       fade = TRUE, directed = TRUE, curveAll = 0.5, label.color = "white", label.cex = 1.25,
       title = 'Competition at the Riverfront', title.cex = 1)

# Facilitation only
qgraph(mean_interactions_rf_matrix,
       layout = 'circle',
       posCol = 'royalblue4',                                 # Facilitation = blue
       negCol = rgb(red = 0, green = 0, blue = 0, alpha = 0), # Competition = transparent
       color = all.sp,
       labels = rownames(mean_interactions_rf_matrix),
       fade = TRUE, directed = TRUE, curveAll = 0.5, label.color = "white", label.cex = 1.25,
       title = 'Facilitation at the Riverfront', title.cex = 1)



# both facil and comp together
qgraph(mean_interactions_rf_matrix,
       layout = 'circle',           # Circular layout
       posCol = 'royalblue4',       # Facilitation = blue
       negCol = 'orange',           # Competition = orange
       color = all.sp,              # Node colors
       labels = rownames(mean_interactions_rf_matrix),  # Node labels
       fade = TRUE, curveAll = 0.5,
       title = 'Competition and Facilitation at the Riverfront', title.cex = 2)


# just looking at our focals from larger plot #####
# Define foundation species
natives <- c('CLAPUR', 'COLLOM', 'COLLIN', "GILCAP", "NAVSQU", 'EPIDEN', "PLAFIG", "PLECON")


# Extract the submatrix: interactions toward foundation species
# Extract submatrix for foundation-to-foundation interactions
native_interactions_rf <- mean_interactions_rf_matrix[natives, natives]

# Plot foundation-to-foundation interactions
qgraph(native_interactions_rf,
       layout = 'circle',
       negCol = 'orange',  # competition = orange
       posCol = 'royalblue4',      # facilitation = green
       color = 'purple',       # All nodes are foundations
       labels = foundation,
       fade = TRUE, directed = TRUE,label.color = "white", label.cex = 1.25,
       title = 'Interactions Among Natives at the Riverfront', title.cex = 1)



# Wild Iris Ridge ########
## Reading in the data ####
mean_interactions_wir_df <- read.csv("mean_interaction_matrix_wir.csv")


rownames(mean_interactions_wir_df) <- mean_interactions_wir_df$RowNames
mean_interactions_wir_df <- mean_interactions_wir_df[, -which(colnames(mean_interactions_wir_df) == "RowNames")]

all_species_wir <- unique(c(rownames(mean_interactions_wir_df), colnames(mean_interactions_wir_df)))

# based on Bimler figure code!! #######
# Convert to matrix
mean_interactions_wir_matrix <- as.matrix(mean_interactions_wir_df)

# Ensure matrix is correctly oriented
mean_interactions_wir_matrix <- t(mean_interactions_wir_matrix)  # Arrows point towards focal species

# Create a square matrix
square_matrix_wir <- matrix(0, nrow = length(all_species_wir), ncol = length(all_species_wir))
rownames(square_matrix_wir) <- all_species_wir
colnames(square_matrix_wir) <- all_species_wir

# Fill the square matrix with values from mean_interactions_wir_matrix
for (i in seq_len(nrow(mean_interactions_wir_df))) {
  for (j in seq_len(ncol(mean_interactions_wir_df))) {
    square_matrix_wir[rownames(mean_interactions_wir_df)[i], colnames(mean_interactions_wir_df)[j]] <- mean_interactions_wir_df[i, j]
  }
}

mean_interactions_wir_matrix <- square_matrix_wir



# Define species groups dynamically
natives_wir <- c('CLAPUR', 'COLLOM', 'COLLIN', "GILCAP", "NAVSQU", 'EPIDEN', "PLAFIG", "PLECON")  # Replace with foundation species at Blanton Ridge
all_species_wir <- rownames(mean_interactions_wir_matrix)  # All species in the dataset
invasives_wir <- setdiff(all_species_wir, natives)         # Invasives = all species except foundation

# Set up colors for nodes
names(all.sp) <- all_species_wir

all.sp[invasives_wir] <- 'forestgreen'  # Color for invasive species
all.sp[natives_wir] <- 'purple'    # Color for natives

# Plot competition and facilitation separately #####
# Competition only
qgraph(mean_interactions_wir_matrix,
       layout = 'circle',
       posCol = rgb(red = 0, green = 0, blue = 0, alpha = 0),  # Facilitation = transparent
       negCol = 'orange',                                     # Competition = orange
       color = all.sp,
       labels = rownames(mean_interactions_wir_matrix),
       fade = TRUE, directed = TRUE, curveAll = 0.5, label.color = "white", label.cex = 1.25,
       title = 'Competition at WIR', title.cex = 1)

# Facilitation only
qgraph(mean_interactions_wir_matrix,
       layout = 'circle',
       posCol = 'royalblue4',                                 # Facilitation = blue
       negCol = rgb(red = 0, green = 0, blue = 0, alpha = 0), # Competition = transparent
       color = all.sp, label.color = "white", label.cex = 1.25,
       labels = rownames(mean_interactions_wir_matrix),
       fade = TRUE, directed = TRUE, curveAll = 0.5,
       title = 'Facilitation at WIR', title.cex = 1)



# both facil and comp together
qgraph(mean_interactions_wir_matrix,
       layout = 'circle',           # Circular layout
       posCol = 'royalblue4',       # Facilitation = blue
       negCol = 'orange',           # Competition = orange
       color = all.sp,              # Node colors
       labels = rownames(mean_interactions_wir_matrix),  # Node labels
       fade = TRUE, curveAll = 0.5, label.color = "white", label.cex = 1.25,
       title = 'Competition and Facilitation at WIR', title.cex = 1)


# just looking at our focals from larger plot #####
# Define foundation species
natives <- c('CLAPUR', 'COLLOM', 'COLLIN', "GILCAP", "NAVSQU", 'EPIDEN', "PLAFIG", "PLECON")


# Extract the submatrix: interactions toward foundation species
# Extract submatrix for foundation-to-foundation interactions
native_interactions_wir <- mean_interactions_wir_matrix[natives, natives]

# Plot foundation-to-foundation interactions
qgraph(native_interactions_wir,
       layout = 'circle',
       negCol = 'orange',  # competition = orange
       posCol = 'royalblue4',      # facilitation = green
       color = 'purple',       # All nodes are foundations
       labels = foundation,
       fade = TRUE, directed = TRUE, label.color = "white", label.cex = 1.25,
       title = 'Interactions Among Natives at WIR', title.cex = 1)



# Plotting something that compares between the parks? #######
library(pheatmap)

# Assuming you have matrices for 4 sites
# keeping the color scheme consistent, blue = facil = positive, red/orange = comp = negative
custom_colors <- colorRampPalette(c("red", "beige", "blue"))(50)
pheatmap(t(mean_interactions_br_df), main = "Blanton Ridge", cluster_rows = FALSE, cluster_cols = FALSE, color = custom_colors)
pheatmap(t(mean_interactions_sem_df), main = "South Eugene Meadows", cluster_rows = FALSE, cluster_cols = FALSE, color = custom_colors)
pheatmap(t(mean_interactions_wir_df), main = "Wild Iris Ridge", cluster_rows = FALSE, cluster_cols = FALSE, color = custom_colors)
pheatmap(t(mean_interactions_rf_df), main = "Riverfront", cluster_rows = FALSE, cluster_cols = FALSE, color = custom_colors)


# Summary Statistics #####
prop_positive_BR <- sum(mean_interactions_br_df > 0) / length(mean_interactions_br_df)
prop_negative_BR <- sum(mean_interactions_br_df < 0) / length(mean_interactions_br_df)

prop_positive_RF <- sum(mean_interactions_rf_df > 0) / length(mean_interactions_rf_df)
prop_negative_RF <- sum(mean_interactions_rf_df < 0) / length(mean_interactions_rf_df)

prop_positive_WIR <- sum(mean_interactions_wir_df > 0) / length(mean_interactions_wir_df)
prop_negative_WIR <- sum(mean_interactions_wir_df < 0) / length(mean_interactions_wir_df)

prop_positive_SEM <- sum(mean_interactions_sem_df > 0) / length(mean_interactions_sem_df)
prop_negative_SEM <- sum(mean_interactions_sem_df < 0) / length(mean_interactions_sem_df)


summary_table <- data.frame(
  Category = c("BR", "RF", "WIR", "SEM"),
  Positive_Proportion = c(prop_positive_BR, prop_positive_RF, prop_positive_WIR, prop_positive_SEM),
  Negative_Proportion = c(prop_negative_BR, prop_negative_RF, prop_negative_WIR, prop_negative_SEM)
)
write.csv(summary_table, "summary_table.csv", row.names = FALSE)

# trying out gt package to create a table png
library(gt)

# Create the table
table_gt <- summary_table %>%
  gt() %>%
  tab_header(title = "Proportion Summary")

# Save as an image
gtsave(data = table_gt, filename = "summary_table.png")



# looking at common species across ALL parks ########
# Get species from each park (assuming row names contain species names)
species_br <- colnames(mean_interactions_br_df)
species_sem <- colnames(mean_interactions_sem_df)
species_rf <- colnames(mean_interactions_rf_df)
species_wir <- colnames(mean_interactions_wir_df)

# Find common species across all parks
common_species <- Reduce(intersect, list(species_br, species_sem, species_rf, species_wir))

# View shared species
print(common_species)

# heatmaps using only the common species

# Filter each dataframe to include only common species
filtered_br_df <- mean_interactions_br_df[colnames(mean_interactions_br_df) %in% common_species, ]
filtered_sem_df <- mean_interactions_sem_df[colnames(mean_interactions_sem_df) %in% common_species, ]
filtered_wir_df <- mean_interactions_wir_df[colnames(mean_interactions_wir_df) %in% common_species, ]
filtered_rf_df <- mean_interactions_rf_df[colnames(mean_interactions_rf_df) %in% common_species, ]


# Generate heatmaps
pheatmap(t(filtered_br_df), main = "Blanton Ridge (Filtered)", cluster_rows = FALSE, cluster_cols = FALSE, color = custom_colors)
pheatmap(t(filtered_sem_df), main = "South Eugene Meadows (Filtered)", cluster_rows = FALSE, cluster_cols = FALSE, color = custom_colors)
pheatmap(t(filtered_wir_df), main = "Wild Iris Ridge (Filtered)", cluster_rows = FALSE, cluster_cols = FALSE, color = custom_colors)
pheatmap(t(filtered_rf_df), main = "Riverfront (Filtered)", cluster_rows = FALSE, cluster_cols = FALSE, color = custom_colors)




# Subset common_species to include only focal species
valid_focal_species <- intersect(common_species, rownames(mean_interactions_br))

# Create interaction data frame
interaction_data <- data.frame(
  Species_A = rep(valid_focal_species, each = length(common_species) * 4), # Focals
  Species_B = rep(common_species, times = length(valid_focal_species) * 4), # Neighbors
  Park = rep(c("BR", "SEM", "RF", "WIR"), each = length(valid_focal_species) * length(common_species)),
  Interaction = c(
    as.vector(mean_interactions_br[valid_focal_species, common_species]),
    as.vector(mean_interactions_sem[valid_focal_species, common_species]),
    as.vector(mean_interactions_rf[valid_focal_species, common_species]),
    as.vector(mean_interactions_wir[valid_focal_species, common_species])
  )
)
# View the interaction data frame
(interaction_data)





# How does the strength of interspecific vs. intraspecific competition differ between each species? #####

# Are there general patterns in species interactions that are the same across the different parks?#####


# Do we see two species with roughly the same interaction coefficients across different parks/environments? ######



# Something I'm just curious about #####
# getting the interactions between daucus carota and each of the natives at WIR
# Define the invasive species and native species
invasive_species <- "DAUCAR"
native_species <- c("COLLOM", "COLLIN", "NAVSQU", "GILCAP", "EPIDEN", "PLAFIG", "PLECON", "CLAPUR") # List of native species

# Subset the interactions of the invasive species with native species
interactions_invasive <- t(mean_interactions_wir)[invasive_species, native_species, drop = FALSE]

# View the isolated interactions
print(interactions_invasive)



# looking at 2 specific species  ######
species_A_B <- data.frame(
  park = c("BR", "SEM", "RF", "WIR"),
  interaction = c(
    mean_interactions_br["GILCAP", "LACTUCA"],
    mean_interactions_sem["GILCAP", "LACTUCA"],
    mean_interactions_rf["GILCAP", "LACTUCA"],
    mean_interactions_wir["GILCAP", "LACTUCA"]
  )
)

print(common_species)


# Plot interaction coefficients across parks
ggplot(species_A_B, aes(x = park, y = interaction)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  ggtitle("Interaction Between GILCAP and LACTUCA Across Parks") +
  ylab("Interaction Coefficient") +
  xlab("Park")


#looking at common species interactions facetplot?
common_species_interactions <- bind_rows(
  mean_interactions_br_df %>%
    rownames_to_column(var = "focal") %>% # Add row names as a column
    pivot_longer(
      cols = -focal,
      names_to = "neighbor", # Second species in interaction
      values_to = "value"               # Interaction strength
    ) %>%
    mutate(park = "BR"),                # Add park identifier
  
  mean_interactions_sem_df %>%
    rownames_to_column(var = "focal") %>%
    pivot_longer(
      cols = -focal,
      names_to = "neighbor",
      values_to = "value"
    ) %>%
    mutate(park = "SEM"),
  
  mean_interactions_rf_df %>%
    rownames_to_column(var = "focal") %>%
    pivot_longer(
      cols = -focal,
      names_to = "neighbor",
      values_to = "value"
    ) %>%
    mutate(park = "RF"),
  
  mean_interactions_wir_df %>%
    rownames_to_column(var = "focal") %>%
    pivot_longer(
      cols = -focal,
      names_to = "neighbor",
      values_to = "value"
    ) %>%
    mutate(park = "WIR")
)

# Check the structure of the combined data
print(head(common_species_interactions))



common_species_interactions <- common_species_interactions %>%
  mutate(interaction = paste(focal, neighbor, sep = " - "))

common_species_interactions <- common_species_interactions %>%
  filter(focal %in% common_species & neighbor %in% common_species)

# Split data into chunks (e.g., 50 interactions per chunk)
chunks <- common_species_interactions %>%
  group_by(interaction) %>%
  group_split() %>%
  split(., ceiling(seq_along(.) / 25))  # Adjust 50 to set number of interactions per plot



# Loop through chunks and create separate plots
for (i in seq_along(chunks)) {
  plot_data <- bind_rows(chunks[[i]])  # Combine rows in the chunk
  
  # Debugging: Print the data passed to the plot
  print(plot_data %>% arrange(interaction, park))
  
  # Create the plot for the current chunk
  p <- ggplot(plot_data, aes(x = park, y = value, fill = park)) +
    geom_bar(stat = "identity", position = "dodge") +   # Bars for each park
    facet_wrap(~ interaction, scales = "free_y") +  # Facet by abbreviated interactions
    theme_minimal() +                                  # Clean theme
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis labels
      strip.text = element_text(size = 10)              # Adjust facet label size
    ) +
    labs(
      title = paste("Interaction Strength Across Parks - Chunk", i),
      x = "Park",
      y = "Interaction Strength",
      fill = "Park"
    )
  
  ggsave(
    filename = paste0("facet_plot_chunk_", i, ".png"),
    plot = p,
    bg = "white",
    width = 16,
    height = 12
  )
}

# running ANOVA #######
# ANOVA test
interactions <- data.frame(
  park = factor(rep(c("BR", "SEM", "RF", "WIR"), each = 1)),
  interaction = c(
    mean_interactions_br["PLAFIG", "PLECON"],
    mean_interactions_sem["PLAFIG", "PLECON"],
    mean_interactions_rf["PLAFIG", "PLECON"],
    mean_interactions_wir["PLAFIG", "PLECON"]
  )
)
aov_result <- aov(interaction ~ park, data = interactions)
summary(aov_result)

# Post hoc test (if significant differences are found)
TukeyHSD(aov_result)


# correlation coeffs across parks?  #######
# looking at posteriors and comparing those across parks? there's so many of them
