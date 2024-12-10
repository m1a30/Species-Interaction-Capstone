
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
invasives <- setdiff(all_species, foundation)         # Invasives = all species except foundation

# Set up colors for nodes
names(all.sp) <- all_species

all.sp[invasives] <- 'forestgreen'  # Color for invasive species
all.sp[foundation] <- 'purple'    # Color for foundation species

# Plot competition and facilitation separately #####
# Competition only
qgraph(mean_interactions_br_matrix,
       layout = 'circle',
       posCol = rgb(red = 0, green = 0, blue = 0, alpha = 0),  # Facilitation = transparent
       negCol = 'red3',                                     # Competition = red
       color = all.sp,
       labels = rownames(mean_interactions_br_matrix),
       fade = TRUE, directed = TRUE, curve = 2,
       title = 'Competition at Blanton Ridge', title.cex = 2)

# Facilitation only
qgraph(mean_interactions_br_matrix,
       layout = 'circle',
       posCol = 'royalblue4',                                 # Facilitation = blue
       negCol = rgb(red = 0, green = 0, blue = 0, alpha = 0), # Competition = transparent
       color = all.sp,
       labels = rownames(mean_interactions_br_matrix),
       fade = TRUE, directed = TRUE,
       title = 'Facilitation at Blanton Ridge', title.cex = 1)



# both facil and comp together
qgraph(mean_interactions_br_matrix,
       layout = 'circle',           # Circular layout
       posCol = 'royalblue4',       # Facilitation = blue
       negCol = 'orange',           # Competition = orange
       color = all.sp,              # Node colors
       labels = rownames(mean_interactions_br_matrix),  # Node labels
       fade = TRUE, curve = 2,
       title = 'Competition and Facilitation at Blanton Ridge', title.cex = 2)


# just looking at our focals from larger plot #####
# Define foundation species
foundation <- c('CLAPUR', 'COLLOM', 'COLLIN', "GILCAP", "NAVSQU", 'EPIDEN', "PLAFIG", "PLECON")


# Extract the submatrix: interactions toward foundation species
# Extract submatrix for foundation-to-foundation interactions
native_interactions <- mean_interactions_br_matrix[foundation, foundation]

# Plot foundation-to-foundation interactions
qgraph(native_interactions,
       layout = 'circle',
       negCol = 'red',  # competition = red
       posCol = 'forestgreen',      # facilitation = green
       color = 'purple',       # All nodes are foundations
       labels = foundation,
       fade = TRUE, directed = TRUE,
       title = 'Interactions Among Natives at BR', title.cex = 2)



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

all.sp[invasives] <- 'forestgreen'  # Color for invasive species
all.sp[natives_sem] <- 'purple'    # Color for natives

# Plot competition and facilitation separately #####
# Competition only
qgraph(mean_interactions_sem_matrix,
       layout = 'circle',
       posCol = rgb(red = 0, green = 0, blue = 0, alpha = 0),  # Facilitation = transparent
       negCol = 'orange',                                     # Competition = orange
       color = all.sp,
       labels = rownames(mean_interactions_sem_matrix),
       fade = TRUE, directed = TRUE, curveAll = 0.5,
       title = 'Competition at South Eugene Meadows', title.cex = 2)

# Facilitation only
qgraph(mean_interactions_sem_matrix,
       layout = 'circle',
       posCol = 'royalblue4',                                 # Facilitation = blue
       negCol = rgb(red = 0, green = 0, blue = 0, alpha = 0), # Competition = transparent
       color = all.sp,
       labels = rownames(mean_interactions_sem_matrix),
       fade = TRUE, directed = TRUE, curveAll = 0.5,
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
       negCol = 'red',  # competition = red
       posCol = 'forestgreen',      # facilitation = green
       color = 'purple',       # All nodes are foundations
       labels = foundation,
       fade = TRUE, directed = TRUE,
       title = 'Interactions Among Natives at SEM', title.cex = 2)



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

all.sp[invasives] <- 'forestgreen'  # Color for invasive species
all.sp[natives_rf] <- 'purple'    # Color for natives

# Plot competition and facilitation separately #####
# Competition only
qgraph(mean_interactions_rf_matrix,
       layout = 'circle',
       posCol = rgb(red = 0, green = 0, blue = 0, alpha = 0),  # Facilitation = transparent
       negCol = 'orange',                                     # Competition = orange
       color = all.sp,
       labels = rownames(mean_interactions_rf_matrix),
       fade = TRUE, directed = TRUE, curveAll = 0.5,
       title = 'Competition at South Eugene Meadows', title.cex = 2)

# Facilitation only
qgraph(mean_interactions_rf_matrix,
       layout = 'circle',
       posCol = 'royalblue4',                                 # Facilitation = blue
       negCol = rgb(red = 0, green = 0, blue = 0, alpha = 0), # Competition = transparent
       color = all.sp,
       labels = rownames(mean_interactions_rf_matrix),
       fade = TRUE, directed = TRUE, curveAll = 0.5,
       title = 'Facilitation at South Eugene Meadows', title.cex = 1)



# both facil and comp together
qgraph(mean_interactions_rf_matrix,
       layout = 'circle',           # Circular layout
       posCol = 'royalblue4',       # Facilitation = blue
       negCol = 'orange',           # Competition = orange
       color = all.sp,              # Node colors
       labels = rownames(mean_interactions_rf_matrix),  # Node labels
       fade = TRUE, curveAll = 0.5,
       title = 'Competition and Facilitation at South Eugene Meadows', title.cex = 2)


# just looking at our focals from larger plot #####
# Define foundation species
natives <- c('CLAPUR', 'COLLOM', 'COLLIN', "GILCAP", "NAVSQU", 'EPIDEN', "PLAFIG", "PLECON")


# Extract the submatrix: interactions toward foundation species
# Extract submatrix for foundation-to-foundation interactions
native_interactions_rf <- mean_interactions_rf_matrix[natives, natives]

# Plot foundation-to-foundation interactions
qgraph(native_interactions_rf,
       layout = 'circle',
       negCol = 'red',  # competition = red
       posCol = 'forestgreen',      # facilitation = green
       color = 'purple',       # All nodes are foundations
       labels = foundation,
       fade = TRUE, directed = TRUE,
       title = 'Interactions Among Natives at rf', title.cex = 2)



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

all.sp[invasives] <- 'forestgreen'  # Color for invasive species
all.sp[natives_wir] <- 'purple'    # Color for natives

# Plot competition and facilitation separately #####
# Competition only
qgraph(mean_interactions_wir_matrix,
       layout = 'circle',
       posCol = rgb(red = 0, green = 0, blue = 0, alpha = 0),  # Facilitation = transparent
       negCol = 'orange',                                     # Competition = orange
       color = all.sp,
       labels = rownames(mean_interactions_wir_matrix),
       fade = TRUE, directed = TRUE, curveAll = 0.5,
       title = 'Competition at South Eugene Meadows', title.cex = 1)

# Facilitation only
qgraph(mean_interactions_wir_matrix,
       layout = 'circle',
       posCol = 'royalblue4',                                 # Facilitation = blue
       negCol = rgb(red = 0, green = 0, blue = 0, alpha = 0), # Competition = transparent
       color = all.sp,
       labels = rownames(mean_interactions_wir_matrix),
       fade = TRUE, directed = TRUE, curveAll = 0.5,
       title = 'Facilitation at South Eugene Meadows', title.cex = 1)



# both facil and comp together
qgraph(mean_interactions_wir_matrix,
       layout = 'circle',           # Circular layout
       posCol = 'royalblue4',       # Facilitation = blue
       negCol = 'orange',           # Competition = orange
       color = all.sp,              # Node colors
       labels = rownames(mean_interactions_wir_matrix),  # Node labels
       fade = TRUE, curveAll = 0.5,
       title = 'Competition and Facilitation at South Eugene Meadows', title.cex = 1)


# just looking at our focals from larger plot #####
# Define foundation species
natives <- c('CLAPUR', 'COLLOM', 'COLLIN', "GILCAP", "NAVSQU", 'EPIDEN', "PLAFIG", "PLECON")


# Extract the submatrix: interactions toward foundation species
# Extract submatrix for foundation-to-foundation interactions
native_interactions_wir <- mean_interactions_wir_matrix[natives, natives]

# Plot foundation-to-foundation interactions
qgraph(native_interactions_wir,
       layout = 'circle',
       negCol = 'red',  # competition = red
       posCol = 'forestgreen',      # facilitation = green
       color = 'purple',       # All nodes are foundations
       labels = foundation,
       fade = TRUE, directed = TRUE,
       title = 'Interactions Among Natives at wir', title.cex = 1)



# Plotting something that compares between the parks? #######
library(pheatmap)

# Assuming you have matrices for 4 sites
pheatmap(mean_interactions_br_df, main = "Blanton Ridge")
pheatmap(mean_interactions_sem_df, main = "South Eugene Meadows")
pheatmap(mean_interactions_wir_df, main = "Wild Iris Ridge")
pheatmap(mean_interactions_rf_df, main = "Riverfront")

# Combined heatmap, not working ########
combined_matrix <- rbind(
  as.vector(mean_interactions_br_df),
  as.vector(mean_interactions_sem_df),
  as.vector(mean_interactions_wir_df),
  as.vector(mean_interactions_rf_df)
)
rownames(combined_matrix) <- c("BR", "SEM", "WIR", "RF")
pheatmap(combined_matrix, cluster_rows = TRUE, cluster_cols = TRUE)

# combined qgraph, not working#####
layout <- qgraph::averageLayout(mean_interactions_br_df, mean_interactions_sem_df, mean_interactions_wir_df, mean_interactions_rf_df)
qgraph(mean_interactions_br_df, layout = layout, title = "Site 1")
qgraph(mean_interactions_sem_df, layout = layout, title = "Site 2")
qgraph(mean_interactions_wir_df, layout = layout, title = "Site 3")
qgraph(mean_interactions_rf_df, layout = layout, title = "Site 4")

# Summary Statistics #####
prop_positive_BR <- sum(mean_interactions_br_df > 0) / length(mean_interactions_br_df)
prop_negative_BR <- sum(mean_interactions_br_df < 0) / length(mean_interactions_br_df)

prop_positive_WIR <- sum(mean_interactions_wir_df > 0) / length(mean_interactions_wir_df)
prop_negative_WIR <- sum(mean_interactions_wir_df < 0) / length(mean_interactions_wir_df)

# looking at 2 specific species  ######
species_A_B <- data.frame(
  park = c("BR", "SEM", "RF", "WIR"),
  interaction = c(
    mean_interactions_br["PLAFIG", "PLECON"],
    mean_interactions_sem["PLAFIG", "PLECON"],
    mean_interactions_rf["PLAFIG", "PLECON"],
    mean_interactions_wir["PLAFIG", "PLECON"]
  )
)

# Plot interaction coefficients across parks
ggplot(species_A_B, aes(x = park, y = interaction)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  ggtitle("Interaction Between PLAFIG and PLECON Across Parks") +
  ylab("Interaction Coefficient") +
  xlab("Park")

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
