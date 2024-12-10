
## Visualizations ####
library(qgraph)

# Blanton Ridge #######
## Reading in the data ####
mean_interactions_br_df <- read.csv("mean_interaction_matrix_br.csv")

### trying out a visualization  ##
rownames(mean_interactions_br_df) <- mean_interactions_br_df$RowNames
mean_interactions_br_df <- mean_interactions_br_df[, -which(colnames(mean_interactions_br_df) == "RowNames")]

all_species_br <- unique(c(rownames(mean_interactions_br_df), colnames(mean_interactions_br_df)))

# based on Bimler figure code!! #######
# Load your data
mean_interactions_br_df <- read.csv("mean_interaction_matrix_br.csv")

# Convert to matrix
rownames(mean_interactions_br_df) <- mean_interactions_br_df$RowNames
mean_interactions_br_df <- mean_interactions_br_df[, -which(colnames(mean_interactions_br_df) == "RowNames")]
mean_interactions_br_matrix <- as.matrix(mean_interactions_br_df)

# Ensure matrix is correctly oriented
mean_interactions_br_matrix <- t(mean_interactions_br_matrix)  # Arrows point towards focal species
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
all.sp <- rep('white', length(all_species))
names(all.sp) <- all_species

all.sp[invasives] <- 'firebrick3'  # Color for invasive species
all.sp[foundation] <- 'purple'    # Color for foundation species

# Plot competition and facilitation separately
library(qgraph)

png('Blanton_Ridge_Interactions_No_Keystone.png', width = 1875, height = 3750, units = 'px')
par(mfrow = c(2, 1))

# Competition only
qgraph(mean_interactions_br_matrix,
       layout = 'circle',
       negCol = rgb(red = 0, green = 0, blue = 0, alpha = 0),  # Facilitation = transparent
       posCol = 'orange',                                     # Competition = orange
       color = all.sp,
       labels = rownames(mean_interactions_br_matrix),
       fade = TRUE, directed = TRUE,
       title = 'A: Competition at Blanton Ridge', title.cex = 2)

# Facilitation only
qgraph(mean_interactions_br_matrix,
       layout = 'circle',
       negCol = 'royalblue4',                                 # Facilitation = blue
       posCol = rgb(red = 0, green = 0, blue = 0, alpha = 0), # Competition = transparent
       color = all.sp,
       labels = rownames(mean_interactions_br_matrix),
       fade = TRUE, directed = TRUE,
       title = 'B: Facilitation at Blanton Ridge', title.cex = 1)



# both facil and comp together
qgraph(mean_interactions_br_matrix,
       layout = 'circle',           # Circular layout
       negCol = 'royalblue4',       # Facilitation = blue
       posCol = 'orange',           # Competition = orange
       color = all.sp,              # Node colors
       labels = rownames(mean_interactions_br_matrix),  # Node labels
       fade = TRUE, directed = TRUE, # Directed edges
       title = 'Competition and Facilitation at Blanton Ridge', title.cex = 2)


# just looking at our focals from larger plot #####
# Define foundation species
foundation <- c('CLAPUR', 'COLLOM', 'COLLIN', "GILCAP", "NAVSQU", 'EPIDEN', "PLAFIG", "PLECON")


# Extract the submatrix: interactions toward foundation species
# Extract submatrix for foundation-to-foundation interactions
foundation_interactions <- mean_interactions_br_matrix[foundation, foundation]

# Plot foundation-to-foundation interactions
qgraph(foundation_interactions,
       layout = 'circle',
       negCol = 'red',  # Facilitation = blue
       posCol = 'forestgreen',      # Competition = orange
       color = 'purple',       # All nodes are foundations
       labels = foundation,
       fade = TRUE, directed = TRUE,
       title = 'Interactions Among Foundation Species', title.cex = 2)


effects_on_foundation <- mean_interactions_br_matrix[, foundation]
qgraph(effects_on_foundation,
       layout = 'circle',
       negCol = 'royalblue4',  # Facilitation = blue
       posCol = 'orange',      # Competition = orange
       color = all.sp,
       labels = rownames(effects_on_foundation),
       fade = TRUE, directed = TRUE,
       title = 'Effects on Foundation Species', title.cex = 2)


# plotting competition and facilitation on same graph 
competition_plot_br <- qgraph(
  mean_interactions_br_df,  # use the mean_interactions matrix
  layout = 'circle', 
  posCol = "forestgreen",  # Facilitation = green
  negCol = 'red3',  # Competition = red
 # color = rep("skyblue", length(all_species)),  # Node colors
  labels = all_species_br,  # Use row names of the matrix for labels
  title = 'Competition and Facilitation at Blanton Ridge', 
  title.cex = 1.5
)


# South Eugene Meadows #######
## Reading in the data ####
mean_interactions_sem_df <- read.csv("mean_interaction_matrix_sem.csv")

### trying out a visualization  ##
rownames(mean_interactions_sem_df) <- mean_interactions_sem_df$RowNames
mean_interactions_sem_df <- mean_interactions_sem_df[, -which(colnames(mean_interactions_sem_df) == "RowNames")]

all_species_sem <- unique(c(rownames(mean_interactions_sem_df), colnames(mean_interactions_sem_df)))


# plotting competition and facilitation on same graph 
competition_plot_sem <- qgraph(
  mean_interactions_sem_df,  # use the mean_interactions matrix
  layout = 'circle', 
  posCol = "forestgreen",  # Facilitation = green
  negCol = 'red3',  # Competition = red
  # color = rep("skyblue", length(all_species)),  # Node colors
  labels = all_species_sem,  # Use row names of the matrix for labels
  title = 'Competition and Facilitation at South Eugene Meadows', 
  title.cex = 1.5
)

# Riverfront ########
## Reading in the data ####
mean_interactions_rf_df <- read.csv("mean_interaction_matrix_rf.csv")

### trying out a visualization  ##
rownames(mean_interactions_rf_df) <- mean_interactions_rf_df$RowNames
mean_interactions_rf_df <- mean_interactions_rf_df[, -which(colnames(mean_interactions_rf_df) == "RowNames")]

all_species_rf <- unique(c(rownames(mean_interactions_rf_df), colnames(mean_interactions_rf_df)))


# plotting competition and facilitation on same graph 
competition_plot_rf <- qgraph(
  mean_interactions_rf_df,  # use the mean_interactions matrix
  layout = 'circle', 
  posCol = "forestgreen",  # Facilitation = green
  negCol = 'red3',  # Competition = red
  # color = rep("skyblue", length(all_species)),  # Node colors
  labels = all_species_rf,  # Use row names of the matrix for labels
  title = 'Competition and Facilitation at Riverfront', 
  title.cex = 1.5
)

# Wild Iris Ridge ########
## Reading in the data ####
mean_interactions_wir_df <- read.csv("mean_interaction_matrix_wir.csv")

### trying out a visualization  ##
rownames(mean_interactions_wir_df) <- mean_interactions_wir_df$RowNames
mean_interactions_wir_df <- mean_interactions_wir_df[, -which(colnames(mean_interactions_wir_df) == "RowNames")]

all_species_wir <- unique(c(rownames(mean_interactions_wir_df), colnames(mean_interactions_wir_df)))


# plotting competition and facilitation on same graph 
competition_plot_wir <- qgraph(
  mean_interactions_wir_df,  # use the mean_interactions matrix
  layout = 'circle', 
  posCol = "forestgreen",  # Facilitation = green
  negCol = 'red3',  # Competition = red
  # color = rep("skyblue", length(all_species)),  # Node colors
  labels = all_species_wir,  # Use row names of the matrix for labels
  title = 'Competition and Facilitation at Wild Iris Ridge', 
  title.cex = 1.5
)



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
