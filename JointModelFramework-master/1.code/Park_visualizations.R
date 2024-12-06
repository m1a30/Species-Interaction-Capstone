
## Visualizations ####
library(qgraph)

# Blanton Ridge #######
## Reading in the data ####
mean_interactions_br_df <- read.csv("mean_interaction_matrix_br.csv")

### trying out a visualization  ##
rownames(mean_interactions_br_df) <- mean_interactions_br_df$RowNames
mean_interactions_br_df <- mean_interactions_br_df[, -which(colnames(mean_interactions_br_df) == "RowNames")]

all_species_br <- unique(c(rownames(mean_interactions_br_df), colnames(mean_interactions_br_df)))


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


# correlation coeffs across parks? #######
