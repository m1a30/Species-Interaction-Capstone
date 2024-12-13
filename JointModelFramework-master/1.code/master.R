# Running the joint model on simulated data 

# set up R environment
require(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()) 
library(tidybayes)
library(bayesplot)
library(devtools)

library(rethinking)
library(reshape2)
library(tidyverse)
library(ggplot2)
library(dplyr)

# load required functions
source('data_prep.R')
source('simul_data.R')

# load simulated data 
set.seed(54)


# reading in the different data ####

rf_df <- read.csv("../../Parks data cleaning/data_cleaned/rf_dat.csv")
br_df <- read.csv("../../Parks data cleaning/data_cleaned/br_dat.csv")
sem_df <- read.csv("../../Parks data cleaning/data_cleaned/sem_dat.csv")
wir_df <- read.csv("../../Parks data cleaning/data_cleaned/wir_dat.csv")


## Visualizing what the seed distributions are and Calculating variance of seed counts #####
br_plot <- ggplot(br_df, aes(x = seeds)) +
  geom_histogram(binwidth = 10) +
  ggtitle("Distribution of Seeds BR")

br_len <- length(br_df$seeds)
br_seed_variance <- var(br_df$seeds)
br_seed_summary <- summary(br_df$seeds)

# sem
sem_plot <- ggplot(sem_df, aes(x = seeds)) +
  geom_histogram(binwidth = 10) +
  ggtitle("Distribution of Seeds SEM")

sem_len <- length(sem_df$seeds)
sem_seed_variance <- var(sem_df$seeds)
sem_seed_summary <- summary(sem_df$seeds)

# wir 
wir_plot <- ggplot(wir_df, aes(x = seeds)) +
  geom_histogram(binwidth = 10) +
  ggtitle("Distribution of Seeds WIR")

wir_len <- length (wir_df$seeds)
wir_seed_variance <- var(wir_df$seeds)
wir_seed_summary <- summary(wir_df$seeds)


# rf 
rf_plot <- ggplot(rf_df, aes(x = seeds)) +
  geom_histogram(binwidth = 10) +
  ggtitle("Distribution of Seeds RF")

rf_len <- length (rf_df$seeds)
rf_seed_variance <- var(rf_df$seeds)
rf_seed_summary <- summary(rf_df$seeds)


# Data cleaning to get columns we need #####

df_br <- br_df %>%
  select(3, 14, 17:ncol(br_df))
# #keeping plot species seeds as is, then alphabetizing the neighbors
df_br <- df_br %>%
  select(1:10, sort(names(df_br)[-(1:3)]))

sorted_cols_br <- names(df_br)[3:10] %>% sort()

df_br <- df_br %>%
  select(1:2, all_of(sorted_cols_br), everything())

# # dropping gallium and others that are mostly 0, to up the number to inferrable interactions
drops_br <- c("GALIUM", "SCARLET.PIMPERNEL", "H","J", "D", "G", "I")
df_br <- df_br[ , !(names(df_br) %in% drops_br)]


# SEM data cleaning  ####
df_sem <- sem_df %>%
  select(3, 14, 17:ncol(sem_df))
# #keeping plot species seeds as is, then alphabetizing the neighbors
df_sem <- df_sem %>%
  select(1:10, sort(names(df_sem)[-(1:3)]))

sorted_cols_sem <- names(df_sem)[3:10] %>% sort()

df_sem <- df_sem %>%
  select(1:2, all_of(sorted_cols_sem), everything())
# # dropping gallium
drops_sem <- c("GALIUM")
df_sem <- df_sem[ , !(names(df_sem) %in% drops_sem)]

# RF data cleaning ######

df_rf <- rf_df %>%
  select(3, 14, 16:ncol(rf_df))
# #keeping plot species seeds as is, then alphabetizing the neighbors
df_rf <- df_rf %>%
  select(1:10, sort(names(df_rf)[-(1:3)]))

sorted_cols_rf <- names(df_rf)[3:10] %>% sort()

df_rf <- df_rf %>%
  select(1:2, all_of(sorted_cols_rf), everything())

# decided not to drop these columns as it didn't converge well when I did
# drops_rf <- c("A", "B", "C")
# df_rf <- df_rf[ , !(names(df_rf) %in% drops_rf)]

# # cleaning for wir
df_wir <- wir_df %>%
  select(3, 14, 16:ncol(wir_df))
# #keeping plot species seeds as is, then alphabetizing the neighbors
df_wir <- df_wir %>%
  select(1:10, sort(names(df_wir)[-(1:3)]))

sorted_cols <- names(df_wir)[3:10] %>% sort()

df_wir <- df_wir %>%
  select(1:2, all_of(sorted_cols), everything())

# decided not to drop these columns, this also didn't converge well without these
# drops_wir <- c("A", "B", "C", "D", "E", "F", "G")
# df_wir <- df_wir[ , !(names(df_wir) %in% drops_wir)]


# Data Cleaning: need the seeds to be ints and not doubles ####
df_br$seeds <- as.integer(df_br$seeds)
df_sem$seeds <- as.integer(df_sem$seeds)
df_rf$seeds <- as.integer(df_rf$seeds)
df_wir$seeds <- as.integer(df_wir$seeds)


#filling the alone neighbors that are NAs to 0s 
df_br[is.na(df_br)] <- 0
df_sem[is.na(df_sem)] <- 0
df_rf[is.na(df_rf)] <- 0
df_wir[is.na(df_wir)] <- 0


# Making the dataframes alphabetical by both row and species #########
# NB: if using real data or named species, ensure they are ordered alphabetically in the dataset

# rearranging row names alphabetically 
df_br <- df_br[order(df_br$species), ]
df_sem <- df_sem[order(df_sem$species), ]
df_rf <- df_rf[order(df_rf$species), ]
df_wir <- df_wir[order(df_wir$species), ]


# rearranging column names alphabetically... don't need to do when subsetting down to just the 8 natives
  # but it would be df <- df[order(colnames(df))], then also want to keep the plot, species, seeds in the front still

# 11/24 Changing df$focal to df$species to match our dataframe
# identify focal and neighbouring species to be matched to parameter estimates
focalID_br <- unique(df_br$species)  # this should return the names of unique focal groups in the order
focalID_sem <- unique(df_sem$species)
focalID_rf <- unique(df_rf$species)
focalID_wir <- unique(df_wir$species)


# 11/10: changed the -c(1:2) to -c(1:3), check df for how many non neighbor cols
neighbourID_br <- colnames(df_br[ , -c(1:2)]) # should be ordered focal first (alphabetically), then
neighbourID_sem <- colnames(df_sem[ , -c(1:2)])
neighbourID_rf <- colnames(df_rf[ , -c(1:2)])
neighbourID_wir <- colnames(df_wir[ , -c(1:2)])
# non-focals in alphabetical order

# ensure neighbours are linearly independent across the whole dataset (see S1.2)
N_all_br <- df_br[ , neighbourID_br]
N_all_sem <- df_sem[ , neighbourID_sem]
N_all_rf <- df_rf[ , neighbourID_rf]
N_all_wir <- df_wir[ , neighbourID_wir]


N_all_br <- apply(N_all_br, c(1,2), as.numeric)
N_all_sem <- apply(N_all_sem, c(1,2), as.numeric)
N_all_rf <- apply(N_all_rf, c(1,2), as.numeric)
N_all_wir <- apply(N_all_wir, c(1,2), as.numeric)


X_all_br <- cbind(model.matrix(~as.factor(df_br$species)), N_all_br)
X_all_sem <- cbind(model.matrix(~as.factor(df_sem$species)), N_all_sem)
X_all_rf <- cbind(model.matrix(~as.factor(df_rf$species)), N_all_rf)
X_all_wir <- cbind(model.matrix(~as.factor(df_wir$species)), N_all_wir)

R_all_br <- pracma::rref(X_all_br)
R_all_sem <- pracma::rref(X_all_sem)
R_all_rf <- pracma::rref(X_all_rf)
R_all_wir <- pracma::rref(X_all_wir)


Z_all_br <- t(R_all_br) %*% R_all_br
Z_all_sem <- t(R_all_sem) %*% R_all_sem
Z_all_rf <- t(R_all_rf) %*% R_all_rf
Z_all_wir <- t(R_all_wir) %*% R_all_wir



indep_br <- sapply(seq(1, dim(Z_all_br)[1], 1), function(k){ 
  ifelse(Z_all_br[k, k] == 1 & sum(Z_all_br[k, -k]) == 0, 1, 0)
}) #
all(indep_br == 1) # if TRUE then neighbours are linearly independent and we can continue
if(!all(indep_br == 1)) warning('WARNING neighbours are not linearly independent') 

indep_sem <- sapply(seq(1, dim(Z_all_sem)[1], 1), function(k){ 
  ifelse(Z_all_sem[k, k] == 1 & sum(Z_all_sem[k, -k]) == 0, 1, 0)
}) #
all(indep_sem == 1) # if TRUE then neighbours are linearly independent and we can continue
if(!all(indep_sem == 1)) warning('WARNING neighbours are not linearly independent') 

indep_rf <- sapply(seq(1, dim(Z_all_rf)[1], 1), function(k){ 
  ifelse(Z_all_rf[k, k] == 1 & sum(Z_all_rf[k, -k]) == 0, 1, 0)
}) #
all(indep_rf == 1) # if TRUE then neighbours are linearly independent and we can continue
if(!all(indep_rf == 1)) warning('WARNING neighbours are not linearly independent') 

indep_wir <- sapply(seq(1, dim(Z_all_wir)[1], 1), function(k){ 
  ifelse(Z_all_wir[k, k] == 1 & sum(Z_all_wir[k, -k]) == 0, 1, 0)
}) #
all(indep_wir == 1) # if TRUE then neighbours are linearly independent and we can continue
if(!all(indep_wir == 1)) warning('WARNING neighbours are not linearly independent') 


# Notes: 10/29: changed the nonNcols from 2 to 3 (bc we have the X: 1,2,3 etc column, focal, and seeds column)
# prepare the data into the format required by STAN and the model code
stan.data_br <- data_prep(perform = 'seeds', 
                       focal = 'species', 
                       nonNcols = 2, # number of columns that aren't neighbour abundances
                       df = df_br)

stan.data_sem <- data_prep(perform = 'seeds', 
                          focal = 'species', 
                          nonNcols = 2, # number of columns that aren't neighbour abundances
                          df = df_sem)

stan.data_rf <- data_prep(perform = 'seeds', 
                           focal = 'species', 
                           nonNcols = 2, # number of columns that aren't neighbour abundances
                           df = df_rf)

stan.data_wir <- data_prep(perform = 'seeds', 
                          focal = 'species', 
                          nonNcols = 2, # number of columns that aren't neighbour abundances
                          df = df_wir)


message(paste0('Data dimensions = ', dim(df_br)[1], ', ', dim(df_br)[2]))
message(paste0('Data dimensions = ', dim(df_sem)[1], ', ', dim(df_sem)[2]))
message(paste0('Data dimensions = ', dim(df_rf)[1], ', ', dim(df_rf)[2]))
message(paste0('Data dimensions = ', dim(df_wir)[1], ', ', dim(df_wir)[2]))


message(paste0('Number of focal groups = ', length(focalID_br)))
message(paste0('Number of focal groups = ', length(focalID_sem)))
message(paste0('Number of focal groups = ', length(focalID_rf)))
message(paste0('Number of focal groups = ', length(focalID_wir)))


message(paste0('Number of neighbour groups = ', length(neighbourID_br)))
message(paste0('Number of neighbour groups = ', length(neighbourID_sem)))
message(paste0('Number of neighbour groups = ', length(neighbourID_rf)))
message(paste0('Number of neighbour groups = ', length(neighbourID_wir)))

# for some reason blanton ridge has 0 inferrable interactions?
message(paste0('Proportion of inferrable interactions = ', sum(stan.data_br$Q)/(stan.data_br$S*stan.data_br$`T`)))


message(paste0('Proportion of inferrable interactions = ', sum(stan.data_sem$Q)/(stan.data_sem$S*stan.data_sem$`T`)))
message(paste0('Proportion of inferrable interactions = ', sum(stan.data_rf$Q)/(stan.data_rf$S*stan.data_rf$`T`)))
message(paste0('Proportion of inferrable interactions = ', sum(stan.data_wir$Q)/(stan.data_wir$S*stan.data_wir$`T`)))

stan.seed <- 1234

# BR model fitting ####

fit_br <- stan(file = 'joint_model.stan', 
            data =  stan.data_br,               # named list of data
            chains = 3,
            cores = 3, # added so runs in parallel
            warmup = 2000,          # number of warmup iterations per chain
            iter = 3000,            # total number of iterations per chain
            refresh = 100,         # show progress every 'refresh' iterations
            control = list(max_treedepth = 19,
                           adapt_delta = 0.99), 
            seed = stan.seed
)


# check convergence
print(summary(fit_br, pars=c("gamma_i","ndd_betaij","ri_betaij"))$summary)

# save traceplot
ggsave(path = "../../figures_tables",
  filename = paste0("BR_traceplot_resp_eff", ".png"),
  plot = rstan::traceplot(fit_br, pars=c("response", "effect")),
  width = 16,
  height = 12
)
#save rhat plot
ggsave(path = "../../figures_tables",
         filename = paste0("BR_rhat_plot", ".png"),
       plot = rstan::stan_rhat(fit_br),
       width = 16,
       height = 12
)


# BR posteriors  ####
joint.post.draws_br <- extract.samples(fit_br)
print(str(joint.post.draws_br))

## BR interaction estimates ###
inter_mat_br <- aperm(joint.post.draws_br$ndd_betaij, c(2, 3, 1))
rownames(inter_mat_br) <- focalID_br
colnames(inter_mat_br) <- neighbourID_br

## getting the mean of all the posteriors 
mean_interactions_br <- apply(inter_mat_br, c(1, 2), mean)
mean_interactions_br_df <- as.data.frame(mean_interactions_br)
# keeping row names as separate column for clarity
mean_interactions_br_df$RowNames <- rownames(mean_interactions_br)
# reorder to rownames is 1st col
mean_interactions_br_df <- mean_interactions_br_df[, c("RowNames", setdiff(colnames(mean_interactions_br_df), "RowNames"))]

write.csv(mean_interactions_br_df, "mean_interaction_matrix_br.csv", row.names = FALSE)

# Uncertainty for BR #########
# Reshape ndd_betaij for summarization
ndd_betaij_long_br <- melt(joint.post.draws_br$ndd_betaij, varnames = c("seeds", "species", "neighbor"), value.name = "interaction")

ndd_betaij_long_br <- ndd_betaij_long_br %>%
  mutate(
    species_name_br = focalID_br[species],
    neighbor_name_br = neighbourID_br[neighbor]
  )



# Calculate credible intervals for each focal-neighbor pair
credible_intervals_br <- ndd_betaij_long_br %>%
  group_by(species_name_br, neighbor_name_br) %>%
  summarize(
    lower = quantile(interaction, 0.025),
    upper = quantile(interaction, 0.975),
    mean = mean(interaction)
  )

# Plot credible intervals
BR_credible_intervals <- ggplot(credible_intervals_br, aes(x = neighbor_name_br, y = mean, ymin = lower, ymax = upper, color = as.factor(species_name_br))) +
  geom_point() +
  geom_errorbar(width = 0.2) +
  labs(title = "Mean and 95% Credible Intervals for Interactions at BR",
       x = "Neighbor Species", y = "Interaction Strength",
       color = "Focal Species") + 
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)  # Rotate x-axis labels vertically
  )

# save credible intervals plot
ggsave(path = "../../figures_tables",
       filename = paste0("BR_credible_int_plot", ".png"),
       plot = BR_credible_intervals,
       width = 16,
       height = 12
)


# Summarize interaction strengths by species and neighbor
summary_stats_br <- ndd_betaij_long_br %>%
  group_by(species_name_br, neighbor_name_br) %>%
  summarize(
    mean = mean(interaction),
    lower = quantile(interaction, 0.025),
    upper = quantile(interaction, 0.975),
    .groups = "drop"
  )

print(summary_stats_br)
write.csv(summary_stats_br, "summary_stats_br.csv")



# SEM model fit ####
fit_sem <- stan(file = 'joint_model.stan', 
               data =  stan.data_sem,               # named list of data
               chains = 3,
               cores = 3, # added so runs in parallel
               warmup = 2000,          # number of warmup iterations per chain
               iter = 3000,            # total number of iterations per chain
               refresh = 100,         # show progress every 'refresh' iterations
               control = list(max_treedepth = 19,
                              adapt_delta = 0.99), 
               seed = stan.seed
               )



print(summary(fit_sem, pars=c("gamma_i","ndd_betaij","ri_betaij"))$summary)
SEM_traceplot  <- rstan::traceplot(fit_sem, pars=c("response","effect"))
SEM_rhat <- rstan::stan_rhat(fit_sem, pars=c("response"))

# save traceplot
ggsave(path = "../../figures_tables",
       filename = paste0("SEM_traceplot_resp_eff", ".png"),
       plot = SEM_traceplot,
       width = 16,
       height = 12
)
#save rhat plot
ggsave(path = "../../figures_tables",
       filename = paste0("BR_rhat_plot", ".png"),
       plot = SEM_rhat,
       width = 16,
       height = 12
)


# SEM posteriors ###### 
joint.post.draws_sem <- extract.samples(fit_sem)

## ---  Get interaction estimates ###
inter_mat_sem <- aperm(joint.post.draws_sem$ndd_betaij, c(2, 3, 1))
rownames(inter_mat_sem) <- focalID_sem
colnames(inter_mat_sem) <- neighbourID_sem

## Getting the mean of all the posteriors  #####
mean_interactions_sem <- apply(inter_mat_sem, c(1, 2), mean)
# want to save this to create the visualizations
mean_interactions_sem_df <- as.data.frame(mean_interactions_sem)
# keeping row names as separate column for clarity
mean_interactions_sem_df$RowNames <- rownames(mean_interactions_sem)
# reorder to rownames is 1st col
mean_interactions_sem_df <- mean_interactions_sem_df[, c("RowNames", setdiff(colnames(mean_interactions_sem_df), "RowNames"))]

write.csv(mean_interactions_sem_df, "mean_interaction_matrix_sem.csv", row.names = FALSE)

# Uncertainty for SEM #####

# Reshape ndd_betaij for summarization
ndd_betaij_long_sem <- melt(joint.post.draws_sem$ndd_betaij, varnames = c("seeds", "species", "neighbor"), value.name = "interaction")

ndd_betaij_long_sem <- ndd_betaij_long_sem %>%
  mutate(
    species_name_sem = focalID_sem[species],
    neighbor_name_sem = neighbourID_sem[neighbor]
  )



# Calculate credible intervals for each focal-neighbor pair
credible_intervals_sem <- ndd_betaij_long_sem %>%
  group_by(species_name_sem, neighbor_name_sem) %>%
  summarize(
    lower = quantile(interaction, 0.025),
    upper = quantile(interaction, 0.975),
    mean = mean(interaction)
  )

# Plot credible intervals
SEM_credible_intervals <- ggplot(credible_intervals_sem, aes(x = neighbor_name_sem, y = mean, ymin = lower, ymax = upper, color = as.factor(species_name_sem))) +
  geom_point() +
  geom_errorbar(width = 0.2) +
  labs(title = "Mean and 95% Credible Intervals for Interactions at SEM",
       x = "Neighbor Species", y = "Interaction Strength",
       color = "Focal Species") + 
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)  # Rotate x-axis labels vertically
  )


# save credible intervals plot
ggsave(path = "../../figures_tables",
       filename = paste0("SEM_credible_int_plot", ".png"),
       plot = SEM_credible_intervals,
       width = 16,
       height = 12
)

# Summarize interaction strengths by species and neighbor
summary_stats_sem <- ndd_betaij_long_sem %>%
  group_by(species_name_sem, neighbor_name_sem) %>%
  summarize(
    mean = mean(interaction),
    lower = quantile(interaction, 0.025),
    upper = quantile(interaction, 0.975),
    .groups = "drop"
  )

print(summary_stats_sem)
write.csv(summary_stats_sem, "summary_stats_sem.csv")




# River front model fit ###########
fit_rf <- stan(file = 'joint_model.stan', 
               data =  stan.data_rf,               # named list of data
               chains = 3,
               cores = 3, # added so runs in parallel
               warmup = 2000,          # number of warmup iterations per chain
               iter = 3000,            # total number of iterations per chain
               refresh = 100,         # show progress every 'refresh' iterations
               control = list(max_treedepth = 19,
                              adapt_delta = 0.99), 
               seed = stan.seed
)


print(summary(fit_rf, pars=c("gamma_i","ndd_betaij","ri_betaij"))$summary)
RF_traceplot <- rstan::traceplot(fit_rf, pars=c("response","effect"))
RF_rhat <- rstan::stan_rhat(fit_rf, pars=c("effect"))

# save traceplot
ggsave(path = "../../figures_tables",
       filename = paste0("RF_traceplot_resp_eff", ".png"),
       plot = RF_traceplot,
       width = 16,
       height = 12
)
#save rhat plot
ggsave(path = "../../figures_tables",
       filename = paste0("RF_rhat_plot", ".png"),
       plot = RF_rhat,
       width = 16,
       height = 12
)



# RF posteriors #####
joint.post.draws_rf <- extract.samples(fit_rf)

## ---  Get interaction estimates ###
inter_mat_rf <- aperm(joint.post.draws_rf$ndd_betaij, c(2, 3, 1))
rownames(inter_mat_rf) <- focalID_rf
colnames(inter_mat_rf) <- neighbourID_rf

## RF mean of posteriors  #####
mean_interactions_rf <- apply(inter_mat_rf, c(1, 2), mean)
# want to save this to create the visualizations
mean_interactions_rf_df <- as.data.frame(mean_interactions_rf)
# keeping row names as separate column for clarity
mean_interactions_rf_df$RowNames <- rownames(mean_interactions_rf)
# reorder to rownames is 1st col
mean_interactions_rf_df <- mean_interactions_rf_df[, c("RowNames", setdiff(colnames(mean_interactions_rf_df), "RowNames"))]

write.csv(mean_interactions_rf_df, "mean_interaction_matrix_rf.csv", row.names = FALSE)


# uncertainty ########

# Reshape ndd_betaij for summarization
ndd_betaij_long_rf <- melt(joint.post.draws_rf$ndd_betaij, varnames = c("seeds", "species", "neighbor"), value.name = "interaction")

ndd_betaij_long_rf <- ndd_betaij_long_rf %>%
  mutate(
    species_name_rf = focalID_rf[species],
    neighbor_name_rf = neighbourID_rf[neighbor]
  )


# Calculate credible intervals for each focal-neighbor pair
credible_intervals_rf <- ndd_betaij_long_rf %>%
  group_by(species_name_rf, neighbor_name_rf) %>%
  summarize(
    lower = quantile(interaction, 0.025),
    upper = quantile(interaction, 0.975),
    mean = mean(interaction)
  )

# Plot credible intervals
RF_credible_intervals <- ggplot(credible_intervals_rf, aes(x = neighbor_name_rf, y = mean, ymin = lower, ymax = upper, color = as.factor(species_name_rf))) +
  geom_point() +
  geom_errorbar(width = 0.2) +
  labs(title = "Mean and 95% Credible Intervals for Interactions at RF",
       x = "Neighbor Species", y = "Interaction Strength",
       color = "Focal Species") + 
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)  # Rotate x-axis labels vertically
  )


# save credible intervals plot
ggsave(path = "../../figures_tables",
       filename = paste0("RF_credible_int_plot", ".png"),
       plot = RF_credible_intervals,
       width = 16,
       height = 12
)


# Summarize interaction strengths by species and neighbor
summary_stats_rf <- ndd_betaij_long_rf %>%
  group_by(species_name_rf, neighbor_name_rf) %>%
  summarize(
    mean = mean(interaction),
    lower = quantile(interaction, 0.025),
    upper = quantile(interaction, 0.975),
    .groups = "drop"
  )

print(summary_stats_rf)
write.csv(summary_stats_rf, "summary_stats_rf.csv")




# WIR model fit #########

fit_wir <- stan(file = 'joint_model.stan', 
                data =  stan.data_wir,               # named list of data
                chains = 3,
                cores = 3, # added so runs in parallel
                warmup = 2000,          # number of warmup iterations per chain
                iter = 3000,            # total number of iterations per chain
                refresh = 100,         # show progress every 'refresh' iterations
                control = list(max_treedepth = 19,
                               adapt_delta = 0.99), 
                seed = stan.seed
)


#pairs(fit_wir, pars = c("response[1]", "effect[1]"))

print(summary(fit_wir, pars=c("gamma_i","ndd_betaij","ri_betaij"))$summary)
WIR_traceplot <- rstan::traceplot(fit_wir, pars=c("response","effect"))

ggsave(path = "../../figures_tables",
       filename = paste0("WIR_traceplot", ".png"),
       plot = WIR_traceplot,
       width = 16,
       height = 12
)

WIR_rhat <- rstan::stan_rhat(fit_wir, pars=c("effect"))
ggsave(path = "../../figures_tables",
       filename = paste0("WIR_rhat_plot", ".png"),
       plot = WIR_rhat,
       width = 16,
       height = 12
)


# WIR posteriors ###### 
joint.post.draws_wir <- extract.samples(fit_wir)

## ---  Get interaction estimates ###
inter_mat_wir <- aperm(joint.post.draws_wir$ndd_betaij, c(2, 3, 1))
rownames(inter_mat_wir) <- focalID_wir
colnames(inter_mat_wir) <- neighbourID_wir

## WIR mean of all the posteriors  #####
mean_interactions_wir <- apply(inter_mat_wir, c(1, 2), mean)
# want to save this to create the visualizations
mean_interactions_wir_df <- as.data.frame(mean_interactions_wir)
# keeping row names as separate column for clarity
mean_interactions_wir_df$RowNames <- rownames(mean_interactions_wir)
# reorder to rownames is 1st col
mean_interactions_wir_df <- mean_interactions_wir_df[, c("RowNames", setdiff(colnames(mean_interactions_wir_df), "RowNames"))]

write.csv(mean_interactions_wir_df, "mean_interaction_matrix_wir.csv", row.names = FALSE)



# uncertainty ########

# Reshape ndd_betaij for summarization
ndd_betaij_long_wir <- melt(joint.post.draws_wir$ndd_betaij, varnames = c("seeds", "species", "neighbor"), value.name = "interaction")

ndd_betaij_long_wir <- ndd_betaij_long_wir %>%
  mutate(
    species_name_wir = focalID_wir[species],
    neighbor_name_wir = neighbourID_wir[neighbor]
  )


# Calculate credible intervals for each focal-neighbor pair
credible_intervals_wir <- ndd_betaij_long_wir %>%
  group_by(species_name_wir, neighbor_name_wir) %>%
  summarize(
    lower = quantile(interaction, 0.025),
    upper = quantile(interaction, 0.975),
    mean = mean(interaction)
  )

# Plot credible intervals
wir_credible_intervals <- ggplot(credible_intervals_wir, aes(x = neighbor_name_wir, y = mean, ymin = lower, ymax = upper, color = as.factor(species_name_wir))) +
  geom_point() +
  geom_errorbar(width = 0.2) +
  labs(title = "Mean and 95% Credible Intervals for Interactions at WIR",
       x = "Neighbor Species", y = "Interaction Strength",
       color = "Focal Species") + 
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)  # Rotate x-axis labels vertically
  )


# save credible intervals plot
ggsave(path = "../../figures_tables",
       filename = paste0("WIR_credible_int_plot", ".png"),
       plot = wir_credible_intervals,
       width = 16,
       height = 12
)


# Summarize interaction strengths by species and neighbor
summary_stats_wir <- ndd_betaij_long_wir %>%
  group_by(species_name_wir, neighbor_name_wir) %>%
  summarize(
    mean = mean(interaction),
    lower = quantile(interaction, 0.025),
    upper = quantile(interaction, 0.975),
    .groups = "drop"
  )

print(summary_stats_wir)
write.csv(summary_stats_wir, "summary_stats_WIR.csv")



