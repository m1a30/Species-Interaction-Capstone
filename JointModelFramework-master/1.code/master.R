# Running the joint model on simulated data 

# set up R environment
require(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()) 
#rstan_options(auto_write = TRUE)

#install.packages('tidybayes')
# getting the rethinking package to download on Jeff's laptop
# https://www.rdocumentation.org/packages/rethinking/versions/1.59
# https://github.com/rmcelreath/rethinking/issues/174
library(tidybayes)
library(bayesplot)
library(devtools)
#devtools::install_github("rmcelreath/rethinking")
library(rethinking)
library(reshape2)
library(tidyverse)
library(ggplot2)

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


# Code for getting rid of the max seed values... ####

#  thresholds
# upper_bound <- quantile(sem_df$seeds, 0.94)  # 98th percentile
# 
# # Filter out outliers
# df_sem <- sem_df %>%
#   filter(seeds <= upper_bound)


# seeing how different the output is when getting rid of large seed counts
#SEM Identify the 20 largest values in the seeds column  #######
top_20_max_sem <- sem_df %>%
  arrange(desc(seeds)) %>%
  slice_head(n = 20)
#Filter out the rows with the 10 largest values
sem_df <- sem_df %>%
  filter(!seeds %in% top_20_max_sem$seeds)

# RF removing the max 20 vaues of seed counts ######
top_20_max_rf <- rf_df %>%
  arrange(desc(seeds)) %>%
  slice_head(n = 20)
#Filter out the rows with the 10 largest values
rf_df <- rf_df %>%
  filter(!seeds %in% top_20_max_rf$seeds)


# WIR removing the max 10 vaues of seed counts ######
top_10_max_wir <- wir_df %>%
  arrange(desc(seeds)) %>%
  slice_head(n = 10)
#Filter out the rows with the 10 largest values
wir_df <- wir_df %>%
  filter(!seeds %in% top_10_max_wir$seeds)

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

# Running the model with only focal species ####
df_br <- br_df %>% dplyr::select(species, plot, seeds, CLAPUR, COLLIN, COLLOM, EPIDEN, GILCAP, NAVSQU, PLAFIG, PLECON)
df_sem <- sem_df %>% dplyr::select(species, plot, seeds, CLAPUR, COLLIN, COLLOM, EPIDEN, GILCAP, NAVSQU, PLAFIG, PLECON)
df_rf <- rf_df %>% dplyr::select(species, plot, seeds, CLAPUR, COLLIN, COLLOM, EPIDEN, GILCAP, NAVSQU, PLAFIG, PLECON)
df_wir <- wir_df %>% dplyr::select(species, plot, seeds, CLAPUR, COLLIN, COLLOM, EPIDEN, GILCAP, NAVSQU, PLAFIG, PLECON)

# Code for if running ALL neighbors, subsetting down to the columns we want ####
# If we run the data with all of the weeds, then consider getting rid of the gallium (the 1-10 scale) and grass species, not many of them? ####
# modifying data frame so it's just plot species seeds and all neighbors
# df <- df %>%
#   select(2:3, 14, 17:ncol(df))
# keeping plot species seeds as is, then alphabetizing the neighbors
# df<- df %>%
#   select(1:3, sort(names(df)[-(1:3)]))

# Data Cleaning: need the seeds to be ints and not doubles? I think that's what's causing the error ####
df_br$seeds <- as.integer(df_br$seeds)
df_sem$seeds <- as.integer(df_sem$seeds)
df_rf$seeds <- as.integer(df_rf$seeds)
df_wir$seeds <- as.integer(df_wir$seeds)


sum(is.na(df_sem))

#filling the alone neighbors that are NAs to 0s 
df_br[is.na(df_br)] <- 0
df_sem[is.na(df_sem)] <- 0
df_rf[is.na(df_rf)] <- 0
df_wir[is.na(df_wir)] <- 0



# NB: if using real data or named species, ensure they are ordered alphabetically in the dataset
# 11/24 Changing df$focal to df$species to match our dataframe
# identify focal and neighbouring species to be matched to parameter estimates
focalID_br <- unique(df_br$species)  # this should return the names of unique focal groups in the order
focalID_sem <- unique(df_sem$species)
focalID_rf <- unique(df_rf$species)
focalID_wir <- unique(df_wir$species)


# extra data cleaning step, need the focals to be alphabetized ####
focalID_br <- sort(focalID_br)
focalID_sem <- sort(focalID_sem)
focalID_rf <- sort(focalID_rf)
focalID_wir <- sort(focalID_wir)
# in which they are encountered in the dataframe - must be alphabetical

# 11/10: changed the -c(1:2) to -c(1:3), check df for how many non neighbor cols
neighbourID_br <- colnames(df_br[ , -c(1:3)]) # should be ordered focal first (alphabetically), then
neighbourID_sem <- colnames(df_sem[ , -c(1:3)])
neighbourID_rf <- colnames(df_rf[ , -c(1:3)])
neighbourID_wir <- colnames(df_wir[ , -c(1:3)])
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
                       nonNcols = 3, # number of columns that aren't neighbour abundances
                       df = df_br)

stan.data_sem <- data_prep(perform = 'seeds', 
                          focal = 'species', 
                          nonNcols = 3, # number of columns that aren't neighbour abundances
                          df = df_sem)

stan.data_rf <- data_prep(perform = 'seeds', 
                           focal = 'species', 
                           nonNcols = 3, # number of columns that aren't neighbour abundances
                           df = df_rf)

stan.data_wir <- data_prep(perform = 'seeds', 
                          focal = 'species', 
                          nonNcols = 3, # number of columns that aren't neighbour abundances
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

# extracting the values for br bc we know that the model likes this data! #####
br_fit <- rstan::extract(fit_br)


## Code to look at where the high rhat values are ####
fit_summary <- as.data.frame(summary(fit)$summary)
# high_rhat <- fit_summary[fit_summary$Rhat > 1.1, ]
# low_ess <- fit_summary[fit_summary$n_eff < 100, ]
# View(low_ess)
# View(high_rhat)
# was finding that response and effect were high values 
bayesplot::mcmc_pairs(fit_sem, pars = c("response[1]", "response[2]", "effect[1]", "effect[2]"))


# check convergence
print(summary(fit_br, pars=c("gamma_i","ndd_betaij","ri_betaij"))$summary)
rstan::traceplot(fit_br, pars=c("gamma_i","ndd_betaij"))
rstan::stan_rhat(fit_br)

traceplot(fit_br, pars=c("effect"))


# Get the full posteriors 
joint.post.draws_br <- extract.samples(fit_br)

## ---  Get interaction estimates ###
inter_mat_br <- aperm(joint.post.draws_br$ndd_betaij, c(2, 3, 1))
rownames(inter_mat_br) <- focalID_br
colnames(inter_mat_br) <- neighbourID_br

## getting the mean of all the posteriors 
mean_interactions_br <- apply(inter_mat_br, c(1, 2), mean)
print(dim(mean_interactions)) 
mean_interactions_br_df <- as.data.frame(mean_interactions_br)
# keeping row names as separate column for clarity
mean_interactions_br_df$RowNames <- rownames(mean_interactions_br)
# reorder to rownames is 1st col
mean_interactions_br_df <- mean_interactions_br_df[, c("RowNames", setdiff(colnames(mean_interactions_br_df), "RowNames"))]

write.csv(mean_interactions_br_df, "mean_interaction_matrix_br.csv", row.names = FALSE)


# trying out saving the beta_ij values from br's fit to plug into other models as initializations #########
# looking at parameter names
# str(joint.post.draws_br) 
# # Extract beta_ij samples from the first model
# beta_ij_samples_br <- joint.post.draws_br$ndd_betaij
# beta_ij_init_br <- apply(beta_ij_samples_br, c(2, 3), median)  # Median for each (S, T) pair
# # need to match the sem data (52, > stan.data_sem$I) with the 64 from br$I
# beta_ij_init_sem <- beta_ij_init_br[cbind(stan.data_sem$irow, stan.data_sem$icol)]
# #length(beta_ij_init_sem)
# #stan.data_br$I
# # creating a function that can be plugged into the next models
# init_function <- function() {
#   list(
#     beta_ij = beta_ij_init_sem
#   )
# }


# figuring out how to get the posterior distributions (getting posteriors of ndd_betaij) ####
# getting posterior distributions (http://mc-stan.org/bayesplot/)
# posterior <- as.matrix(fit_br)
# # # Check parameter names in the posterior object
# str(posterior)
# colnames(posterior)
# # # Extract all ndd_betaij values
# beta_ij_values <- as.data.frame(posterior[, grep("ndd_betaij", colnames(posterior))])
# # # Create a faceted density plot
# summary(beta_ij_values)
# nrow(posterior)
# subset_beta_ij <- beta_ij_values[, 1:20]
# mcmc_dens(subset_beta_ij, facet_args = list(ncol = 5), adjust = 2)



#length(beta_ij_init_br) == stan.data_sem$I 
#str(init_list)

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
               #init = init_function
               )


pairs(fit_sem, pars = c("response[1]", "effect[1]"))

print(summary(fit_sem, pars=c("gamma_i","ndd_betaij","ri_betaij"))$summary)
rstan::traceplot(fit_sem, pars=c("response","effect"))
rstan::stan_rhat(fit_sem, pars=c("effect"))


# Get the full posteriors ###### 
joint.post.draws_sem <- extract.samples(fit_sem)

## ---  Get interaction estimates ###
inter_mat_sem <- aperm(joint.post.draws_sem$ndd_betaij, c(2, 3, 1))
rownames(inter_mat_sem) <- focalID_sem
colnames(inter_mat_sem) <- neighbourID_sem

## getting the mean of all the posteriors  #####
mean_interactions_sem <- apply(inter_mat_sem, c(1, 2), mean)
# want to save this to create the visualizations
mean_interactions_sem_df <- as.data.frame(mean_interactions_sem)
# keeping row names as separate column for clarity
mean_interactions_sem_df$RowNames <- rownames(mean_interactions_sem)
# reorder to rownames is 1st col
mean_interactions_sem_df <- mean_interactions_sem_df[, c("RowNames", setdiff(colnames(mean_interactions_sem_df), "RowNames"))]

write.csv(mean_interactions_sem_df, "mean_interaction_matrix_sem.csv", row.names = FALSE)



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

pairs(fit_rf, pars = c("response[1]", "effect[1]"))

print(summary(fit_rf, pars=c("gamma_i","ndd_betaij","ri_betaij"))$summary)
rstan::traceplot(fit_rf, pars=c("response","effect"))
rstan::stan_rhat(fit_rf, pars=c("effect"))


# Get the full posteriors ###### 
joint.post.draws_rf <- extract.samples(fit_rf)

## ---  Get interaction estimates ###
inter_mat_rf <- aperm(joint.post.draws_rf$ndd_betaij, c(2, 3, 1))
rownames(inter_mat_rf) <- focalID_rf
colnames(inter_mat_rf) <- neighbourID_rf

## getting the mean of all the posteriors  #####
mean_interactions_rf <- apply(inter_mat_rf, c(1, 2), mean)
# want to save this to create the visualizations
mean_interactions_rf_df <- as.data.frame(mean_interactions_rf)
# keeping row names as separate column for clarity
mean_interactions_rf_df$RowNames <- rownames(mean_interactions_rf)
# reorder to rownames is 1st col
mean_interactions_rf_df <- mean_interactions_rf_df[, c("RowNames", setdiff(colnames(mean_interactions_rf_df), "RowNames"))]

write.csv(mean_interactions_rf_df, "mean_interaction_matrix_rf.csv", row.names = FALSE)


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


pairs(fit_wir, pars = c("response[1]", "effect[1]"))

print(summary(fit_wir, pars=c("gamma_i","ndd_betaij","ri_betaij"))$summary)
rstan::traceplot(fit_wir, pars=c("response","effect"))
rstan::stan_rhat(fit_wir, pars=c("effect"))


# Get the full posteriors ###### 
joint.post.draws_wir <- extract.samples(fit_wir)

## ---  Get interaction estimates ###
inter_mat_wir <- aperm(joint.post.draws_wir$ndd_betaij, c(2, 3, 1))
rownames(inter_mat_wir) <- focalID_wir
colnames(inter_mat_wir) <- neighbourID_wir

## getting the mean of all the posteriors  #####
mean_interactions_wir <- apply(inter_mat_wir, c(1, 2), mean)
# want to save this to create the visualizations
mean_interactions_wir_df <- as.data.frame(mean_interactions_wir)
# keeping row names as separate column for clarity
mean_interactions_wir_df$RowNames <- rownames(mean_interactions_wir)
# reorder to rownames is 1st col
mean_interactions_wir_df <- mean_interactions_wir_df[, c("RowNames", setdiff(colnames(mean_interactions_wir_df), "RowNames"))]

write.csv(mean_interactions_wir_df, "mean_interaction_matrix_wir.csv", row.names = FALSE)






# ORIGINAL CODE
# Draw 1000 samples from the 80% posterior interval for each parameter of interest
# p.samples <- list()
# p.samples <- sapply(param.vec[param.vec %in% c('ndd_betaij')], function(p) {
#   p.samples[[p]] <- apply(joint.post.draws[[p]], 2, function(x){
#     sample(x[x > quantile(x, 0.1) & x < quantile(x, 0.9)], size = 100)
#   })  # this only works for parameters which are vectors
# })


# ORIGINAL CODE
# WARNING: in the STAN model for annual wildflowers, parameter 'gamma_i' lies within an exponential,
# 'gamma_i' estimates must thus be exponentiated to return estimates of intrinsic performance
# intrinsic.perf <- exp(p.samples$gamma_i)
# colnames(intrinsic.perf) <- focalID



# inter_mat is now a 3 dimensional array, where rows = focals, columns = neighbours and 3rd dim = samples from the posterior
# inter_mat[ , , 1] should return a matrix consisting of one sample for every interaction 
# apply(inter_mat, c(1, 2), mean) will return the mean estimate for every interaction (NB: this is the 
# mean of the 80% posterior interval, so will be slightly different to the mean value returned from 
# summary(fit), which is calculated from the full posterior distribution) 

# Interactions can now be divided by the appropriate scaling (intrinsic performance, and demographic rates
# if a population dynamics model is used) in order to return per capita interaction strengths. 


