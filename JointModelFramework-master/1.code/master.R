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
#install.packages("reshape2")
library(reshape2)
library(tidyverse)
library(ggplot2)

# load required functions
source('data_prep.R')
source('simul_data.R')

# load simulated data 
set.seed(54)


#sem data
df <- read.csv("../../Parks data cleaning/data_cleaned/wir_dat.csv")

# seeing how different the output is when getting rid of large seed counts
# Identify the 10 largest values in the seeds column
top_10_max <- df %>%
  arrange(desc(seeds)) %>%
  slice_head(n = 10)

# Filter out the rows with the 10 largest values
df <- df %>%
  filter(!seeds %in% top_10_max$seeds)

ggplot(df, aes(x = seeds)) +
  geom_histogram(binwidth = 10) +
  ggtitle("Distribution of Seeds")

#start by trying only with the focal species (add in other non-focals later once this one is working)
df <- df %>% dplyr::select(species, plot, seeds, CLAPUR, COLLIN, COLLOM, EPIDEN, GILCAP, NAVSQU, PLAFIG, PLECON)

# modifying data frame so it's just plot species seeds and all neighbors
# df <- df %>%
#   select(2:3, 14, 17:ncol(df))

# keeping plot species seeds as is, then alphabetizing the neighbors
# df<- df %>%
#   select(1:3, sort(names(df)[-(1:3)]))

# TODO: need the seeds to be ints and not doubles? I think that's what's causing the error
df$seeds <- as.integer(df$seeds)

#filling the alone neighbors that are NAs to 0s 
df[is.na(df)] <- 0

# NB: if using real data or named species, ensure they are ordered alphabetically in the dataset

# 11/24 Changing df$focal to df$species to match our dataframe
# identify focal and neighbouring species to be matched to parameter estimates

focalID <- unique(df$species)  # this should return the names of unique focal groups in the order
# TODO: extra data cleaning step, need the focals to be alphabetized 
focalID <- sort(focalID)
# in which they are encountered in the dataframe - must be alphabetical

# 11/10: changed the -c(1:2) to -c(1:3), check df for how many non neighbor cols
neighbourID <- colnames(df[ , -c(1:3)]) # should be ordered focal first (alphabetically), then
# non-focals in alphabetical order

# ensure neighbours are linearly independent across the whole dataset (see S1.2)
N_all <- df[ , neighbourID]
N_all <- apply(N_all, c(1,2), as.numeric)
X_all <- cbind(model.matrix(~as.factor(df$species)), N_all)
R_all <- pracma::rref(X_all)
Z_all <- t(R_all) %*% R_all
indep <- sapply(seq(1, dim(Z_all)[1], 1), function(k){ 
  ifelse(Z_all[k, k] == 1 & sum(Z_all[k, -k]) == 0, 1, 0)
}) #
all(indep == 1) # if TRUE then neighbours are linearly independent and we can continue
if(!all(indep == 1)) warning('WARNING neighbours are not linearly independent') 


# Notes: 10/29: changed the nonNcols from 2 to 3 (bc we have the X: 1,2,3 etc column, focal, and seeds column)
# prepare the data into the format required by STAN and the model code
stan.data <- data_prep(perform = 'seeds', 
                       focal = 'species', 
                       nonNcols = 3, # number of columns that aren't neighbour abundances
                       df = df)


message(paste0('Data dimensions = ', dim(df)[1], ', ', dim(df)[2]))
message(paste0('Number of focal groups = ', length(focalID)))
message(paste0('Number of neighbour groups = ', length(neighbourID)))
message(paste0('Proportion of inferrable interactions = ', sum(stan.data$Q)/(stan.data$S*stan.data$`T`)))

# Run the model! 
stan.seed <- 1234
fit <- stan(file = 'joint_model.stan', 
            data =  stan.data,               # named list of data
            chains = 3,
            cores = 3, # added so runs in parallel
            warmup = 2000,          # number of warmup iterations per chain
            iter = 3000,            # total number of iterations per chain
            refresh = 100,         # show progress every 'refresh' iterations
            control = list(max_treedepth = 19,
                           adapt_delta = 0.99), 
            seed = stan.seed
)

## code to look at where the high rhat values are
# fit_summary <- as.data.frame(summary(fit)$summary)
# high_rhat <- fit_summary[fit_summary$Rhat > 1.1, ]
# low_ess <- fit_summary[fit_summary$n_eff < 100, ]
# View(low_ess)
# View(high_rhat)
# was finding that response and effect were high values 
#bayesplot::mcmc_pairs(fit, pars = c("response[1]", "response[2]", "effect[1]", "effect[2]"))


# check convergence
print(summary(fit, pars=c("gamma_i","ndd_betaij","ri_betaij"))$summary)
rstan::traceplot(fit, pars=c("gamma_i","ndd_betaij"))
rstan::stan_rhat(fit, pars=c("gamma_i","ndd_betaij"))


# figuring out how to get the posterior distributions (idk which thing to get the posteriors of) ####
# getting posterior distributions (http://mc-stan.org/bayesplot/)
# posterior <- as.matrix(fit)
# # Check parameter names in the posterior object
# str(posterior)
# # Extract all log_lik_nddm values
# log_lik_nddm_values <- as.data.frame(posterior[, grep("log_lik_nddm", colnames(posterior))])
# # Create a faceted density plot
# mcmc_dens(log_lik_nddm_long, facet_args = list(ncol = 5))



# Get the full posteriors 
joint.post.draws <- extract.samples(fit)

## ---  Get interaction estimates ###
inter_mat <- aperm(joint.post.draws$ndd_betaij, c(2, 3, 1))
rownames(inter_mat) <- focalID
colnames(inter_mat) <- neighbourID

## getting the mean of all the posteriors 
mean_interactions <- apply(inter_mat, c(1, 2), mean)
library(qgraph)

### trying out a visualization  ##
all_sp <- rep("lightblue", nrow(mean_interactions)) # Example: all nodes are lightblue

# plotting for competition 
competition_plot <- qgraph(
  mean_interactions,  # use the mean_interactions matrix
  layout = 'circle', 
  posCol = rgb(red = 0, green = 0, blue = 0, alpha = 0),  # Facilitation = transparent
  negCol = 'orange',  # Competition = orange
  color = all_sp,  # Node colors
  labels = rownames(mean_interactions),  # Use row names of the matrix for labels
  fade = TRUE,
  directed = TRUE,
  title = 'Competition', 
  title.cex = 1.5
)

# plotting for facilitation 
facilitation_plot <- qgraph(
  mean_interactions,  # use the mean_interactions matrix
  layout = 'circle', 
  posCol = 'royalblue4',  # Facilitation = blue
  negCol = rgb(red = 0, green = 0, blue = 0, alpha = 0),  # Competition = transparent
  color = all_sp,  # Node colors
  labels = rownames(mean_interactions),  # Use row names of the matrix for labels
  fade = TRUE,
  directed = TRUE,
  title = 'Facilitation', 
  title.cex = 1.4
)


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


