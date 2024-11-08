# Running the joint model on simulated data 

# set up R environment ----
require(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()) 
#rstan_options(auto_write = TRUE)
library(tidybayes)
library(bayesplot)
library(tidyverse)
library(rethinking)
#install.packages("reshape2")
library(reshape2)

# load required functions
source('data_prep.R')
# instead of the simul_data, want our cleaned dataframe
# source('simul_data.R')  (don't need for running own data)

# instead of simdat, usingobserved values 

# all of this below is for the sem data ()
# this is the more simple version of the code, is only the diverse values
  # ! TODO: double check how I created this dataframe ----

# load raw data----
sem_data <- read.csv("../../Parks data cleaning/only_diverse_sem.csv")
head(sem_data)

#start by trying only with the focal species (add in other non-focals later once this one is working)
df <- sem_data %>% dplyr::select(focal, PLOT, seeds, CLAPUR, COLLIN, COLLOM, EPIDEN, GILCAP, NAVSQU, PLAFIG, PLECON)
head(df)

# identify focal and neighbouring species to be matched to parameter estimates
focalID <- unique(df$focal)  # this should return the names of unique focal groups in the order in which they are encountered in the dataframe - must be alphabetical

neighbourID <- colnames(df[ , -c(1:3)]) # should be ordered focal first (alphabetically), then non-focals in alphabetical order

# ensure neighbours are linearly independent across the whole dataset (see S1.2)
N_all <- df[ , neighbourID]
N_all <- apply(N_all, c(1,2), as.numeric)
X_all <- cbind(model.matrix(~as.factor(df$focal)), N_all)
anyNA(X_all)
any(is.nan(X_all))
# replacing the na values with 0!! 
X_all[is.na(X_all)] <- 0
# double checking for na values again
anyNA(X_all)

R_all <- pracma::rref(X_all)
Z_all <- t(R_all) %*% R_all
indep <- sapply(seq(1, dim(Z_all)[1], 1), function(k){ 
  ifelse(Z_all[k, k] == 1 & sum(Z_all[k, -k]) == 0, 1, 0)
}) #
all(indep == 1) # if TRUE then neighbours are linearly independent and we can continue
if(!all(indep == 1)) warning('WARNING neighbours are not linearly independent') 


# prepare the data into the format required by STAN and the model code
stan.data <- data_prep(perform = 'seeds', 
                       focal = 'focal', 
                       nonNcols = 3, # number of columns that aren't neighbour abundances
                       df = df)

# now also have to force the doubles in perform to be ints
stan.data$perform <- as.integer(round(stan.data$perform))

message(paste0('Data dimensions = ', dim(df)[1], ', ', dim(df)[2]))
message(paste0('Number of focal groups = ', length(focalID)))
message(paste0('Number of neighbour groups = ', length(neighbourID)))
message(paste0('Proportion of inferrable interactions = ', sum(stan.data$Q)/(stan.data$S*stan.data$`T`)))


# Run the model ----
stan.seed <- 1234
fit <- stan(file = 'joint_model.stan', 
            data =  stan.data,               # named list of data
            chains = 3,
            cores = 3, # added so runs in parallel
            warmup = 2000,          # number of warmup iterations per chain
            iter = 5000,            # total number of iterations per chain
            refresh = 100,         # show progress every 'refresh' iterations
            control = list(max_treedepth = 19,
                           adapt_delta = 0.99), 
            seed = stan.seed
)

#saveRDS(fit, file = "../../Parks interaction SHARED/saved_models/fit.RDS")

# check convergence
print(summary(fit, pars=c("gamma_i","ndd_betaij",""))$summary)
rstan::traceplot(fit, pars=c("gamma_i","ndd_betaij"))
rstan::stan_rhat(fit, pars=c("gamma_i","ndd_betaij"))

# save traceplots 
#pdf(file="../../Parks interaction SHARED/saved_models/fit_traceplots.pdf", width = 12, height = 12)
print(mcmc_trace(fit, regex_pars = "ndd_betaij") + xlab("Post-warmup iteration"))
print(mcmc_trace(fit, regex_pars = "ri_betaij") + xlab("Post-warmup iteration"))
print(mcmc_trace(fit, regex_pars = "gamma_i") + xlab("Post-warmup iteration"))
dev.off()


# Get the full posteriors 


# using the rethinking package

joint.post.draws <- extract.samples(fit)

# Select parameters of interest
param.vec <- fit@model_pars[!fit@model_pars %in% c('lp__')]
# Draw 1000 samples from the 80% posterior interval for each parameter of interest
p.samples <- list()

# this sapply gives this error 'MARGIN' does not match dim(X) -----
p.samples <- sapply(param.vec[!param.vec %in% c('ri_betaij', 'ndd_betaij')], function(p) {
  p.samples[[p]] <- apply(joint.post.draws[[p]], 2, function(x){
    sample(x[x > quantile(x, 0.1) & x < quantile(x, 0.9)], size = 100)
  })  # this only works for parameters which are vectors
})

# alternate function for the one above (confused on what exactly this is doing/check) ----
p.samples <- sapply(param.vec[!param.vec %in% c('ri_betaij', 'ndd_betaij')], function(p) {
  
  cat("Inspecting parameter:", p, "\n")
  print(str(joint.post.draws[[p]]))  # Check the structure of the current element
  
  # If it's a matrix, use apply
  if (is.matrix(joint.post.draws[[p]])) {
    return(apply(joint.post.draws[[p]], 2, function(x){
      sample(x[x > quantile(x, 0.1) & x < quantile(x, 0.9)], size = 100)
    }))
    
    # If it's a vector, sample directly
  } else if (is.vector(joint.post.draws[[p]])) {
    return(sample(joint.post.draws[[p]][joint.post.draws[[p]] > quantile(joint.post.draws[[p]], 0.1) & joint.post.draws[[p]] < quantile(joint.post.draws[[p]], 0.9)], size = 100))
    
    # For unsupported structures
  } else {
    cat("Unsupported structure for parameter:", p, "\n")
    return(NULL)  # Return NULL for unsupported cases
  }
})

# WARNING: in the STAN model for annual wildflowers, parameter 'gamma_i' lies within an exponential,
# 'gamma_i' estimates must thus be exponentiated to return estimates of intrinsic performance
intrinsic.perf <- exp(p.samples$gamma_i)
colnames(intrinsic.perf) <- focalID

# Get interaction estimates
#--------------------------
inter_mat <- aperm(joint.post.draws$ndd_betaij, c(2, 3, 1))
rownames(inter_mat) <- focalID
colnames(inter_mat) <- neighbourID



## instead of using the rethinking package ----
n.draws <- dim(as.matrix(fit))[[1]]
rand.draws <- sample(1:n.draws, 1000)
# # Get 
get_variables(fit)

joint.post.draws <- fit %>%
  spread_draws(ndd_betaij[i,j]) %>%
  filter(.draw %in% rand.draws)%>%
  mutate(param = "beta")
head(joint.post.draws); dim(joint.post.draws)
with(joint.post.draws, table(i,j)) # check

# so just need to wrangle joint.post.draws into array format if you want to get it in same output form as what Malyon is using below, from the sound of it...



# Get interaction estimates
#--------------------------
inter_mat <- aperm(joint.post.draws$ndd_betaij, c(2, 3, 1))
rownames(inter_mat) <- focalID
colnames(inter_mat) <- neighbourID

# inter_mat is now a 3 dimensional array, where rows = focals, columns = neighbours and 3rd dim = samples from the posterior
# inter_mat[ , , 1] should return a matrix consisting of one sample for every interaction 
# apply(inter_mat, c(1, 2), mean) will return the mean estimate for every interaction (NB: this is the 
# mean of the 80% posterior interval, so will be slightly different to the mean value returned from 
# summary(fit), which is calculated from the full posterior distribution) 

# Interactions can now be divided by the appropriate scaling (intrinsic performance, and demographic rates
# if a population dynamics model is used) in order to return per capita interaction strengths. 

# general code on some ways to look at outputs

# inter_sample_1 <- inter_mat[, , 1]  # Extract the first sample
# print(inter_sample_1)               # View the first interaction matrix
# 
# mean_inter_mat <- apply(inter_mat, c(1, 2), mean)
# print(mean_inter_mat)
# heatmap(mean_inter_mat, Rowv = NA, Colv = NA, scale = "none", 
#         main = "Mean Interaction Matrix", xlab = "Neighbors", ylab = "Focals")





