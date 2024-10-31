# Running the joint model on simulated data 

# set up R environment
require(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()) 
#rstan_options(auto_write = TRUE)

library(rethinking)
#install.packages("reshape2")
library(reshape2)

# load required functions
# make sure ur in the correct working directory
source('data_prep.R')
# instead of the simul_data, want our cleaned dataframe
source('simul_data.R')
#source("~/Documents/capstone/Parks cleaning data/sem_neighbor_focal_alone_diverse_merge.csv")


# load simulated data 
set.seed(54)

# instead of simdat, using list(observed values df, sim_a/intrinsic performance data,
# sim_true_alpha or interaction strengths 

#simdat <- simul_data(S=3, T=3, p=0.4)   # p = proportion of interactions which are NOT observed 


# all of this below is for the sem data ()

#sem_data <- read.csv("~/Documents/capstone/Parks data cleaning/sem_neighbor_focal_alone_diverse_merge.csv")

# this is the more simple version of the code, is only the diverse values
  # TODO: double check how I created this dataframe
sem_data <- read.csv("~/Documents/capstone/Parks data cleaning/bimler_sem.csv")

# 
# # getting the alpha values/intrinsic performance data 
avg_seeds <- aggregate(seeds ~ focal, data = sem_data, mean)
sim_a <- log(avg_seeds$seeds)

# 
# # using a regression model to find the simulated true alphas, the interaction strengths
# 
# # Example using a linear model (might need a more complex model?)
interaction_model <- lm(seeds ~ APOAND + CERGLO + CIRSIUM + CLAPUR + COLLIN + COLLOM + DAUCAR + EPIDEN + ERILAN + GALIUM + GERANIUM + GILCAP + HYPERNICUM + LACTUCA + LATHVIC + LEUVUL + LUPBIC + MYOSOTIS + NAVSQU + PLAFIG + PLECON + SHEARV + TRIFOLIUM + VERARV, data = sem_data)
sim_truealpha <- coef(interaction_model)[-1]  # Exclude the intercept

# converting the sim_truealpha into a matrix like they had
# Assuming you have S species and T neighbor species
S <- length(unique(sem_data$focal))
# getting the number of unique values
# TODO: confused about the simulated values?
T <- length(sim_truealpha)/ S  # Adjust this based on your model
sim_truealpha_matrix <- matrix(sim_truealpha, nrow = S, ncol = T)
# 
simdat <- list(sem_data, sim_a, sim_truealpha_matrix)





df <- simdat[[1]]
sim_gamma <- simdat[[2]]
sim_interactions <- simdat[[3]]
# NB: if using real data or named species, ensure they are ordered alphabetically in the dataset

# identify focal and neighbouring species to be matched to parameter estimates
focalID <- unique(df$focal)  # this should return the names of unique focal groups in the order
# in which they are encountered in the dataframe - must be alphabetical
neighbourID <- colnames(df[ , -c(1:2)]) # should be ordered focal first (alphabetically), then
# non-focals in alphabetical order

# ensure neighbours are linearly independent across the whole dataset (see S1.2)
N_all <- df[ , neighbourID]
N_all <- apply(N_all, c(1,2), as.numeric)
X_all <- cbind(model.matrix(~as.factor(df$focal)), N_all)
# note about using the bimler_sem.csv, it seems like when it's trying to get these into T/F there are some NA values, messing it up 
  #  checking for na values, na values for some of the focals at sem 
  # 
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


# Notes: 10/29: changed the nonNcols from 2 to 3 (bc we have the X: 1,2,3 etc column, focal, and seeds column)
# prepare the data into the format required by STAN and the model code
stan.data <- data_prep(perform = 'seeds', 
                       focal = 'focal', 
                       nonNcols = 3, # number of columns that aren't neighbour abundances
                       df = df)


message(paste0('Data dimensions = ', dim(df)[1], ', ', dim(df)[2]))
message(paste0('Number of focal groups = ', length(focalID)))
message(paste0('Number of neighbour groups = ', length(neighbourID)))
message(paste0('Proportion of inferrable interactions = ', sum(stan.data$Q)/(stan.data$S*stan.data$`T`)))

# checking again for na values bc the fit/model is giving an error about na values
sapply(stan.data, anyNA)
# 10/29 I don't know exactly why this happening
  # TODO: look more into this/data cleaning
# # replacing the na values with 0!! 
# X_all[is.na(X_all)] <- 0
stan.data$perform[is.na(stan.data$perform)] <- 0 
# now also have to force the doubles in perform to be ints
stan.data$perform <- as.integer(round(stan.data$perform))

# Run the model! 
stan.seed <- 1234
#library(V8)
#install.packages("V8")
fit <- stan(file = 'joint_model.stan', 
            data =  stan.data,               # named list of data
            chains = 4,
            warmup = 2000,          # number of warmup iterations per chain
            iter = 3000,            # total number of iterations per chain
            refresh = 100,         # show progress every 'refresh' iterations
            control = list(max_treedepth = 10,
                           adapt_delta = 0.9), 
            seed = stan.seed
)

# check convergence
print(summary(fit, pars=c("gamma_i","ndd_betaij","ri_betaij"))$summary)
rstan::traceplot(fit, pars=c("gamma_i","ndd_betaij"))
rstan::stan_rhat(fit, pars=c("gamma_i","ndd_betaij"))

# Get the full posteriors 
# joint.post.draws <- extract.samples(fit)
# instead of using the rethinking package
n.draws <- dim(as.matrix(fit))

# Select parameters of interest
param.vec <- fit@model_pars[!fit@model_pars %in% c('lp__')]

# Draw 1000 samples from the 80% posterior interval for each parameter of interest
p.samples <- list()


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


# 
# p.samples <- sapply(param.vec[!param.vec %in% c('ri_betaij', 'ndd_betaij')], function(p) {
#   p.samples[[p]] <- apply(joint.post.draws[[p]], 2, function(x){
#     sample(x[x > quantile(x, 0.1) & x < quantile(x, 0.9)], size = 100)
#   })  # this only works for parameters which are vectors
# })



# WARNING: in the STAN model for annual wildflowers, parameter 'gamma_i' lies within an exponential,
# 'gamma_i' estimates must thus be exponentiated to return estimates of intrinsic performance
intrinsic.perf <- exp(p.samples$gamma_i)
colnames(intrinsic.perf) <- focalID

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




