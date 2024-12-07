# Running the joint model on simulated data 

# set up R environment
require(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()) 
#rstan_options(auto_write = TRUE)
theme_set(theme_minimal())

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
library(patchwork)
library(kableExtra)

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

head(sem_df)
head(br_df)

# Deleted:
# Code for getting rid of the max seed values... ####

## Plot fitness data ----
br_plot <- ggplot(br_df, aes(x = seeds)) +
  geom_histogram() +
  facet_wrap(~species, scales='free')+
  ggtitle("Distribution of Seeds BR")
# it's just NAVSQU with super igh seed values

# sem
ggplot(sem_df, aes(x = seeds)) +
  geom_histogram() +
  facet_wrap(~species, scales='free')+
  ggtitle("Distribution of Seeds SEM")

# wir 
ggplot(wir_df, aes(x = seeds)) +
  geom_histogram() +
  facet_wrap(~species, scales='free')+
  ggtitle("Distribution of Seeds WIR")

# rf 
ggplot(rf_df, aes(x = seeds)) +
  geom_histogram() +
  facet_wrap(~species, scales='free')+
  ggtitle("Distribution of Seeds RF")




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

# BR ----

## Plot raw curves ----
head(temp)

# get intra and inter specific densities to plot
species_map <- sort(unique(br_df$species))
br.raw <- data.frame(seeds = stan.data_br$perform, species=stan.data_br$species_ID,stan.data_br$X) %>%
  mutate(species_name = species_map[species]) %>%
  rowwise() %>%
  mutate(inter=sum(across(CLAPUR:PLECON)[names(across(CLAPUR:PLECON)) != species_name]),
         intra=sum(across(CLAPUR:PLECON)[names(across(CLAPUR:PLECON)) == species_name]),
         total=sum(across(CLAPUR:PLECON)))%>%
  ungroup()
head(br.raw)

head(filter(br.raw, total==0))
  
# Plot raw data as a check
ggplot(br.raw, aes(x=total, y=seeds)) +          
  geom_jitter(width=.1) + 
  facet_wrap(~ species_name, scales = "free") +
  stat_smooth(method = "nls",
              formula = y ~ a/(1+b*x),
              method.args = list(start = list(a = 1000, b = .01)),
              se = FALSE) +
  labs(x='number of competitors', y='seeds produced', title='Blanton Ridge')
ggsave("../../../Parks interaction SHARED/Figures/BR_raw.pdf", width=9, height=6)

## Run BR model ----

fit_br <- stan(file = 'joint_model.stan', 
            data =  stan.data_br,               # named list of data
            chains = 3,
            cores = 3, # added so runs in parallel
            # warmup = 2000,          # number of warmup iterations per chain
            iter = 5000,            # total number of iterations per chain
            # refresh = 100,         # show progress every 'refresh' iterations
            control = list(max_treedepth = 19,
                           adapt_delta = 0.99), 
            seed = stan.seed
)

fit_br  # converged well!
# extracting the values for br bc we know that the model likes this data
br_fit <- rstan::extract(fit_br)

# save model object
saveRDS(fit_br, file = paste("../../../Parks interaction SHARED/saved_models/BR_Stan_fit.RDS", sep=''))

# save summary
temp <-summary(fit_br)$summary
save_kable(kbl(temp, digits=3), file=paste("../../../Parks interaction SHARED/saved_models/BR.fit.summary.pdf",sep=''))

## save traceplots ----
# pdf(file="../../Parks interaction SHARED/saved_models/BR_fit_traceplots.pdf", width = 12, height = 12)
pdf(file="../../../Parks interaction SHARED/saved_models/BR_fit_traceplots.pdf", width = 12, height = 12)
print(mcmc_trace(fit_br, regex_pars = "ndd_betaij") + xlab("Post-warmup iteration"))
print(mcmc_trace(fit_br, regex_pars = "ri_betaij") + xlab("Post-warmup iteration"))
print(mcmc_trace(fit_br, regex_pars = "gamma_i") + xlab("Post-warmup iteration"))
print(mcmc_trace(fit_br, regex_pars = "effect") + xlab("Post-warmup iteration"))
print(mcmc_trace(fit_br, regex_pars = "response") + xlab("Post-warmup iteration"))
dev.off()


# Get the full posteriors 
joint.post.draws_br <- extract.samples(fit_br)

## ---  Get interaction estimates ###
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


# figuring out how to get the posterior distributions (idk which thing to get the posteriors of) 
# getting posterior distributions (http://mc-stan.org/bayesplot/)
# posterior <- as.matrix(fit)
# # Check parameter names in the posterior object
# str(posterior)
# # Extract all log_lik_nddm values
# log_lik_nddm_values <- as.data.frame(posterior[, grep("log_lik_nddm", colnames(posterior))])
# # Create a faceted density plot
# mcmc_dens(log_lik_nddm_long, facet_args = list(ncol = 5))



# SEM ----

## Plot raw curves ----
# get intra and inter specific densities to plot
species_map <- sort(unique(sem_df$species))
sem.raw <- data.frame(seeds = stan.data_sem$perform, species=stan.data_sem$species_ID,stan.data_sem$X) %>%
  mutate(species_name = species_map[species]) %>%
  rowwise() %>%
  mutate(inter=sum(across(CLAPUR:PLECON)[names(across(CLAPUR:PLECON)) != species_name]),
         intra=sum(across(CLAPUR:PLECON)[names(across(CLAPUR:PLECON)) == species_name]),
         total=sum(across(CLAPUR:PLECON)))%>%
  ungroup()
head(sem.raw)

head(filter(sem.raw, total==0))

# Plot raw data as a check
ggplot(sem.raw, aes(x=total, y=seeds)) +          
  geom_jitter(width=.1) + 
  facet_wrap(~ species_name, scales = "free") +
  stat_smooth(method = "nls",
              formula = y ~ a/(1+b*x),
              method.args = list(start = list(a = 1000, b = .01)),
              se = FALSE) +
  labs(x='number of competitors', y='seeds produced', title='SEM')
ggsave("../../../Parks interaction SHARED/Figures/sem_raw.pdf", width=9, height=6)

## Run sem model ----

fit_sem <- stan(file = 'joint_model.stan', 
               data =  stan.data_sem,               # named list of data
               chains = 3,
               cores = 3, # added so runs in parallel
               # warmup = 2000,          # number of warmup iterations per chain
               iter = 5000,            # total number of iterations per chain
               # refresh = 100,         # show progress every 'refresh' iterations
               control = list(max_treedepth = 19,
                              adapt_delta = 0.99), 
               seed = stan.seed
)

fit_sem  #
# extracting the values for sem bc we know that the model likes this data
# sem_fit <- rstan::extract(fit_sem)

# save summary
temp <- summary(fit_sem)$summary
save_kable(kbl(temp, digits=3), file=paste("../../../Parks interaction SHARED/saved_models/SEM.fit.summary.pdf",sep=''))

## save traceplots ----
pdf(file="../../../Parks interaction SHARED/saved_models/sem_fit_traceplots.pdf", width = 12, height = 12)
print(mcmc_trace(fit_sem, regex_pars = "ndd_betaij") + xlab("Post-warmup iteration"))
print(mcmc_trace(fit_sem, regex_pars = "ri_betaij") + xlab("Post-warmup iteration"))
print(mcmc_trace(fit_sem, regex_pars = "gamma_i") + xlab("Post-warmup iteration"))
print(mcmc_trace(fit_sem, regex_pars = "effect") + xlab("Post-warmup iteration"))
print(mcmc_trace(fit_sem, regex_pars = "response") + xlab("Post-warmup iteration"))
dev.off()

# save model object
fit_sem
saveRDS(fit_sem, file = paste("../../../Parks interaction SHARED/saved_models/SEM_Stan_fit.RDS", sep=''))



# WIR ----

## Plot raw curves ----
# get intra and inter specific densities to plot
species_map <- sort(unique(wir_df$species))
wir.raw <- data.frame(seeds = stan.data_wir$perform, species=stan.data_wir$species_ID,stan.data_wir$X) %>%
  mutate(species_name = species_map[species]) %>%
  rowwise() %>%
  mutate(inter=sum(across(CLAPUR:PLECON)[names(across(CLAPUR:PLECON)) != species_name]),
         intra=sum(across(CLAPUR:PLECON)[names(across(CLAPUR:PLECON)) == species_name]),
         total=sum(across(CLAPUR:PLECON)))%>%
  ungroup()
head(wir.raw)

head(filter(wir.raw, total==0))

# Plot raw data as a check
ggplot(wir.raw, aes(x=total, y=seeds)) +          
  geom_jitter(width=.1) + 
  facet_wrap(~ species_name, scales = "free") +
  stat_smooth(method = "nls",
              formula = y ~ a/(1+b*x),
              method.args = list(start = list(a = 1000, b = .01)),
              se = FALSE) +
  labs(x='number of competitors', y='seeds produced', title='WIR')
ggsave("../../../Parks interaction SHARED/Figures/wir_raw.pdf", width=9, height=6)

## Run wir model ----

fit_wir <- stan(file = 'joint_model.stan', 
                data =  stan.data_wir,               # named list of data
                chains = 3,
                cores = 3, # added so runs in parallel
                # warmup = 2000,          # number of warmup iterations per chain
                iter = 5000,            # total number of iterations per chain
                # refresh = 100,         # show progress every 'refresh' iterations
                control = list(max_treedepth = 19,
                               adapt_delta = 0.99), 
                seed = stan.seed
)

fit_wir
# extracting the values for wir bc we know that the model likes this data
# wir_fit <- rstan::extract(fit_wir)

## save traceplots ----
pdf(file="../../../Parks interaction SHARED/saved_models/wir_fit_traceplots.pdf", width = 12, height = 12)
print(mcmc_trace(fit_wir, regex_pars = "ndd_betaij") + xlab("Post-warmup iteration"))
print(mcmc_trace(fit_wir, regex_pars = "ri_betaij") + xlab("Post-warmup iteration"))
print(mcmc_trace(fit_wir, regex_pars = "gamma_i") + xlab("Post-warmup iteration"))
print(mcmc_trace(fit_wir, regex_pars = "effect") + xlab("Post-warmup iteration"))
print(mcmc_trace(fit_wir, regex_pars = "response") + xlab("Post-warmup iteration"))
dev.off()

# save summary
temp <-summary(fit_wir)$summary
save_kable(kbl(temp, digits=3), file=paste("../../../Parks interaction SHARED/saved_models/WIR.fit.summary.pdf",sep=''))

# save model object
fit_wir
saveRDS(fit_wir, file = paste("../../../Parks interaction SHARED/saved_models/WIR_Stan_fit.RDS", sep=''))


# RF ----

## Plot raw curves ----
# get intra and inter specific densities to plot
species_map <- sort(unique(rf_df$species))
rf.raw <- data.frame(seeds = stan.data_rf$perform, species=stan.data_rf$species_ID,stan.data_rf$X) %>%
  mutate(species_name = species_map[species]) %>%
  rowwise() %>%
  mutate(inter=sum(across(CLAPUR:PLECON)[names(across(CLAPUR:PLECON)) != species_name]),
         intra=sum(across(CLAPUR:PLECON)[names(across(CLAPUR:PLECON)) == species_name]),
         total=sum(across(CLAPUR:PLECON)))%>%
  ungroup()
head(rf.raw)

head(filter(rf.raw, total==0))

# Plot raw data as a check
ggplot(rf.raw, aes(x=total, y=seeds)) +          
  geom_jitter(width=.1) + 
  facet_wrap(~ species_name, scales = "free") +
  stat_smooth(method = "nls",
              formula = y ~ a/(1+b*x),
              method.args = list(start = list(a = 1000, b = .01)),
              se = FALSE) +
  labs(x='number of competitors', y='seeds produced', title='RF')
ggsave("../../../Parks interaction SHARED/Figures/RF_raw.pdf", width=9, height=6)

## Run rf model ----

fit_rf <- stan(file = 'joint_model.stan', 
                data =  stan.data_rf,               # named list of data
                chains = 3,
                cores = 3, # added so runs in parallel
                # warmup = 2000,          # number of warmup iterations per chain
                iter = 5000,            # total number of iterations per chain
                # refresh = 100,         # show progress every 'refresh' iterations
                control = list(max_treedepth = 19,
                               adapt_delta = 0.99)
                # seed = stan.seed
)

fit_rf
# extracting the values for rf bc we know that the model likes this data
# rf_fit <- rstan::extract(fit_rf)

# save summary
temp <-summary(fit_rf)$summary
save_kable(kbl(temp, digits=3), file=paste("../../../Parks interaction SHARED/saved_models/RF.fit.summary.pdf",sep=''))

## save traceplots ----
pdf(file="../../../Parks interaction SHARED/saved_models/rf_fit_traceplots.pdf", width = 12, height = 12)
print(mcmc_trace(fit_rf, regex_pars = "ndd_betaij") + xlab("Post-warmup iteration"))
print(mcmc_trace(fit_rf, regex_pars = "ri_betaij") + xlab("Post-warmup iteration"))
print(mcmc_trace(fit_rf, regex_pars = "gamma_i") + xlab("Post-warmup iteration"))
print(mcmc_trace(fit_br, regex_pars = "effect") + xlab("Post-warmup iteration"))
print(mcmc_trace(fit_br, regex_pars = "response") + xlab("Post-warmup iteration"))
dev.off()

# save model object
saveRDS(fit_rf, file = paste("../../../Parks interaction SHARED/saved_models/RF_Stan_fit.RDS", sep=''))


