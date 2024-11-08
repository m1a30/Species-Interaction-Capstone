
  
sem_data <- read.csv("~/Documents/capstone/Parks data cleaning/sem_neighbor_focal_alone_diverse_merge.csv")


# getting the alpha values/intrinsic performance data 
avg_seeds <- aggregate(seeds ~ focal, data = sem_data, mean)
sim_a <- log(avg_seeds$seeds)

# using a regression model to find the simulated true alphas, the interaction strengths

# Example using a linear model (might need a more complex model?)
interaction_model <- lm(seeds ~ APOAND + CERGLO + CIRSIUM + CLAPUR + COLLIN + COLLOM + DAUCAR + EPIDEN + ERILAN + GALIUM + GERANIUM + GILCAP + HYPERNICUM + LACTUCA + LATHVIC + LEUVUL + LUPBIC + MYOSOTIS + NAVSQU + PLAFIG + PLECON + SHEARV + TRIFOLIUM + VERARV, data = sem_data)
sim_truealpha <- coef(interaction_model)[-1]  # Exclude the intercept

# converting the sim_truealpha into a matrix like they had
# Assuming you have S species and T neighbor species
S <- length(unique(sem_data$focal))
# getting the number of 
  # THIS DOESNT MAKE SENSE RN BUT IM GONNA COME BACK TO IT 
T <- length(sim_truealpha)/ S  # Adjust this based on your model
sim_truealpha_matrix <- matrix(sim_truealpha, nrow = S, ncol = T)

final_data <- list(sem_data,sim_a,sim_truealpha_matrix)



