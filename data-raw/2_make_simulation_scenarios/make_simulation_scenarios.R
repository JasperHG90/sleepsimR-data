# Make simulation scenarios

# Load sleepsimR
library(sleepsimR)

# Make scenarios
scen <- generate_scenarios()

# Subset scenarios
unique(scen$zeta)
# Rank
library(dplyr)
scen <- scen %>%
  arrange(desc(n), desc(n_t), scenario_id)

# Get a subset of 32 scenarios each
scen_subs <- scen %>%
  group_by(scenario_id) %>%
  slice(1:48)

# Remove from other scenarios
scen <- scen %>%
  filter(!iteration_id %in% scen_subs$iteration_id)

# Test a model
library(jsonlite)
sim <- as.list(scen[29000,])
# Make dataset
data_simulated <- simulate_dataset(sim$n, sim$n_t, sim$zeta, sim$Q, sim$dsim_seed)
# To data frame
tdf <- data.frame(
  id = data_simulated$obs[,1],
  EEG_mean_beta = data_simulated$obs[,2],
  EOG_median_theta = data_simulated$obs[,3],
  EOG_min_beta = data_simulated$obs[,4]
)
states <- data_simulated$states[,2]
# Get summary statistics for each
em1 <- tapply(tdf[,-1]$EEG_mean_beta, states, mean)
em2 <- tapply(tdf[,-1]$EOG_median_theta, states, mean)
em3 <- tapply(tdf[,-1]$EOG_min_beta, states, mean)
hyp_priors <- list(
  as.vector(em1[sort.list(as.numeric(names(em1)))]),
  as.vector(em2[sort.list(as.numeric(names(em2)))]),
  as.vector(em3[sort.list(as.numeric(names(em3)))])
)
# Reshape start values
sim$start_gamma <- jsonlite::parse_json(sim$start_gamma)
sim$start_emiss <- jsonlite::parse_json(sim$start_emiss)
m <- sqrt(length(sim$start_gamma$tpm))
start_values <- list(
  matrix(unlist(sim$start_gamma$tpm), nrow=m, ncol=m,
         byrow = TRUE),
  matrix(unlist(sim$start_emiss$EEG_Fpz_Cz_mean_beta),
         ncol=2, byrow = TRUE),
  matrix(unlist(sim$start_emiss$EOG_median_theta),
         ncol=2, byrow = TRUE),
  matrix(unlist(sim$start_emiss$EOG_min_beta),
         ncol=2, byrow=TRUE)
)
# Run model
mod <- sleepsimR::run_mHMM(tdf, start_values = start_values, hyperprior_means = hyp_priors,
                           model_seed = sim$model_seed,mcmc_iterations=3250, mcmc_burn_in = 1250,
                           order_data = FALSE)

# Arrange subset
scen_subs <- scen_subs %>%
  arrange(n, n_t, zeta, scenario_id)

# Save scenarios
write.csv(scen, file=gzfile("make_scenarios/scenarios.csv.gz"))
write.csv(scen_subs, file=gzfile("make_scenarios/scenarios_subs.csv.gz"))

# Make a set with scenarios to run again (models are saved to disk, need to do this again to check convergence)
scen_rerun <- bind_rows(scen %>% filter(save_model), scen_subs %>% filter(save_model)) %>%
  group_by(scenario_id) %>%
  slice(1:6)

# Make new starting values & seeds for each model
set.seed(22244214)
scen_rerun$dsim_seed <- sample(1:5678623, nrow(scen_rerun), replace=FALSE)
scen_rerun$model_seed <- sample(1:5678623, nrow(scen_rerun), replace=FALSE)

# Create start values
## TPM gamma
diag_value <- runif(nrow(scen_rerun), 0.5, 0.8)
assertthat::are_equal(length(unique(diag_value)), nrow(scen_rerun))
start_gamma <- lapply(diag_value, function(x) {
  gam <- diag(x, 3)
  gam[lower.tri(gam) | upper.tri(gam)] <- (1-x) / 2
  # Assert row sums equal to one
  assertthat::are_equal(rowSums(gam), c(1,1,1))
  # Return json array with flat vector
  return(jsonlite::toJSON(list("tpm" = as.vector(gam)),
                          pretty = TRUE,
                          flatten = TRUE))
})
## Emission distributions
# For the 3 continuous emission distributions
start_emiss <- lapply(1:nrow(scen_rerun), function(x) {
  jsonlite::toJSON(list(
    #EEG_Fpz_Cz_max_gamma
    EEG_Fpz_Cz_mean_beta = c( -0.1 + runif(1, -.2, .2), 0.2 + runif(1, -.1,.1),
                              -0.65 + runif(1, -.2, .2), 0.2 + runif(1, -.1,.1),
                              0.8 + runif(1, -.2, .2), 0.2 + runif(1, -.1,.1)),
    # EOG_median_theta
    EOG_median_theta = c( 0.8 + runif(1, -.2, .2), 0.2 + runif(1, -.1,.1),
                          -1 + runif(1, -.2, .2), 0.2 + runif(1, -.1,.1),
                          0 + runif(1, -.2, .2), 0.2 + runif(1, -.1,.1)),
    # EOG_min_beta
    EOG_min_beta = c(  0.9 + runif(1, -.2, .2), 0.2 + runif(1, -.1,.1),
                       -1 + runif(1, -.2, .2), 0.2 + runif(1, -.1,.1),
                       0 + runif(1, -.2, .2), 0.2 + runif(1, -.1,.1))
  ), pretty = TRUE, flatten = TRUE)
})
# Add to data
scen_rerun$start_gamma <- as.character(start_gamma)
scen_rerun$start_emiss <- as.character(start_emiss)

# Save
write.csv(scen_rerun, file=gzfile("make_scenarios/scenarios_rerun.csv.gz"))
