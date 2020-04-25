# Make simulation scenarios for my simulation study
# Written by: Jasper Ginn

# Load sleepsimR
library(sleepsimR)

opts <- getOption("sleepsimR_simulate")
opts$emission_bar <- opts$emission_bar_original_failed
options("sleepsimR_simulate" = opts)

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
sim$n <- 20
sim$zeta <- 1
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

# View
library(ggplot2)
library(tidyr)
tdf %>%
  mutate(state = states) %>%
  gather(var, val, -id, -state) %>%
  ggplot(aes(x=val, fill = as.factor(state))) +
    geom_density(alpha=0.4) +
    facet_wrap(". ~ var")

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
# start_values
# start_values[[2]] <- start_values[[2]] + runif(6, 0.01, 0.05)
# start_values[[3]] <- start_values[[3]] + runif(6, 0.01, 0.05)
# start_values[[4]] <- start_values[[4]] + runif(6, 0.01, 0.05)
# gamma_diag <- runif(1, 0.5, 0.8)
# gamma_off_diag <- (1-gamma_diag) / 2
# diag(start_values[[1]]) <- gamma_diag
# start_values[[1]][upper.tri(start_values[[1]])] <- gamma_off_diag
# start_values[[1]][lower.tri(start_values[[1]])] <- gamma_off_diag
# sim$model_seed <- sample.int(99999999, 1)

# Run model
mod2 <- sleepsimR::run_mHMM(tdf, start_values = start_values,
                            mprop = list("m" = 3, "n_dep" = 3),
                            hyperprior_means = hyp_priors,
                           model_seed = sim$model_seed,mcmc_iterations=3250, mcmc_burn_in = 1250,
                           order_data = FALSE)

sleepsimReval::tpp(mod2, mod1, "emiss_mu_bar", var=2)
sleepsimReval::compute_grs(mod1, mod2, "emiss_mu_bar", var=2)

sleepsimReval::dens_plot(mod1, mod2, "emiss_mu_bar", var=2)

# Arrange subset
scen_subs <- scen_subs %>%
  arrange(n, n_t, zeta, scenario_id)

# Save scenarios
write.csv(scen, file=gzfile("data-raw/2_make_simulation_scenarios/scenarios.csv.gz"))
write.csv(scen_subs, file=gzfile("data-raw/2_make_simulation_scenarios/scenarios_subs.csv.gz"))

# Make a set with scenarios to run again (models are saved to disk, need to do this again to check convergence)
scen_rerun <- bind_rows(scen %>% filter(save_model), scen_subs %>% filter(save_model)) %>%
  group_by(scenario_id) %>%
  slice(1:3)

# Make new starting values & seeds for each model
set.seed(22244214)
# scen_rerun$dsim_seed <- sample(1:5678623, nrow(scen_rerun), replace=FALSE)
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
    EEG_Fpz_Cz_mean_beta = c( -0.36 + runif(1, -.2, .2), 0.2 + runif(1, -.1,.1),
                              -0.60 + runif(1, -.2, .2), 0.2 + runif(1, -.1,.1),
                              0.7 + runif(1, -.2, .2), 0.2 + runif(1, -.1,.1)),
    # EOG_median_theta
    EOG_median_theta = c( 1.01 + runif(1, -.2, .2), 0.2 + runif(1, -.1,.1),
                          -1.31 + runif(1, -.2, .2), 0.2 + runif(1, -.1,.1),
                          -.24 + runif(1, -.2, .2), 0.2 + runif(1, -.1,.1)),
    # EOG_min_beta
    EOG_min_beta = c(  0.75 + runif(1, -.2, .2), 0.2 + runif(1, -.1,.1),
                       -1.3 + runif(1, -.2, .2), 0.2 + runif(1, -.1,.1),
                       0 + runif(1, -.2, .2), 0.2 + runif(1, -.1,.1))
  ), pretty = TRUE, flatten = TRUE)
})
# Add to data
scen_rerun$start_gamma <- as.character(start_gamma)
scen_rerun$start_emiss <- as.character(start_emiss)

# Save
write.csv(scen_rerun, file=gzfile("data-raw/2_make_simulation_scenarios/scenarios_rerun.csv.gz"))

# Save as datasets in this library
usethis::use_data(scen)
usethis::use_data(scen_subs)
usethis::use_data(scen_rerun)
