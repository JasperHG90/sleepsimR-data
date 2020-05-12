## Test

rm(list=ls())
run_single_varmHMM <- function(data, start_values, mprop, hyperprior_means, model_seed, mcmc_iterations = 2000, mcmc_burn_in = 1000, show_progress = TRUE) {
  # Mcmc options
  mcmcOpts <- list(
    "J"=mcmc_iterations,
    "burn_in"=mcmc_burn_in
  )
  # Order!
  state_orders <- vector("list", mprop$n_dep)
  names(state_orders) <- colnames(data[,-1])
  # Set hyper-prior values
  hyper_priors <- list(
    # Hyperprior on intercepts of dependent variables
    emiss_mu0 = list(matrix(hyperprior_means[[1]], nrow=1)),
    # Hyperprior on the number of subjects in each state
    # Hypothetical subjects --> c(1,1,1)
    emiss_K0 = list(1,1,1),
    # Degrees of freedom on the emission dist. means
    # Hypothetical degrees of freedom
    emiss_nu = list(1,1,1),
    # Variances between subjects
    # Prior on hyperparameter between-subject variance
    # Hypothetical variances between hypothetical subjects
    emiss_V = list(
      # Depvar1
      c(.1,.1,.1)
    ),
    # shape values. Fixed variances of normal emission distributions
    # SUbject-fixed normal emission distribution shape/scale parameters
    # This is a regular inverse gamma
    emiss_a0 = list(
      # Depvar1
      c(.01,.01,.01)
    ),
    # Hyperprior on scale values of inverse gamma
    emiss_b0 = list(
      # Depvar1
      c(.01,.01,.01)
    )
  )
  # Set seed
  set.seed(model_seed)
  # Run model
  mod <- mHMMbayes::mHMM_cont(as.matrix(data),
                              gen=mprop,
                              start_val=start_values,
                              mcmc = mcmcOpts,
                              emiss_hyp_prior = hyper_priors,
                              show_progress=show_progress)
  # Add state orders
  mod$state_orders <- state_orders
  # Return model
  return(mod)
}

# Load sleepsimR
library(sleepsimR)

opts <- getOption("sleepsimR_simulate")
# Spread means
opts$emission_bar <- opts$emission_bar_original_failed
opts$emission_bar[[1]][,1] <- c(-3.9, -1, 2.4)
opts$emission_bar[[1]][,2] <- 0.5
opts$emission_bar[[2]][,1] <- c(3.05, -3.4, -0.5)
opts$emission_bar[[2]][,2] <- 0.5
opts$emission_bar[[3]][,1] <- c(0.4, 3.5, -2.8)
opts$emission_bar[[3]][,2] <- 0.5
options("sleepsimR_simulate" = opts)

# Make scenarios
scen <- generate_scenarios(seed = 912326)

# Subset scenarios
unique(scen$zeta)
# Rank
library(dplyr)
scen <- scen %>%
  arrange(desc(n), desc(n_t), scenario_id)

# Subset for baseline
library(jsonlite)
scen_baseline <- scen %>%
  filter(n %in% c(40),
         zeta %in% c(0.25),
         Q == 0.1,
         n_t %in% c(800))

set.seed(46203)
# Copy models that should be run twice (for convergence checks)
scen_baseline_rerun <- scen_baseline %>%
  filter(save_model) %>%
  mutate(model_seed = sample.int(99999999, n()))

# Bind
scen_baseline <- scen_baseline %>%
  bind_rows(scen_baseline_rerun)

# Make new start values
# For the 3 continuous emission distributions

start_emiss <- unlist(lapply(1:nrow(scen_baseline), function(x) {
  jsonlite::toJSON(list(
    #EEG_Fpz_Cz_max_gamma
    EEG_Fpz_Cz_mean_beta = c( -3.9 + runif(1, -.2, .2), 0.2 + runif(1, -.1,.1),
                              -1 + runif(1, -.2, .2), 0.2 + runif(1, -.1,.1),
                              2.4 + runif(1, -.2, .2), 0.2 + runif(1, -.1,.1)),
    # EOG_median_theta
    EOG_median_theta = c( 3.05 + runif(1, -.2, .2), 0.2 + runif(1, -.1,.1),
                          -3.4 + runif(1, -.2, .2), 0.2 + runif(1, -.1,.1),
                          -0.5 + runif(1, -.2, .2), 0.2 + runif(1, -.1,.1)),
    # EOG_min_beta
    EOG_min_beta = c(  0.4 + runif(1, -.2, .2), 0.2 + runif(1, -.1,.1),
                       3.5 + runif(1, -.2, .2), 0.2 + runif(1, -.1,.1),
                       -2.8 + runif(1, -.2, .2), 0.2 + runif(1, -.1,.1))
  ), pretty = TRUE, flatten = TRUE)
}))

# Set baseline emission starts
scen_baseline$start_emiss <- start_emiss

# Test a model ----

# Get one
sim <- scen_baseline %>%
  filter(n_t == 800) %>%
  slice(1) %>%
  as.list()

#sim$n <- 20
#sim$zeta <- 1
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
hyp_priors <- list(
  as.vector(em1[sort.list(as.numeric(names(em1)))])
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
# Remove
start_values[[3]] <- start_values[[4]] <- NULL
#start_values
#start_values[[2]] <- start_values[[2]] + runif(6, 0.01, 0.05)
#start_values[[3]] <- start_values[[3]] + runif(6, 0.01, 0.05)
#start_values[[4]] <- start_values[[4]] + runif(6, 0.01, 0.05)
#gamma_diag <- runif(1, 0.5, 0.8)
#gamma_off_diag <- (1-gamma_diag) / 2
#diag(start_values[[1]]) <- gamma_diag
#start_values[[1]][upper.tri(start_values[[1]])] <- gamma_off_diag
#start_values[[1]][lower.tri(start_values[[1]])] <- gamma_off_diag
#sim$model_seed <- sample.int(99999999, 1)

tdf[,3] <- tdf[,4] <- NULL

# Run model
mod <- run_single_varmHMM(tdf, start_values = start_values,
                            mprop = list("m" = 3, "n_dep" = 1),
                            hyperprior_means = hyp_priors,
                            model_seed = sim$model_seed,
                            mcmc_iterations=3250, mcmc_burn_in = 1250)

mapval <- sleepsimR::MAP(mod)
