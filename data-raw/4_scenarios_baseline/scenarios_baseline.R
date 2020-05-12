# Make simulation scenarios for my simulation study
# Written by: Jasper Ginn

rm(list=ls())
# Load sleepsimR
library(sleepsimR)

opts <- getOption("sleepsimR_simulate")
# Spread means
opts$emission_bar <- opts$emission_bar_original_failed
opts$emission_bar[[1]][,1] <- c(-3.9, -1, 2.4)
opts$emission_bar[[2]][,1] <- c(3.05, -3.4, -0.5)
opts$emission_bar[[3]][,1] <- c(0.4, 3.5, -2.8)
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
  filter(n %in% c(40, 80),
         zeta %in% c(0.25, 0.5),
         Q == 0.1,
         n_t %in% c(800, 1600)) %>%
  filter(!(n == 40 & n_t == 1600))

# Add modesl with n=140
scen_largeN <- scen_baseline %>%
  filter(n == 40) %>%
  mutate(n = 140)

# Add to baseline
scen_baseline <- scen_baseline %>%
  bind_rows(scen_largeN)

set.seed(46203)
# Copy models that should be run twice (for convergence checks)
scen_baseline_rerun <- scen_baseline %>%
  filter(save_model) %>%
  mutate(model_seed = sample.int(99999999, n())) %>%
  mutate(iteration_id = paste0(iteration_id, "_chain2"))

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

sim$n <- 20
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
tdf <- tdf %>%
  mutate(state= states) %>%
  group_by(state) %>%
  mutate(EEG_mean_beta = EEG_mean_beta / sd(EEG_mean_beta),
         EOG_median_theta = EOG_median_theta / sd(EOG_median_theta),
         EOG_min_beta = EOG_min_beta / sd(EOG_min_beta)) %>%
  ungroup() %>%
  select(-state)

# View
library(ggplot2)
library(tidyr)
tdf %>%
  mutate(state = states) %>%
  gather(var, val, -id, -state) %>%
  #group_by(state, var) %>%
  #mutate(val = val / sd(val)) #%>%
  #group_by(var) %>%
  #mutate(val = (val - mean(val)) / sd(val)) %>%
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
start_values
start_values[[2]] <- start_values[[2]] + runif(6, 0.01, 0.05)
start_values[[3]] <- start_values[[3]] + runif(6, 0.01, 0.05)
start_values[[4]] <- start_values[[4]] + runif(6, 0.01, 0.05)
gamma_diag <- runif(1, 0.5, 0.8)
gamma_off_diag <- (1-gamma_diag) / 2
diag(start_values[[1]]) <- gamma_diag
start_values[[1]][upper.tri(start_values[[1]])] <- gamma_off_diag
start_values[[1]][lower.tri(start_values[[1]])] <- gamma_off_diag
sim$model_seed <- sample.int(99999999, 1)

# Run model
mod2 <- sleepsimR::run_mHMM(tdf, start_values = start_values,
                            mprop = list("m" = 3, "n_dep" = 3),
                            hyperprior_means = hyp_priors,
                            model_seed = sim$model_seed,
                            mcmc_iterations=3250, mcmc_burn_in = 1250,
                            order_data = FALSE)

sleepsimReval::tpp(mod2, mod1, "gamma_int_bar", var=3)
sleepsimReval::compute_grs(mod1, mod2, "gamma_int_bar", var=3)

sleepsimReval::dens_plot(mod1, mod2, "emiss_varmu_bar", var=2)

# Write output

# Baseline scenarios:
#  scen 1: emission distributions spread but high self transition
#  scen 2: same as 1 but lower self-transition.
#           NB. the values to generate data are set in sleepsimR-run program.
scen_baseline_n40 <- scen_baseline %>%
  filter(n == 40)
#  scen 3: very high N
scen_baseline_n140 <- scen_baseline %>%
  filter(n == 140)
#  scen 4: high N
scen_baseline_n80nt800 <- scen_baseline %>%
  filter(n == 80, n_t == 800)
#  scen 5: high N, very high n_t
scen_baseline_n80nt3200 <- scen_baseline %>%
  filter(n == 80, n_t == 1600) %>%
  mutate(n_t = 3200) %>%
  group_by(zeta) %>%
  slice(1:140)

# Save
write.csv(scen_baseline_n80nt3200, file=gzfile("data-raw/4_scenarios_baseline/scenarios_largeNlargeNt.csv.gz"))

# Create a second set
scen_baseline_2 <- scen_baseline
scen_baseline_2$n <- 140

write.csv(scen_baseline, file=gzfile("data-raw/4_scenarios_baseline/scenarios_baseline_n80.csv.gz"))
write.csv(scen_baseline, file=gzfile("data-raw/4_scenarios_baseline/scenarios_baseline.csv.gz"))
write.csv(scen_baseline_2, file=gzfile("data-raw/4_scenarios_baseline/scenarios_baseline_largeN.csv.gz"))

# Legacy ----

# Arrange subset
scen_baseline <- scen_baseline %>%
  arrange(n, n_t, zeta, scenario_id)

# Save scenarios
write.csv(scen, file=gzfile("data-raw/2_make_simulation_scenarios/scenarios_large_N.csv.gz"))
write.csv(scen_subs, file=gzfile("data-raw/2_make_simulation_scenarios/scenarios_subs.csv.gz"))

# Make a set with scenarios to run again (models are saved to disk, need to do this again to check convergence)
scen_rerun <- bind_rows(scen %>% filter(save_model), scen_subs %>% filter(save_model)) %>%
  group_by(scenario_id) %>%
  slice(1:3)

# Make new starting values & seeds for each model
set.seed(22244214)
scen_rerun$model_seed <- sample(1:99999999, nrow(scen_rerun), replace=FALSE)

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
  jsonlite::toJSON(
    lapply(jsonlite::fromJSON(scen_rerun$start_emiss[x]), function(y) y + runif(6, -0.03, 0.03)),
    pretty = TRUE,
    flatten = TRUE
  )
})
# Add to data
scen_rerun$start_gamma <- as.character(start_gamma)
scen_rerun$start_emiss <- as.character(start_emiss)

# Test two models ----

# Get scenario
sim_rer <- scen_rerun %>%
  filter(n==20,
         n_t==800,
         zeta == 0.25,
         Q == 0.1) %>%
  slice(1) %>%
  as.list()

# Model 2
sim_orig <- scen %>%
  filter(iteration_id == sim_rer$iteration_id) %>%
  as.list()

# Make dataset
data_simulated <- simulate_dataset(sim_rer$n, sim_rer$n_t, sim_rer$zeta, sim_rer$Q, sim_rer$dsim_seed)
# To data frame
tdf <- data.frame(
  id = data_simulated$obs[,1],
  EEG_mean_beta = data_simulated$obs[,2],
  EOG_median_theta = data_simulated$obs[,3],
  EOG_min_beta = data_simulated$obs[,4]
)
states <- data_simulated$states[,2]

# Plot
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
# Reshape start values (model 1)
sim_rer$start_gamma <- jsonlite::parse_json(sim_rer$start_gamma)
sim_rer$start_emiss <- jsonlite::parse_json(sim_rer$start_emiss)
m <- sqrt(length(sim_rer$start_gamma$tpm))
start_values1 <- list(
  matrix(unlist(sim_rer$start_gamma$tpm), nrow=m, ncol=m,
         byrow = TRUE),
  matrix(unlist(sim_rer$start_emiss$EEG_Fpz_Cz_mean_beta),
         ncol=2, byrow = TRUE),
  matrix(unlist(sim_rer$start_emiss$EOG_median_theta),
         ncol=2, byrow = TRUE),
  matrix(unlist(sim_rer$start_emiss$EOG_min_beta),
         ncol=2, byrow=TRUE)
)
# Run model
mod1 <- sleepsimR::run_mHMM(tdf, start_values = start_values1,
                            mprop = list("m" = 3, "n_dep" = 3),
                            hyperprior_means = hyp_priors,
                            model_seed = sim_rer$model_seed,mcmc_iterations=3250, mcmc_burn_in = 1250,
                            order_data = FALSE)

# Reshape start values
sim_orig$start_gamma <- jsonlite::parse_json(sim_orig$start_gamma)
sim_orig$start_emiss <- jsonlite::parse_json(sim_orig$start_emiss)
m <- sqrt(length(sim_orig$start_gamma$tpm))
start_values2 <- list(
  matrix(unlist(sim_orig$start_gamma$tpm), nrow=m, ncol=m,
         byrow = TRUE),
  matrix(unlist(sim_orig$start_emiss$EEG_Fpz_Cz_mean_beta),
         ncol=2, byrow = TRUE),
  matrix(unlist(sim_orig$start_emiss$EOG_median_theta),
         ncol=2, byrow = TRUE),
  matrix(unlist(sim_orig$start_emiss$EOG_min_beta),
         ncol=2, byrow=TRUE)
)

# Run model
mod2 <- sleepsimR::run_mHMM(tdf, start_values = start_values,
                            mprop = list("m" = 3, "n_dep" = 3),
                            hyperprior_means = hyp_priors,
                            model_seed = sim_orig$model_seed,mcmc_iterations=3250, mcmc_burn_in = 1250,
                            order_data = FALSE)

sleepsimReval::tpp(mod1, mod2, "emiss_mu_bar", var=1)
sleepsimReval::compute_grs(mod1, mod2, "emiss_mu_bar", var=1)
sleepsimReval::tpp(mod1, mod2, "emiss_mu_bar", var=2)
sleepsimReval::compute_grs(mod1, mod2, "emiss_mu_bar", var=2)
sleepsimReval::tpp(mod1, mod2, "emiss_mu_bar", var=3)
sleepsimReval::compute_grs(mod1, mod2, "emiss_mu_bar", var=3)

sleepsimReval::tpp(mod1, mod2, "emiss_varmu_bar", var=1)
sleepsimReval::compute_grs(mod1, mod2, "emiss_varmu_bar", var=1)
sleepsimReval::tpp(mod1, mod2, "emiss_varmu_bar", var=2)
sleepsimReval::compute_grs(mod1, mod2, "emiss_varmu_bar", var=2)
sleepsimReval::tpp(mod1, mod2, "emiss_varmu_bar", var=3)
sleepsimReval::compute_grs(mod1, mod2, "emiss_varmu_bar", var=3)

sleepsimReval::tpp(mod1, mod2, "gamma_int_bar", var=1)
sleepsimReval::compute_grs(mod1, mod2, "gamma_int_bar", var=1)

# Save
write.csv(scen_rerun, file=gzfile("data-raw/2_make_simulation_scenarios/scenarios_rerun.csv.gz"))

# Save as datasets in this library
usethis::use_data(scen)
usethis::use_data(scen_subs)
usethis::use_data(scen_rerun, overwrite = TRUE)
