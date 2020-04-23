## Scenarios for baseline results

rm(list=ls())
library(sleepsimR)
library(dplyr)
library(tidyr)
library(ggplot2)
library(jsonlite)
set.seed(12658854)

# Get options for parameters
opts <- options()$sleepsimR_simulate
#opts_orig <- opts

# Set new mean values for emission distributions
#opts$emission_bar[[1]][,1] <- c(-.5, -5.7, 6)
#opts$emission_bar[[2]][,1] <- c(5, -5, .1)
#opts$emission_bar[[3]][,1] <- c(5, -0.3, -5)
#opts$emission_bar[[1]][,2] <- 0.5
#opts$emission_bar[[2]][,2] <- 0.5
#opts$emission_bar[[3]][,2] <- 0.5

# Set option
#options(sleepsimR_simulate = opts)

# Make scenarios
scen <- generate_scenarios(seed = 6672945)
scen_select <- scen %>%
  filter(n %in% c(20, 40, 80),
         zeta %in% c(0.5, 1, 2))

# Select fifty scenarios for each
scen_select_sub <- scen_select %>%
  group_by(scenario_id) %>%
  arrange(rank) %>%
  slice(1:50) %>%
  # Save models of 3 iterations
  mutate(scol = sample.int(n(), n()),
         save_model = ifelse(scol <= 2, TRUE, FALSE)) %>%
  select(-scol) %>%
  ungroup()

# Select the models that need to run twice (for model convergence check)
models_run_twice <- scen_select_sub %>%
  filter(save_model) %>%
  mutate(iteration_id = paste0(iteration_id, "_chain2"))
# Update model seeds
models_run_twice$model_seed <- sample.int(999999999, nrow(models_run_twice))
# Row-bind
scen_select_sub <- bind_rows(scen_select_sub,
                             models_run_twice)

# Update the starting values
diag_value <- runif(nrow(scen_select_sub), 0.5, 0.8)
start_emiss_new <- rep("", nrow(scen_select_sub))
start_gamma_new <- rep("", nrow(scen_select_sub))
for(idx in seq_along(start_emiss_new)) {
  start_emiss_new[idx] <- jsonlite::toJSON(list(
    #EEG_Fpz_Cz_max_gamma
    EEG_Fpz_Cz_mean_beta = c( -0.5 + runif(1, -.2, .2), 0.5 + scen_select[idx, "zeta"] + runif(1, -.1,.1),
                              -5.7 + runif(1, -.2, .2), 0.5 + scen_select[idx, "zeta"] + runif(1, -.1,.1),
                              6 + runif(1, -.2, .2), 0.5 + scen_select[idx, "zeta"] + runif(1, -.1,.1)),
    # EOG_median_theta
    EOG_median_theta = c( 5 + runif(1, -.2, .2), 0.5 + scen_select[idx, "zeta"] + runif(1, -.1,.1),
                          -5 + runif(1, -.2, .2), 0.5 + scen_select[idx, "zeta"] + runif(1, -.1,.1),
                          .1 + runif(1, -.2, .2), 0.5 + scen_select[idx, "zeta"] + runif(1, -.1,.1)),
    # EOG_min_beta
    EOG_min_beta = c(  5 + runif(1, -.2, .2), 0.5 + scen_select[idx, "zeta"] + runif(1, -.1,.1),
                       -0.3 + runif(1, -.2, .2), 0.5 + scen_select[idx, "zeta"] + runif(1, -.1,.1),
                       -5 + runif(1, -.2, .2), 0.5 + scen_select[idx, "zeta"] + runif(1, -.1,.1))
  ), pretty = TRUE, flatten = TRUE)
  # Gamma
  gam <- diag(diag_value[idx], 3)
  gam[lower.tri(gam) | upper.tri(gam)] <- (1-diag_value[idx]) / 2
  # Assert row sums equal to one
  assertthat::are_equal(rowSums(gam), c(1,1,1))
  # Return json array with flat vector
  start_gamma_new[idx] <- jsonlite::toJSON(list("tpm" = as.vector(gam)),
                                            pretty = TRUE,
                                            flatten = TRUE)
}
# Add to data
scen_select_sub$start_emiss <- start_emiss_new
scen_select_sub$start_gamma <- start_gamma_new

# Test a model ----

# Test a model
sim <- scen_select_sub %>%
  filter(save_model) %>%
  filter(zeta == 2, n == 20,
         Q == 0.4, n_t == 400) %>%
  filter(iteration_id == "2f335f12498486e59435d0d533d7c80c_chain") %>%
  as.list()
# Make dataset
data_simulated <- simulate_dataset(sim$n, sim$n_t, sim$zeta, sim$Q, sim$dsim_seed)
# To data frame
tdf <- data.frame(
  id = data_simulated$obs[,1],
  EEG_mean_beta = data_simulated$obs[,2],
  EOG_median_theta = data_simulated$obs[,3],
  EOG_min_beta = data_simulated$obs[,4],
  state = data_simulated$states[,2]
)
# Get summary stats
tdf_ss <- tdf %>%
  gather(var, val, -state, -id) %>%
  group_by(var, state) %>%
  summarize(mvar = mean(val),
            tvar = var(val))
tdf$state <- NULL
states <- data_simulated$states[,2]

# View
tdf %>%
  as_tibble() %>%
  mutate(state = states) %>%
  ungroup() %>%
  gather(var, val, -state, -id) %>%
  mutate(state = factor(state, labels=c("state1", "state2", "state3"),
                        levels = c(1,2,3))) %>%
  ggplot(aes(x=val, fill = state)) +
    geom_density(alpha = 0.4) +
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

# Model values
mprop <- list(
  "n_dep" = 3,
  "m" = 3
)
# Run model
mod1 <- sleepsimR::run_mHMM(tdf, start_values = start_values, mprop = mprop, hyperprior_means = hyp_priors,
                           model_seed = sim$model_seed,mcmc_iterations=3250, mcmc_burn_in = 1250,
                           order_data = FALSE)

mm1 <- sleepsimR::MAP(mod1)
mm2 <- sleepsimR::MAP(mod2)
sleepsimReval::tpp(mod1, mod2, "emiss_mu_bar", var=2)
sleepsimReval::dens_plot(mod1, mod2, "emiss_varmu_bar", var=1)
sleepsimReval::compute_grs(mod1, mod2, "gamma_int_bar", var=3)

# Arrange subset
scen_subs <- scen_subs %>%
  arrange(n, n_t, zeta, scenario_id)

# Save scenarios
write.csv(scen, file=gzfile("data-raw/2_make_simulation_scenarios/scenarios.csv.gz"))
write.csv(scen_subs, file=gzfile("data-raw/2_make_simulation_scenarios/scenarios_subs.csv.gz"))

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
write.csv(scen_rerun, file=gzfile("data-raw/2_make_simulation_scenarios/scenarios_rerun.csv.gz"))

# Save as datasets in this library
usethis::use_data(scen)
usethis::use_data(scen_subs)
usethis::use_data(scen_rerun)



# Split into 4 groups ----

scen_select_sub_orig <- scen_select_sub

# First set
scen_1 <- scen_select_sub %>%
  filter(n %in% c(20, 40),
         n_t == 400,
         save_model == FALSE)
# Filter scenarios
out <- !(scen_select_sub$iteration_id %in% scen_1$iteration_id)
scen_select_sub <- scen_select_sub[out,]
# Second set
scen_2 <- scen_select_sub %>%
  filter(save_model)
out <- !(scen_select_sub$iteration_id %in% scen_2$iteration_id)
# Filter
scen_select_sub <- scen_select_sub[out,]
# Third set
scen_3 <- scen_select_sub %>%
  filter(n == 80)
# Split into two
scen_6 <- scen_3 %>%
  slice(1:630)
# Filter
out <- !(scen_select_sub$iteration_id %in% scen_3$iteration_id)
scen_select_sub <- scen_select_sub[out,]
scen_3 <- scen_3 %>%
  filter(!(iteration_id %in% scen_6$iteration_id))
# Fourth set
scen_4 <- scen_select_sub %>%
  filter(n == 40)
out <- !(scen_select_sub$iteration_id %in% scen_4$iteration_id)
scen_select_sub <- scen_select_sub[out,]

# Store
write.csv(scen_1, file=gzfile("data-raw/4_scenarios_baseline/scen1.csv.gz"))
write.csv(scen_2, file=gzfile("data-raw/4_scenarios_baseline/scen2.csv.gz"))
write.csv(scen_3, file=gzfile("data-raw/4_scenarios_baseline/scen3.csv.gz"))
write.csv(scen_4, file=gzfile("data-raw/4_scenarios_baseline/scen4.csv.gz"))
write.csv(scen_select_sub, file=gzfile("data-raw/4_scenarios_baseline/scen5.csv.gz"))
write.csv(scen_6, file=gzfile("data-raw/4_scenarios_baseline/scen6.csv.gz"))
