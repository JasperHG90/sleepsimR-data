## Compute bias, MSE, emperical SE etc.

library(sleepsimReval)
library(tidyr)
library(purrr)
library(stringr)
library(dplyr)

# Emission means -----

# Read data
emiss_means_est <- readRDS("data-raw/3_preprocess_simulation_results/data_preprocessed/emission_means_est.rds")
emiss_means_cci <- readRDS("data-raw/3_preprocess_simulation_results/data_preprocessed/emiss_mu_cci.rds")
emiss_means_true_val <- readRDS("data-raw/3_preprocess_simulation_results/data_preprocessed/emiss_means_true_mapval.rds")

# Reshape CCI data
emiss_means_cci_wide <- emiss_means_cci %>%
  select(-n, -n_t, -zeta, -Q, -scenario_id) %>%
  gather(var, val, -iteration_id) %>%
  separate(var, into = c("var1", "var2", "var3", "statistic", "state", "direction"),
           remove = FALSE) %>%
  # NB: statistic doesn't mean anything here. This is a naming error.
  mutate(emiss_var_short = paste0(var1, "_", var2, "_", var3),
         emiss_var = paste0(emiss_var_short, "_", state)) %>%
  select(-var1,-var2,-var3, -statistic) %>%
  # Spread data --> lower and upper get own columns
  pivot_wider(names_from = direction,
              values_from = val,
              id_cols = c("emiss_var", "emiss_var_short",
                          "state", "iteration_id")) %>%
  # Join with original data for scenarios
  left_join(emiss_means_cci %>%
              select(iteration_id, scenario_id)) %>%
  # Arrange
  arrange(iteration_id, emiss_var, emiss_var_short, state)

# Prep the emission mean data.
prepped_emiss_data <- emiss_means_est %>%
  gather(var, val, -scenario_id, - iteration_id,
         -n, -n_t, -zeta, -Q) %>%
  # Add state // emission distribution name // statistic
  separate(var, into = c("var1", "var2", "var3", "statistic", "state"),
           remove=FALSE) %>%
  mutate(emiss_var_short = paste0(var1, "_", var2, "_", var3),
         emiss_var = var) %>%
  select(-var1, -var2, -var3, -var) %>%
  # Keep median MAP estimates
  filter(statistic %in% c("median", "SE")) %>%
  # Rename SE variables
  #mutate(emiss_var = str_replace(emiss_var, "SE", "median"))
  # Spread median and SE into their own columns
  pivot_wider(names_from = statistic,
              values_from = val,
              id_cols = c("iteration_id", "emiss_var_short", "state")) %>%
  # Join with true values
  left_join(emiss_means_true_val, by = c("emiss_var_short", "state")) %>%
  # Get CCI data
  left_join(emiss_means_cci_wide %>%
              select(-emiss_var),
            by=c("iteration_id", "emiss_var_short", "state")) %>%
  # Get scenario params
  left_join(emiss_means_est %>%
              select(iteration_id, scenario_id, n, n_t, zeta, Q),
            by = c("iteration_id", "scenario_id"))

# Compute evaluation statistics
emiss_means_eval_params <- prepped_emiss_data %>%
  # Group by scenario id & emission distribution
  group_by(scenario_id, emiss_var_short, state) %>%
  # Summarize by bias, MSE, empirical SE etc.
  do(summarize_simulation_scenario(.$true_val, .$median,
                                   .$lower, .$upper, .$SE,
                                   compute_multimodal = TRUE)) %>%
  ungroup() %>%
  # Retrieve simulation scenario parameters
  left_join(prepped_emiss_data %>%
              select(scenario_id, emiss_var_short, state, n, n_t, zeta, Q, true_val) %>%
              distinct(),
            by = c("scenario_id", "emiss_var_short", "state")) %>%
  # Arrange
  arrange(desc(n), desc(n_t), desc(zeta), desc(Q), scenario_id, emiss_var_short, state)

# To list
simulation_data_emiss_means <- list(
  "data_preprocessed" = prepped_emiss_data,
  "summary_by_scenario" = emiss_means_eval_params
)
# Use data in library
usethis::use_data(simulation_data_emiss_means)

# Emission distribution between-subject variances ----

rm(list=ls())
# Read data
emiss_varmu_est <- readRDS("data-raw/3_preprocess_simulation_results/data_preprocessed/emission_varmu_est.rds")
emiss_varmu_cci <- readRDS("data-raw/3_preprocess_simulation_results/data_preprocessed/emiss_varmu_cci.rds")
# True values are in the dataset already (zeta)

# Reshape CCI data
emiss_varmu_cci_wide <- emiss_varmu_cci %>%
  select(-n, -n_t, -zeta, -Q, -scenario_id) %>%
  gather(var, val, -iteration_id) %>%
  separate(var, into = c("var1", "var2", "var3", "statistic", "state", "direction"),
           remove = FALSE) %>%
  # NB: statistic doesn't mean anything here. This is a naming error.
  mutate(emiss_var_short = paste0(var1, "_", var2, "_", var3),
         emiss_var = paste0(emiss_var_short, "_", state)) %>%
  select(-var1,-var2,-var3, -statistic) %>%
  # Spread data --> lower and upper get own columns
  pivot_wider(names_from = direction,
              values_from = val,
              id_cols = c("emiss_var", "emiss_var_short",
                          "state", "iteration_id")) %>%
  # Join with original data for scenarios
  left_join(emiss_varmu_cci %>%
              select(iteration_id, scenario_id)) %>%
  # Arrange
  arrange(iteration_id, emiss_var, emiss_var_short, state)

# Prep the emission mean data.
prepped_emiss_varmu_data <- emiss_varmu_est %>%
  gather(var, val, -scenario_id, - iteration_id,
         -n, -n_t, -zeta, -Q) %>%
  # Add state // emission distribution name // statistic
  separate(var, into = c("var1", "var2", "var3", "statistic", "state"),
           remove=FALSE) %>%
  mutate(emiss_var_short = paste0(var1, "_", var2, "_", var3),
         emiss_var = var) %>%
  select(-var1, -var2, -var3, -var) %>%
  # Keep median MAP estimates
  filter(statistic %in% c("median", "SE")) %>%
  # Rename SE variables
  #mutate(emiss_var = str_replace(emiss_var, "SE", "median"))
  # Spread median and SE into their own columns
  pivot_wider(names_from = statistic,
              values_from = val,
              id_cols = c("iteration_id", "emiss_var_short", "state")) %>%
  # Get CCI data
  left_join(emiss_varmu_cci_wide %>%
              select(-emiss_var),
            by=c("iteration_id", "emiss_var_short", "state")) %>%
  # Get scenario params
  left_join(emiss_varmu_est %>%
              select(iteration_id, scenario_id, n, n_t, zeta, Q),
            by = c("iteration_id", "scenario_id"))

# Compute evaluation statistics
emiss_varmu_eval_params <- prepped_emiss_varmu_data %>%
  # Group by scenario id & emission distribution
  group_by(scenario_id, emiss_var_short, state) %>%
  # Summarize by bias, MSE, empirical SE etc.
  do(summarize_simulation_scenario(.$zeta, .$median,
                                   .$lower, .$upper, .$SE,
                                   compute_multimodal = TRUE)) %>%
  ungroup() %>%
  # Retrieve simulation scenario parameters
  left_join(prepped_emiss_varmu_data %>%
              select(scenario_id, emiss_var_short, state, n, n_t, zeta, Q) %>%
              distinct(),
            by = c("scenario_id", "emiss_var_short", "state")) %>%
  # Arrange
  arrange(desc(n), desc(n_t), desc(zeta), desc(Q), scenario_id, emiss_var_short, state)

# Save data to disk
simulation_data_emiss_varmu <- list(
  "data_preprocessed" = prepped_emiss_varmu_data,
  "summary_by_scenario" = emiss_varmu_eval_params
)
usethis::use_data(simulation_data_emiss_varmu)

# Transition probabilities ----

rm(list=ls())
# Read data
gamma_prob_est <- readRDS("data-raw/3_preprocess_simulation_results/data_preprocessed/gamma_int_prob.rds")
gamma_prob_cci <- readRDS("data-raw/3_preprocess_simulation_results/data_preprocessed/gamma_prob_cci.rds")
gamma_true_val <- readRDS("data-raw/3_preprocess_simulation_results/data_preprocessed/tpm_probs_true_mapval.rds")

# Reshape CCI data
gamma_prob_cci_wide <- gamma_prob_cci %>%
  select(-n, -n_t, -zeta, -Q, -scenario_id) %>%
  gather(var, val, -iteration_id) %>%
  separate(var, into=c("int", "transition", "direction")) %>%
  select(-int) %>%
  # Spread data --> lower and upper get own columns
  pivot_wider(names_from = direction,
              values_from = val,
              id_cols = c("iteration_id", "transition")) %>%
  # Join with original data for scenarios
  left_join(gamma_prob_cci %>%
              select(iteration_id, scenario_id))

# Prep the gamma mean data.
prepped_gamma_prob_data <- gamma_prob_est %>%
  gather(var, val, -scenario_id, -iteration_id,
         -n, -n_t, -zeta, -Q) %>%
  separate(var, into=c("transition", "statistic")) %>%
  rename(probability = val) %>%
  # Turn into wide dataset to put the median/SE in their own columns
  pivot_wider(
    names_from = "statistic",
    values_from = "probability"
  ) %>%
  # Get CCI data
  left_join(gamma_prob_cci_wide,
            by=c("iteration_id", "scenario_id", "transition")) %>%
  # Join with true values
  left_join(gamma_true_val, by="transition")

# Compute evaluation statistics
gamma_prob_eval_params <- prepped_gamma_prob_data %>%
  # Group by scenario id & gammaion distribution
  group_by(scenario_id, transition) %>%
  # Summarize by bias, MSE, empirical SE etc.
  do(summarize_simulation_scenario(.$true_val, .$median,
                                   .$lower, .$upper, .$SE,
                                   compute_multimodal = TRUE)) %>%
  ungroup() %>%
  # Retrieve simulation scenario parameters
  left_join(prepped_gamma_prob_data %>%
              select(scenario_id, transition, n, n_t, zeta, Q) %>%
              distinct(),
            by = c("scenario_id", "transition")) %>%
  # Arrange
  arrange(desc(n), desc(n_t), desc(zeta), desc(Q), scenario_id, transition)

# Save data to disk
simulation_data_gamma_prob <- list(
  "data_preprocessed" = prepped_gamma_prob_data,
  "summary_by_scenario" = gamma_prob_eval_params
)
usethis::use_data(simulation_data_gamma_prob)
