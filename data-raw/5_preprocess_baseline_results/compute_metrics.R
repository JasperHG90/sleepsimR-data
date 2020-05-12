## Compute simulation metrics

## Compute bias, MSE, emperical SE etc.
rm(list=ls())
library(sleepsimReval)
library(tidyr)
library(purrr)
library(stringr)
library(dplyr)

# Base folder
SCEN <- 5
BASE_FOLDER <- "data-raw/5_preprocess_baseline_results/data_preprocessed/"
ZETA <- 0.25

# Rounder function
rounder <- function(x, digits=2) {
  as.numeric(format(round(x, digits=digits), nsmall=digits))
}

# Emission means -----

# Read data
emiss_means_est <- readRDS(paste0(BASE_FOLDER, "emission_means_est_scen", SCEN, ".rds"))
emiss_means_cci <- readRDS(paste0(BASE_FOLDER, "emiss_mu_cci_scen", SCEN, ".rds"))
emiss_means_true_val <- readRDS(paste0(BASE_FOLDER, "emiss_means_true_mapval.rds"))

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
  distinct() %>%
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
  distinct() %>%
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
  distinct() %>%
  # Group by scenario id & emission distribution
  group_by(scenario_id, emiss_var_short, state) %>%
  # Summarize by bias, MSE, empirical SE etc.
  do(summarize_simulation_scenario(.$true_val, .$median,
                                   .$lower, .$upper, .$SE,
                                   compute_multimodal = FALSE)) %>%
  ungroup() %>%
  # Retrieve simulation scenario parameters
  left_join(prepped_emiss_data %>%
              select(scenario_id, emiss_var_short, state, n, n_t, zeta, Q, true_val) %>%
              distinct(),
            by = c("scenario_id", "emiss_var_short", "state")) %>%
  # Arrange
  arrange(desc(n), desc(n_t), desc(zeta), desc(Q), scenario_id, emiss_var_short, state) %>%
  # Add percent bias
  mutate(pbias = ((bias) / true_val) * 100,
         pbias_mcmc_se = ((bias_mcmc_se) / true_val) * 100) %>%
  filter(zeta == ZETA) %>%
  select(emiss_var_short,
         state,
         pbias,
         pbias_mcmc_se,
         empirical_se,
         empirical_se_mcmc_se,
         modSE,
         modSE_mcmc_se,
         MSE,
         MSE_mcmc_se,
         coverage,
         coverage_mcmc_se,
         bias_corr_coverage,
         bias_corr_coverage_mcmc_se) %>%
  # Round
  mutate_at(vars(-any_of(c("emiss_var_short", "state"))), list(~ rounder(., digits=3))) %>%
  mutate_at(vars(contains("pbias")), list(~ rounder(., digits=2))) %>%
  mutate_at(vars(contains("coverage")), list(~ rounder(., digits=2))) %>%
  # Paste values + MC Errors
  mutate(pbias = paste0(format(pbias, nsmall=2), " (", format(pbias_mcmc_se, nsmall=2), ")"),
         empirical_se = paste0(format(empirical_se, nsmall=3), " (", format(empirical_se_mcmc_se, nsmall=3), ")"),
         modSE = paste0(format(modSE, nsmall=3), " (", format(modSE_mcmc_se, nsmall=3), ")"),
         MSE = paste0(format(MSE, nsmall=3), " (", format(MSE_mcmc_se, nsmall=3), ")"),
         coverage = paste0(format(coverage, nsmall=2), " (", format(coverage_mcmc_se, nsmall=2), ")"),
         bias_corr_coverage = paste0(format(bias_corr_coverage, nsmall=2), " (", format(bias_corr_coverage_mcmc_se, nsmall=2), ")")) %>%
  # Remove SE values
  select(-ends_with("mcmc_se")) %>%
  # Add order
  mutate(rank = 1:n())

emiss_means_eval_params$rank[4:6] <- emiss_means_eval_params$rank[4:6] + 3
emiss_means_eval_params$rank[7:9] <- emiss_means_eval_params$rank[7:9] + 6

#emiss_means_eval_params[,3:7] <- round(emiss_means_eval_params[,3:7], 3)

# To list
simulation_data_emiss_means <- list(
  "data_preprocessed" = prepped_emiss_data,
  "summary_by_scenario" = emiss_means_eval_params
)
# Use data in library
#usethis::use_data(simulation_data_emiss_means)

# Emission distribution between-subject variances ----

#rm(list=ls())
# Read data
emiss_varmu_est <- readRDS(paste0(BASE_FOLDER, "emission_varmu_est_scen", SCEN, ".rds"))
emiss_varmu_cci <- readRDS(paste0(BASE_FOLDER, "emiss_varmu_cci_scen", SCEN, ".rds"))
# True values are in the dataset already (zeta)

# Reshape CCI data
emiss_varmu_cci_wide <- emiss_varmu_cci %>%
  distinct() %>%
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
  distinct() %>%
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
  distinct() %>%
  # Group by scenario id & emission distribution
  group_by(scenario_id, emiss_var_short, state) %>%
  # Summarize by bias, MSE, empirical SE etc.
  do(summarize_simulation_scenario(.$zeta, .$median,
                                   .$lower, .$upper, .$SE,
                                   compute_multimodal = FALSE)) %>%
  ungroup() %>%
  # Retrieve simulation scenario parameters
  left_join(prepped_emiss_varmu_data %>%
              select(scenario_id, emiss_var_short, state, n, n_t, zeta, Q) %>%
              distinct(),
            by = c("scenario_id", "emiss_var_short", "state")) %>%
  # Arrange
  arrange(desc(n), desc(n_t), desc(zeta), desc(Q), scenario_id, emiss_var_short, state) %>%
  mutate(pbias = ((bias) / zeta) * 100,
         pbias_mcmc_se = ((bias_mcmc_se) / zeta) * 100) %>%
  filter(zeta == ZETA) %>%
  select(emiss_var_short,
         state,
         pbias,
         pbias_mcmc_se,
         bias,
         bias_mcmc_se,
         empirical_se,
         empirical_se_mcmc_se,
         modSE,
         modSE_mcmc_se,
         MSE,
         MSE_mcmc_se,
         coverage,
         coverage_mcmc_se,
         bias_corr_coverage,
         bias_corr_coverage_mcmc_se) %>%
  # Round
  mutate_at(vars(-any_of(c("emiss_var_short", "state"))), list(~ rounder(., digits=3))) %>%
  mutate_at(vars(contains("pbias")), list(~ rounder(., digits=2))) %>%
  mutate_at(vars(contains("coverage")), list(~ rounder(., digits=2))) %>%
  # Paste values + MC Errors
  mutate(pbias = paste0(format(pbias, nsmall=2), " (", format(pbias_mcmc_se, nsmall=2), ")"),
         empirical_se = paste0(format(empirical_se, nsmall=3), " (", format(empirical_se_mcmc_se, nsmall=3), ")"),
         modSE = paste0(format(modSE, nsmall=3), " (", format(modSE_mcmc_se, nsmall=3), ")"),
         MSE = paste0(format(MSE, nsmall=3), " (", format(MSE_mcmc_se, nsmall=3), ")"),
         coverage = paste0(format(coverage, nsmall=2), " (", format(coverage_mcmc_se, nsmall=2), ")"),
         bias_corr_coverage = paste0(format(bias_corr_coverage, nsmall=2), " (", format(bias_corr_coverage_mcmc_se, nsmall=2), ")")) %>%
  # Remove SE values
  select(-ends_with("mcmc_se")) %>%
  # Add order
  mutate(rank = 1:n() + 3)

emiss_varmu_eval_params$rank[4:6] <- emiss_varmu_eval_params$rank[4:6] + 3
emiss_varmu_eval_params$rank[7:9] <- emiss_varmu_eval_params$rank[7:9] + 6

# Save data to disk
simulation_data_emiss_varmu <- list(
  "data_preprocessed" = prepped_emiss_varmu_data,
  "summary_by_scenario" = emiss_varmu_eval_params
)
#usethis::use_data(simulation_data_emiss_varmu)

# Transition probabilities ----

#rm(list=ls())
# Read data
gamma_prob_est <- readRDS(paste0(BASE_FOLDER, "gamma_int_prob_scen", SCEN, ".rds"))
gamma_prob_cci <- readRDS(paste0(BASE_FOLDER, "gamma_prob_cci_est_scen", SCEN, ".rds"))
if(SCEN == 1) {
  gamma_true_val <- readRDS(paste0(BASE_FOLDER, "tpm_probs_true_mapval_scen", SCEN, ".rds"))
} else {
  gamma_true_val <- readRDS(paste0(BASE_FOLDER, "tpm_probs_true_mapval_scen", 2, ".rds"))
}

# Reshape CCI data
gamma_prob_cci_wide <- gamma_prob_cci %>%
  distinct() %>%
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
  distinct() %>%
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
  distinct() %>%
  # Group by scenario id & gammaion distribution
  group_by(scenario_id, transition) %>%
  # Summarize by bias, MSE, empirical SE etc.
  do(summarize_simulation_scenario(.$true_val, .$median,
                                   .$lower, .$upper, .$SE,
                                   compute_multimodal = FALSE)) %>%
  ungroup() %>%
  # Retrieve simulation scenario parameters
  left_join(prepped_gamma_prob_data %>%
              select(scenario_id, transition, n, n_t, zeta, Q) %>%
              distinct(),
            by = c("scenario_id", "transition")) %>%
  # Arrange
  left_join(gamma_true_val) %>%
  arrange(desc(n), desc(n_t), desc(zeta), desc(Q), scenario_id, transition) %>%
  select(-scenario_id) %>%
  filter(zeta == ZETA) %>%
  mutate(pbias = ((bias) / zeta) * 100,
                  pbias_mcmc_se = ((bias_mcmc_se) / zeta) * 100) %>%
  select(transition,
         pbias,
         pbias_mcmc_se,
         empirical_se,
         empirical_se_mcmc_se,
         modSE,
         modSE_mcmc_se,
         MSE,
         MSE_mcmc_se,
         coverage,
         coverage_mcmc_se,
         bias_corr_coverage,
         bias_corr_coverage_mcmc_se) %>%
  # Round
  mutate_at(vars(-any_of("transition")), list(~ rounder(., digits=3))) %>%
  mutate_at(vars(contains("pbias")), list(~ rounder(., digits=2))) %>%
  mutate_at(vars(contains("coverage")), list(~ rounder(., digits=2))) %>%
  # Paste values + MC Errors
  mutate(pbias = paste0(format(pbias, nsmall=2), " (", format(pbias_mcmc_se, nsmall=2), ")"),
         empirical_se = paste0(format(empirical_se, nsmall=3), " (", format(empirical_se_mcmc_se, nsmall=3), ")"),
         modSE = paste0(format(modSE, nsmall=3), " (", format(modSE_mcmc_se, nsmall=3), ")"),
         MSE = paste0(format(MSE, nsmall=3), " (", format(MSE_mcmc_se, nsmall=3), ")"),
         coverage = paste0(format(coverage, nsmall=2), " (", format(coverage_mcmc_se, nsmall=2), ")"),
         bias_corr_coverage = paste0(format(bias_corr_coverage, nsmall=2), " (", format(bias_corr_coverage_mcmc_se, nsmall=2), ")")) %>%
  # Remove SE values
  select(-ends_with("mcmc_se")) %>%
  # Add order
  mutate(rank = 19:(19+(n()-1)))

# Add emiss_var_short
gamma_prob_eval_params <- cbind(data.frame("emiss_var_short" = "MLR_intercepts"), gamma_prob_eval_params)
colnames(gamma_prob_eval_params)[2] <- "state"

# Add everything in a dataframe
scen <- bind_rows(emiss_means_eval_params, emiss_varmu_eval_params) %>%
  bind_rows(gamma_prob_eval_params) %>%
  arrange(rank) %>%
  select(-rank, -bias)

# Assign to environment
assign(paste0("scen", SCEN), scen)

# When all scenarios are computed, column bind
scen_out <- bind_cols(scen1,
                      scen2 %>% select(-emiss_var_short, -state)) %>%
  bind_cols(., scen3 %>%
              select(-emiss_var_short, -state)) %>%
  bind_cols(., scen4 %>%
              select(-emiss_var_short, -state)) %>%
  bind_cols(., scen5 %>%
              select(-emiss_var_short, -state))

# Write
write.csv(scen_out, file = paste0("data-raw/5_preprocess_baseline_results/scen_zeta_05.csv"), row.names = FALSE)
