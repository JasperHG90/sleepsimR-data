## Preprocess baseline results

library(sleepsimReval)
library(dplyr)
library(purrr)
library(stringr)

# Results for scenario 2 (means are spread out + lower self-probabilities)
# Results for scenario 1 (means are spread out)
# Run this file twice with different SCEN settings (1 or 2)
SCEN <- 5

# Folder paths
fp <- paste0("data-raw/5_preprocess_baseline_results/results_scen",SCEN)

# Full file paths
f <- file.path(fp,
                list.files(fp))

# Read scenario data
scen <- read.csv("data-raw/4_scenarios_baseline/scenarios_baseline.csv.gz", stringsAsFactors = FALSE) %>%
  select(iteration_id, scenario_id,
         n, n_t, zeta, Q)
scen2 <- read.csv("data-raw/4_scenarios_baseline/scenarios_baseline_largeN.csv.gz", stringsAsFactors = FALSE) %>%
  select(iteration_id, scenario_id,
         n, n_t, zeta, Q)
scen3 <- read.csv("data-raw/4_scenarios_baseline/scenarios_baseline_n80.csv.gz", stringsAsFactors = FALSE) %>%
  select(iteration_id, scenario_id,
         n, n_t, zeta, Q)
scen4 <- read.csv("data-raw/4_scenarios_baseline/scenarios_largeNlargeNt.csv.gz", stringsAsFactors = FALSE) %>%
  select(iteration_id, scenario_id,
         n, n_t, zeta, Q)

# Bind rows
scen <- scen %>%
  bind_rows(scen2) %>%
  bind_rows(scen3) %>%
  bind_rows(scen4)

# Read results from disk ----

r <- lapply(f, function(x) parse_sleepsimR_result(x))

# Extract iteration ids
iteration_uid <- map_chr(r, function(x) x$iteration_uid)

# Emission distribution means

# Only emission distributions
em_mu_only <- map(r, function(x) x$emiss_mu_bar)

# Get parameter estimates of emission distributions means
em_mean_est <- lapply(em_mu_only, function(x) {
  sleepsimReval:::postprocess_param_est(x, 3)
}) %>%
  bind_rows() %>%
  mutate(iteration_id = iteration_uid) %>%
  left_join(scen, by="iteration_id")

# Write to disk
saveRDS(em_mean_est, paste0("data-raw/5_preprocess_baseline_results/data_preprocessed/emission_means_est_scen",
                            SCEN,
                            ".rds"))

# Clean
rm(em_mu_only, em_mean_est)

# Emission distribution variances

# Only variance estimates
em_var_only <- map(r, function(x) x$emiss_varmu_bar)

# Retrieve parameter estimates
em_var_est <- lapply(em_var_only,
                            function(x) {
                              sleepsimReval:::postprocess_param_est(x, 3)
                            }) %>%
  bind_rows() %>%
  mutate(iteration_id = iteration_uid) %>%
  left_join(scen, by="iteration_id")

# Write to disk
saveRDS(em_var_est, paste0("data-raw/5_preprocess_baseline_results/data_preprocessed/emission_varmu_est_scen",
                           SCEN,
                           ".rds"))

# Clean
rm(em_var_only, em_var_est)

# Transition probabilities

# Only probabilities
gam_int_prob <- map(r, function(x) x$gamma_prob_bar)

# Retrieve parameter estimates
gam_prob_est <- lapply(gam_int_prob,
                              function(x) {
                                sleepsimReval:::postprocess_gamma_int(x, 3)
                              }) %>%
  bind_rows() %>%
  mutate(iteration_id = iteration_uid) %>%
  left_join(scen, by="iteration_id")

# Write to disk
saveRDS(gam_prob_est, paste0("data-raw/5_preprocess_baseline_results/data_preprocessed/gamma_int_prob_scen",
        SCEN,
        ".rds"))

# Clean
rm(gam_prob_est, gam_int_prob)

# Extract credible intervals ----

# Emission means

# Subset
emiss_mu_cci <- map(r, function(x) x$credible_intervals$emiss_mu_bar)

# Retrieve parameter estimates
emiss_mu_cci_est <- lapply(emiss_mu_cci,
                                  function(x) {
                                    ci <- x
                                    nams <- names(ci)
                                    for(idx in seq_along(nams)) {
                                      names(ci[[idx]]) <- paste0("mean_",
                                                                 c("state1_lower", "state1_upper",
                                                                   "state2_lower", "state2_upper",
                                                                   "state3_lower", "state3_upper"))
                                    }
                                    y <- do.call(cbind.data.frame, ci)
                                    colnames(y) <- gsub("\\.", "_", colnames(y))
                                    y
                                  }) %>%
  bind_rows() %>%
  mutate(iteration_id = iteration_uid) %>%
  left_join(scen, by="iteration_id")

# Write to disk
saveRDS(emiss_mu_cci_est, paste0("data-raw/5_preprocess_baseline_results/data_preprocessed/emiss_mu_cci_scen",
        SCEN, ".rds"))

# Clean
rm(emiss_mu_cci, emiss_mu_cci_est)

# Emission variances -----

# Subset
emiss_varmu_cci <- map(r, function(x) x$credible_intervals$emiss_varmu_bar)

# Retrieve parameter estimates
emiss_varmu_cci_est <- lapply(emiss_varmu_cci,
                                     function(x) {
                                       ci <- x
                                       nams <- names(ci)
                                       for(idx in seq_along(nams)) {
                                         names(ci[[idx]]) <- paste0("varmu_",
                                                                    c("state1_lower", "state1_upper",
                                                                      "state2_lower", "state2_upper",
                                                                      "state3_lower", "state3_upper"))
                                       }
                                       y <- do.call(cbind.data.frame, ci)
                                       colnames(y) <- gsub("\\.", "_", colnames(y))
                                       y
                                     }) %>%
  bind_rows() %>%
  mutate(iteration_id = iteration_uid) %>%
  left_join(scen, by="iteration_id")

# Write to disk
saveRDS(emiss_varmu_cci_est, paste0("data-raw/5_preprocess_baseline_results/data_preprocessed/emiss_varmu_cci_scen",
                                 SCEN, ".rds"))

# Clean
rm(emiss_varmu_cci, emiss_varmu_cci_est)

# Transition probabilities -----

# Subset
gamma_prob_cci <- map(r, function(x) x$credible_intervals$gamma_prob_bar)

# Retrieve parameter estimates
gamma_prob_cci_est <- lapply(gamma_prob_cci,
                                    function(x) {
                                      ci <- x
                                      io <- c("int_S1toS1_lower", "int_S1toS1_upper",
                                              "int_S1toS2_lower", "int_S1toS2_upper",
                                              "int_S1toS3_lower", "int_S1toS3_upper",
                                              "int_S2toS1_lower", "int_S2toS1_upper",
                                              "int_S2toS2_lower", "int_S2toS2_upper",
                                              "int_S2toS3_lower", "int_S2toS3_upper",
                                              "int_S3toS1_lower", "int_S3toS1_upper",
                                              "int_S3toS2_lower", "int_S3toS2_upper",
                                              "int_S3toS3_lower", "int_S3toS3_upper")
                                      y <- data.frame(ci)
                                      colnames(y) <- io
                                      y
                                    }) %>%
  bind_rows() %>%
  mutate(iteration_id = iteration_uid) %>%
  left_join(scen, by="iteration_id")

# Write to disk
saveRDS(gamma_prob_cci_est, paste0("data-raw/5_preprocess_baseline_results/data_preprocessed/gamma_prob_cci_est_scen",
                                 SCEN, ".rds"))

# Clean
rm(gamma_prob_cci, gamma_prob_cci_est)

# Make dataset with true parameter MAP values ----

# Real values (means)
rv <- data.frame(
  "var" = c("EEG_mean_beta_median_state1", "EEG_mean_beta_median_state2",
            "EEG_mean_beta_median_state3", "EOG_median_theta_median_state1",
            "EOG_median_theta_median_state2", "EOG_median_theta_median_state3",
            "EOG_min_beta_median_state1", "EOG_min_beta_median_state2",
            "EOG_min_beta_median_state3"),
  "emiss_var_short" = c("EEG_mean_beta", "EEG_mean_beta",
                        "EEG_mean_beta", "EOG_median_theta",
                        "EOG_median_theta", "EOG_median_theta",
                        "EOG_min_beta", "EOG_min_beta",
                        "EOG_min_beta"),
  "state" = c("state1", "state2",
              "state3", "state1",
              "state2", "state3",
              "state1", "state2",
              "state3"),
  "true_val" = c(-3.9, -1, 2.4, 3.05, -3.4, -0.5, 0.4, 3.5, -2.8),
  stringsAsFactors = FALSE
)
saveRDS(rv, "data-raw/5_preprocess_baseline_results/data_preprocessed/emiss_means_true_mapval.rds")

# True values (tpm)
library(sleepsimR)
opts <- getOption("sleepsimR_simulate")
tpm_scen1 <- data.frame(
  "transition" = c(
    'S1toS1','S1toS2','S1toS3',
    'S2toS1','S2toS2','S2toS3',
    'S3toS1','S3toS2','S3toS3'
  ),
  "true_val" = as.vector(opts$gamma_bar),
  stringsAsFactors = FALSE
)
saveRDS(tpm_scen1, "data-raw/5_preprocess_baseline_results/data_preprocessed/tpm_probs_true_mapval_scen1.rds")
tpm_scen2 <- data.frame(
  "transition" = c(
    'S1toS1','S1toS2','S1toS3',
    'S2toS1','S2toS2','S2toS3',
    'S3toS1','S3toS2','S3toS3'
  ),
  "true_val" = c(0.8, 0.1, 0.1,0.15, 0.7, 0.15, 0.18, 0.18, 0.64),
  stringsAsFactors = FALSE
)
saveRDS(tpm_scen2, "data-raw/5_preprocess_baseline_results/data_preprocessed/tpm_probs_true_mapval_scen2.rds")
