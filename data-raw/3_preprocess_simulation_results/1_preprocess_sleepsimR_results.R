# Preprocessing of the simulation results. Make sure to read the full documentation
# here <https://github.com/JasperHG90/sleepsimR-documentation>

library(sleepsimReval)
library(dplyr)
library(future)
library(future.apply)
library(purrr)
library(stringr)

# Set plan for future library
# This starts up a local PSOCK cluster
# With garbage collection.
plan("future::multisession", gc=TRUE, workers = 4)

# Read the individual files into RAM (large object; +- 4GB) -----

# Load result files
f <- paste0("results/", list.files("results"))
t1 <- Sys.time()
res <- future_lapply(f, function(x) parse_sleepsimR_result(x),
                     future.globals = FALSE,
                     future.packages = "sleepsimReval")
Sys.time() - t1

# Save parsed results
saveRDS(res, "data-raw/3_preprocess_simulation_results/sleepsimR_results_parsed.rds")

# Extract parameter estimates from the parsed results ----

# Load data
res <- readRDS("data-raw/3_preprocess_simulation_results/sleepsimR_results_parsed.rds")

# Extract iteration ids
iteration_uid <- map_chr(res, function(x) x$iteration_uid)

# Read scenario data
scen <- bind_rows(
  read.csv("data-raw/2_make_simulation_scenarios/scenarios_subs.csv.gz", stringsAsFactors = FALSE),
  read.csv("data-raw/2_make_simulation_scenarios/scenarios.csv.gz", stringsAsFactors = FALSE)
) %>%
  select(iteration_id, scenario_id,
         n, n_t, zeta, Q)

# Emission distribution means ----

# Only emission distributions
em_mu_only <- map(res, function(x) x$emiss_mu_bar)

# Get parameter estimates of emission distributions means
t1 <- Sys.time()
em_mean_est <- future_lapply(em_mu_only, function(x) {
                        sleepsimReval:::postprocess_param_est(x, 3)
                     },
                     future.globals = FALSE,
                     future.packages = "sleepsimReval") %>%
  bind_rows() %>%
  mutate(iteration_id = iteration_uid) %>%
  left_join(scen, by="iteration_id")
Sys.time() - t1

# Write to disk
saveRDS(em_mean_est, "data-raw/3_preprocess_simulation_results/data_preprocessed/emission_means_est.rds")

# Clean
rm(em_mu_only, em_mean_est)

# Emission distribution variances ----

# Only variance estimates
em_var_only <- map(res, function(x) x$emiss_varmu_bar)

# Retrieve parameter estimates
t1 <- Sys.time()
em_var_est <- future_lapply(em_var_only,
                             function(x) {
                               sleepsimReval:::postprocess_param_est(x, 3)
                               },
                             future.globals = FALSE,
                             future.packages = "sleepsimReval") %>%
  bind_rows() %>%
  mutate(iteration_id = iteration_uid) %>%
  left_join(scen, by="iteration_id")
Sys.time() - t1

# Write to disk
saveRDS(em_var_est, "data-raw/3_preprocess_simulation_results/data_preprocessed/emission_varmu_est.rds")

# Clean
rm(em_var_only, em_var_est)

# Transition probabilities ----

# Only probabilities
gam_int_prob <- map(res, function(x) x$gamma_prob_bar)

# Retrieve parameter estimates
t1 <- Sys.time()
gam_prob_est <- future_lapply(gam_int_prob,
                            function(x) {
                              sleepsimReval:::postprocess_gamma_int(x, 3)
                            },
                            future.globals = FALSE,
                            future.packages = "sleepsimReval") %>%
  bind_rows() %>%
  mutate(iteration_id = iteration_uid) %>%
  left_join(scen, by="iteration_id")
Sys.time() - t1

# Write to disk
saveRDS(gam_prob_est, "data-raw/3_preprocess_simulation_results/data_preprocessed/gamma_int_prob.rds")

# Clean
rm(gam_prob_est, gam_int_prob)

# Extract credible intervals ----

# Emission means ----

# Subset
emiss_mu_cci <- map(res, function(x) x$credible_intervals$emiss_mu_bar)

# Retrieve parameter estimates
t1 <- Sys.time()
emiss_mu_cci_est <- future_lapply(emiss_mu_cci,
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
                              },
                              future.globals = FALSE) %>%
  bind_rows() %>%
  mutate(iteration_id = iteration_uid) %>%
  left_join(scen, by="iteration_id")
Sys.time() - t1

# Write to disk
saveRDS(emiss_mu_cci_est, "data-raw/3_preprocess_simulation_results/data_preprocessed/emiss_mu_cci.rds")

# Clean
rm(emiss_mu_cci, emiss_mu_cci_est)

# Emission variances -----

# Subset
emiss_varmu_cci <- map(res, function(x) x$credible_intervals$emiss_varmu_bar)

# Retrieve parameter estimates
t1 <- Sys.time()
emiss_varmu_cci_est <- future_lapply(emiss_varmu_cci,
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
                                  },
                                  future.globals = FALSE) %>%
  bind_rows() %>%
  mutate(iteration_id = iteration_uid) %>%
  left_join(scen, by="iteration_id")
Sys.time() - t1

# Write to disk
saveRDS(emiss_varmu_cci_est, "data-raw/3_preprocess_simulation_results/data_preprocessed/emiss_varmu_cci.rds")

# Clean
rm(emiss_varmu_cci, emiss_varmu_cci_est)

# Transition probabilities -----

# Subset
gamma_prob_cci <- map(res, function(x) x$credible_intervals$gamma_prob_bar)

# Retrieve parameter estimates
t1 <- Sys.time()
gamma_prob_cci_est <- future_lapply(gamma_prob_cci,
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
                                     },
                                     future.globals = FALSE) %>%
  bind_rows() %>%
  mutate(iteration_id = iteration_uid) %>%
  left_join(scen, by="iteration_id")
Sys.time() - t1

# Write to disk
saveRDS(gamma_prob_cci_est, "data-raw/3_preprocess_simulation_results/data_preprocessed/gamma_prob_cci.rds")

# Clean
rm(gamma_prob_cci, gamma_prob_cci_est)

# Close PSOCK cluster ----

future:::ClusterRegistry("stop")

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
  "true_val" = c(-0.36, -0.6, 0.7, 1.01, -1.31, -.24, .75, -1.31, 0.005),
  stringsAsFactors = FALSE
)
saveRDS(rv, "data-raw/3_preprocess_simulation_results/data_preprocessed/emiss_means_true_mapval.rds")

# True values (tpm)
library(sleepsimR)
opts <- getOption("sleepsimR_simulate")
tpm <- data.frame(
  "transition" = c(
    'S1toS1','S1toS2','S1toS3',
    'S2toS1','S2toS2','S2toS3',
    'S3toS1','S3toS2','S3toS3'
  ),
  "true_val" = as.vector(opts$gamma_bar),
  stringsAsFactors = FALSE
)
saveRDS(tpm, "data-raw/3_preprocess_simulation_results/data_preprocessed/tpm_probs_true_mapval.rds")
