# Make baseline simulation scenarios for my simulation study
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
  mutate(n = 140) %>%
  mutate(iteration_id = paste0(iteration_id, "_largeN"))

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
write.csv(scen_baseline_n80nt3200, file=gzfile("data-raw/4_scenarios_baseline/scenarios_baseline_n80nt3200.csv.gz"))
write.csv(scen_baseline_n80nt800, file=gzfile("data-raw/4_scenarios_baseline/scenarios_baseline_n80nt800.csv.gz"))
write.csv(scen_baseline_n140, file=gzfile("data-raw/4_scenarios_baseline/scenarios_baseline_n140nt800.csv.gz"))
write.csv(scen_baseline_n40, file=gzfile("data-raw/4_scenarios_baseline/scenarios_baseline_n40nt800.csv.gz"))

# Save as dataset in this library
usethis::use_data(scen_baseline)
