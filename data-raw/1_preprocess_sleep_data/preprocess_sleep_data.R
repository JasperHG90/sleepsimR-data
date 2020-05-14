# Preprocess the sleep dataset obtained by running the docker program
# Written by: Jasper Ginn

library(dplyr)
library(tidyr)

# Select these channels
selected_channels <- c("EEG_Fpz_Cz_mean_theta","EOG_min_beta", "EOG_median_theta", "EEG_Fpz_Cz_mean_beta")
# Load data
ait <- readRDS("data-raw/1_preprocess_sleep_data/EEG_data_final.rds") %>%
  filter(age >= 20 & age < 50) %>%
  # Remove these variables
  select(-patient, -gender, -age) %>%
  # From wide to long --> var == channels, val == value
  gather(var, val, -identifier, -sleep_state, -epoch) %>%
  # Group by person, variable and state
  group_by(identifier, var, sleep_state) %>%
  # Perform logit transformation
  mutate(val = log(val / (1-val))) %>%
  ungroup() %>%
  # Group by person and variable
  group_by(identifier, var) %>%
  # Grand-mean center variables
  mutate(val = ((val - mean(val)) / sd(val))) %>%
  ungroup()
# Filter the data
io <- ait %>%
  filter(var %in% selected_channels) %>%
  # Pivot to wide
  pivot_wider(id_cols = c(identifier, epoch, sleep_state), names_from = var, values_from = val) %>%
  # Remove epochs
  select(-epoch)

# Get unique identifiers
# Map to numeric 1-len(identifer))
uid <- unique(as.character(io$identifier))
ids <- 1:length(uid)
names(ids) <- uid
io$id <- unname(ids[as.character(io$identifier)])
io$identifier <- NULL
io <- io[,c(6, 1, 2 ,3, 4, 5)]
# Rename
sleepdata <- io

# Get summary statistics for each dep var
summary_statistics_sleepdata <- io %>%
  gather(variable, value, -id, -sleep_state) %>%
  group_by(id, variable, sleep_state) %>%
  summarize(mvar = mean(value)) %>%
  ungroup() %>%
  group_by(variable, sleep_state) %>%
  summarize(mmvar = mean(mvar),
            vvvar = var(mvar))
# Total variance
total_variance_sleepdata <- io %>%
  gather(variable, value, -id, -sleep_state) %>%
  group_by(variable, sleep_state) %>%
  summarize(tvar = var(value))

# Save datasets
usethis::use_data(summary_statistics_sleepdata)
usethis::use_data(total_variance_sleepdata)
usethis::use_data(sleepdata)
