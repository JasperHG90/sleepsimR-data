## Compute simulation metrics

## Compute bias, MSE, emperical SE etc.
rm(list=ls())
library(sleepsimReval)
library(tidyr)
library(purrr)
library(stringr)
library(dplyr)

# Base folder
SCEN <- 4
BASE_FOLDER <- "data-raw/5_preprocess_baseline_results/data_preprocessed/"

# Rounder function
rounder <- function(x, digits=2) {
  as.numeric(format(round(x, digits=digits), nsmall=digits))
}

# Get estimates for scenarios 3 and 4
emiss_varmu_est_3 <- readRDS(paste0(BASE_FOLDER, "emission_varmu_est_scen", 3, ".rds")) %>% filter(zeta==0.5) %>% distinct()
emiss_varmu_est_4 <- readRDS(paste0(BASE_FOLDER, "emission_varmu_est_scen", 4, ".rds")) %>% filter(zeta==0.5) %>% distinct() %>%
  slice(1:245)
# True values are in the dataset already (zeta)

# Prep the emission mean data.
prepped_emiss_varmu_data3 <- emiss_varmu_est_3 %>%
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
  filter(statistic %in% c("median")) %>%
  arrange(emiss_var_short, state) %>%
  mutate(pbias = (val - 0.5) / 0.5)

prepped_emiss_varmu_data4 <- emiss_varmu_est_4 %>%
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
  filter(statistic %in% c("median")) %>%
  arrange(emiss_var_short, state) %>%
  mutate(pbias = (val - 0.5) / 0.5)

plot_against <- data.frame(
  x = prepped_emiss_varmu_data3$pbias,
  y = prepped_emiss_varmu_data4$pbias,
  emiss_var_short = prepped_emiss_varmu_data3$emiss_var_short,
  state = prepped_emiss_varmu_data3$state
)

library(ggplot2)
library(ggExtra)
state_map <- c("state1" = "Awake", "state2" = "NREM", "state3" = "REM")
var_map <- c("EEG_mean_beta" = "EEG mean beta", "EOG_median_theta" = "EOG median theta",
             "EOG_min_beta" = "EOG min beta")
ggplot(plot_against %>%
         mutate(state = state_map[state],
                emiss_var_short = var_map[emiss_var_short]),
       aes(x=y,y=x, color=state)) +
  geom_point(alpha=0.4) +
  scale_x_continuous(name = latex2exp::TeX("Percent bias (Scenario 3B, $N=80$)"),
                     limits=c(0,1.5),
                     labels = scales::percent) +
  scale_y_continuous(name = latex2exp::TeX("Percent bias (Scenario 5B, $N=140$)"),
                     limits=c(0,1.5),
                     labels = scales::percent) +
  facet_grid("state ~ emiss_var_short") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_thesis() +
  theme(legend.position = "none",
        axis.text = element_text(size=rel(0.9)))

#%% Plot emp SE and modSE ----

