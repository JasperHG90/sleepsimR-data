# Documentation for datasets go here

#' Preprocessed EEG and EOG channels that I use in my thesis.
#'
#' This dataset contains 4 EEG and EOG channels that I use to conduct my
#' simulation study and empirical application.
#'
#' @format A data frame with 58947 rows and 6 variables. Each subject has 1.440 occasions
#' (=time points), except subject 20, who has 1.347 occasions.:
#' \describe{
#'   \item{id}{Unique id of the subject}
#'   \item{sleep_state}{Label of the state at occasion t}
#'   \item{EEG_Fpz_Cz_mean_theta}{Value of the logit-transformed EEG channel at occasion t}
#'   \item{EEG_Fpz_Cz_mean_beta}{Value of the logit-transformed EEG channel at occasion ts}
#'   \item{EOG_median_theta}{Value of the logit-transformed EOG channel at occasion t}
#'   \item{EOG_min_beta}{Value of the logit-transformed EOG channel at occasion t}
#' }
#'
#' @seealso The docker app used to download the data can be found here:
#' <https://doi.org/10.5281/zenodo.3826066>. The docker app used to preprocess the data can be found here:
#' <https://doi.org/10.5281/zenodo.3826068>. The script used generate this dataset can be found in the folder
#' "data-raw/1_preprocess_sleep_data" of this R library.
#'
#' @source \url{https://physionet.org/content/sleep-edfx/1.0.0/}
"sleepdata"

#' Summary statistics collected for the 'sleepdata' dataset.
#'
#' This dataset contains summary statistics for the \link[sleepsimRdata]{sleepdata} dataset.
#'
#' @format A data frame with 12 rows and 4 variables:
#' \describe{
#'   \item{variable}{Name of the emission distribution variable.}
#'   \item{sleep_state}{Label of the latent sleep state.}
#'   \item{mmvar}{Mean of the emission distribution.}
#'   \item{vvvar}{Between-subject variance of the mean of the emission distribution.}
#' }
#'
#' @source \url{https://physionet.org/content/sleep-edfx/1.0.0/}
"summary_statistics_sleepdata"

#' Total variance found in each emission distribution.
#'
#' This dataset contains the total variance for the emission distributions found in
#' the \link[sleepsimRdata]{sleepdata} dataset.
#'
#' @format A data frame with 12 rows and 3 variables:
#' \describe{
#'   \item{variable}{Name of the emission distribution variable.}
#'   \item{sleep_state}{Label of the latent sleep state.}
#'   \item{tvar}{Total variance for each emission distribution.}
#' }
#'
#' @source \url{https://physionet.org/content/sleep-edfx/1.0.0/}
"total_variance_sleepdata"

#' Specification of the simulation scenarios (6.912 iterations)
#'
#' This dataset contains simulation settings for 6.912 iterations (approx. 48 iterations
#' per simulation scenario) of my simulation study. It was created using the function
#' \link[sleepsimR]{generate_scenarios} with the default seed. This subset was run first
#' to ensure that the simulation architecture worked as desired and to obtain a variance
#' estimate with which I could calculate the expected MCMC SE and required iterations for
#' my simulation study.
#'
#' @format A data frame with 6.912 rows and 12 variables. :
#' \describe{
#'    \item{n}{int. number of subjects}
#'    \item{n_t}{int. number of observed data points for each subject}
#'    \item{zeta}{float. between-subject variance for the means of the emission distributions}
#'    \item{Q}{float. between-subject variance for the intercepts of the transition-probability matrix}
#'    \item{scenario_id}{string .unique id of the simulation scenario}
#'    \item{rank}{int. iteration number of the simulation scenario}
#'    \item{iteration_id}{string. unique id of the iteration}
#'    \item{dsim_seed}{int. random seed used to generate the data}
#'    \item{model_seed}{int. random seed used to run the mHMM}
#'    \item{start_gamma}{json. initial values for the between-subject transition probability matrix}
#'    \item{start_emiss}{json. intiial values for each of the 3 emission distributions}
#' }
#'
#' @seealso Version 1.3.1 <https://github.com/JasperHG90/sleepsimR-api/releases/tag/v1.3.1> of the
#' sleepsimR-api program contains this dataset to manage the simulation study. For more information,
#' visit <https://github.com/JasperHG90/sleepsimR-api>
"scen_subs"

#' Specification of the simulation scenarios (29.088 iterations)
#'
#' This dataset contains simulation settings for 29.088 iterations (approx. 202 iterations
#' per simulation scenario) of my simulation study. It was created using the function
#' \link[sleepsimR]{generate_scenarios} with the default seed. This subset was run after the
#' \link[sleepsimRdata]{scen_subs} dataset.
#'
#' @format A data frame with 6.912 rows and 12 variables. :
#' \describe{
#'    \item{n}{int. number of subjects}
#'    \item{n_t}{int. number of observed data points for each subject}
#'    \item{zeta}{float. between-subject variance for the means of the emission distributions}
#'    \item{Q}{float. between-subject variance for the intercepts of the transition-probability matrix}
#'    \item{scenario_id}{string .unique id of the simulation scenario}
#'    \item{rank}{int. iteration number of the simulation scenario}
#'    \item{iteration_id}{string. unique id of the iteration}
#'    \item{dsim_seed}{int. random seed used to generate the data}
#'    \item{model_seed}{int. random seed used to run the mHMM}
#'    \item{start_gamma}{json. initial values for the between-subject transition probability matrix}
#'    \item{start_emiss}{json. intiial values for each of the 3 emission distributions}
#' }
#'
#' @seealso Version 1.3.2 <https://github.com/JasperHG90/sleepsimR-api/releases/tag/v1.3.2> of the
#' sleepsimR-api program contains this dataset to manage the simulation study. For more information,
#' visit <https://github.com/JasperHG90/sleepsimR-api>
"scen"

#' Specification of the simulation scenarios (862 iterations)
#'
#' This dataset contains simulation settings for 862 iterations (approx. 6 iterations
#' per simulation scenario) of my simulation study. This dataset contains instructions to execute
#' a second MCMC chain for specific models. The models from the original simulations are then compared
#' to the second chain to evaluate whether or not they have converged to their stable posterior
#' distributions. This dataset is an adaption of the dataset found in \link[sleepsimRdata]{subs} and
#' was created using the R file found in "data-raw/2_make_simulation_scenarios/make_simulation_scenarios.R"
#'
#' @format A data frame with 6.912 rows and 12 variables. :
#' \describe{
#'    \item{n}{int. number of subjects}
#'    \item{n_t}{int. number of observed data points for each subject}
#'    \item{zeta}{float. between-subject variance for the means of the emission distributions}
#'    \item{Q}{float. between-subject variance for the intercepts of the transition-probability matrix}
#'    \item{scenario_id}{string .unique id of the simulation scenario}
#'    \item{rank}{int. iteration number of the simulation scenario}
#'    \item{iteration_id}{string. unique id of the iteration}
#'    \item{dsim_seed}{int. random seed used to generate the data}
#'    \item{model_seed}{int. random seed used to run the mHMM}
#'    \item{start_gamma}{json. initial values for the between-subject transition probability matrix}
#'    \item{start_emiss}{json. intiial values for each of the 3 emission distributions}
#' }
#'
#' @seealso Version 1.3.3 <https://github.com/JasperHG90/sleepsimR-api/releases/tag/v1.3.3> of the
#' sleepsimR-api program contains this dataset to manage the simulation study. For more information,
#' visit <https://github.com/JasperHG90/sleepsimR-api>
"scen_rerun"

#' Specification of the baseline simulation scenarios (2.104 iterations)
#'
#' This dataset contains simulation settings for 2.104 baseline iterations  of my simulation study.
#' It was created using the function \link[sleepsimR]{generate_scenarios} with seed value 912326.
#'
#' @format A data frame with 2.104 rows and 12 variables. :
#' \describe{
#'    \item{n}{int. number of subjects}
#'    \item{n_t}{int. number of observed data points for each subject}
#'    \item{zeta}{float. between-subject variance for the means of the emission distributions}
#'    \item{Q}{float. between-subject variance for the intercepts of the transition-probability matrix}
#'    \item{scenario_id}{string .unique id of the simulation scenario}
#'    \item{rank}{int. iteration number of the simulation scenario}
#'    \item{iteration_id}{string. unique id of the iteration}
#'    \item{dsim_seed}{int. random seed used to generate the data}
#'    \item{model_seed}{int. random seed used to run the mHMM}
#'    \item{start_gamma}{json. initial values for the between-subject transition probability matrix}
#'    \item{start_emiss}{json. intiial values for each of the 3 emission distributions}
#' }
#'
#' @seealso Versions 1.5.1 <https://github.com/JasperHG90/sleepsimR-api/releases/tag/v1.5.1> and
#' 1.5.2 <https://github.com/JasperHG90/sleepsimR-api/releases/tag/v1.5.2> of the
#' sleepsimR-api program contain this dataset to manage the simulation study. For more information,
#' visit <https://github.com/JasperHG90/sleepsimR-api>
"scen_baseline"

#' Parsed simulation results for emission distribution between-subject means.
#'
#' This dataset contains two separate datasets in a list. The first dataset is raw data that contains all
#' parameter estimates obtained from the simulation results (median MAP, SE MAP, lower 95\% CCI and upper
#' 95\% CCI). The second dataset summarizes the raw in terms of the simulation metrics (bias, MSE etc.) for
#' each simulation scenario. To reproduce these datasets, follow the instructions in the
#' "data-raw/3_preprocess_simulation_results" found in this R library.
#'
#' @format A list with two entries:
#' data_preprocessed: A data frame with 323.937 rows and 14 variables. :
#' \describe{
#'    \item{iteration_id}{string. unique id of the iteration}
#'    \item{emiss_var_short}{string. One of 'EEG_mean_beta', 'EOG_median_theta' or 'EOG_min_beta'.}
#'    \item{state}{string. Latent state. One of 'state1' (Awake), 'state2' (NREM) or 'state3' (REM)}
#'    \item{median}{numeric. MAP value (median) of the posterior distribution.}
#'    \item{SE}{numeric. MAP value (SE) of the posterior distribution.}
#'    \item{var}{string. Concatenation of 'emiss_var_short' plus the statistic (median/SE) plus the
#'    state name.}
#'    \item{true_val}{numeric. Value of the population parameter used in the simulation study. These values are
#'    the same as those used in the sleepsimR library. See
#'    <https://github.com/JasperHG90/sleepsimR/blob/master/R/zzz.R>}
#'    \item{lower}{numeric. Lower 95\% CCI for the median MAP value.}
#'    \item{upper}{numeric. Upper 95\% CCI for the median MAP value.}
#'    \item{scenario_id}{string .unique id of the simulation scenario}
#'    \item{n}{int. number of subjects}
#'    \item{n_t}{int. number of observed data points for each subject}
#'    \item{zeta}{float. between-subject variance for the means of the emission distributions}
#'    \item{Q}{float. between-subject variance for the intercepts of the transition-probability matrix}
#' }
#' summary_by_scenario: A data frame with 1.296 rows and 21 variables. In my study, I use 144 scenarios. Given
#' that I use three emission variables, each of which has 3 latent states, each scenario is represented 9 times
#' in this dataset.
#' \describe{
#'    \item{scenario_id}{string .unique id of the simulation scenario}
#'    \item{emiss_var_short}{string. One of 'EEG_mean_beta', 'EOG_median_theta' or 'EOG_min_beta'.}
#'    \item{bias}{bias of the parameter estimates in the scenario relative to the population parameter.}
#'    \item{bias_mcmc_se}{MCMC standard error of bias estimate.}
#'    \item{empirical_se}{Empirical standard error of the parameter estimates in the scenario.}
#'    \item{empirical_se_mcmc_se}{MCMC standard error of the empirical standard error.}
#'    \item{MSE}{MSE of the parameter estimates in the scenario relative to the population parameter.}
#'    \item{MSE_mcmc_se}{MCMC standard error of the MSE estimate.}
#'    \item{coverage}{Percentage of cases in a scenario where the population parameter is included in
#'    the 95\% CCI.}
#'    \item{coverage_mcmc_se}{MCMC standard error of the coverage estimate.}
#'    \item{bias_corr_coverage}{Percentage of cases in a scenario where the average estimated parameter
#'    value is included in the 95\% CCI.}
#'    \item{bias_corr_coverage_mcmc_se}{MCMC standard error of the bias corrected coverage estimate.}
#'    \item{multimodal}{P-value of a test for multimodality in the distribution of parameter estimates.
#'    See \link[multimode]{modetest}.}
#'    \item{modSE}{Model standard error of the parameter estimates in the scenario.}
#'    \item{modSE_mcmc_se}{MCMC standard error of the modSE estimate.}
#'    \item{n}{int. number of subjects}
#'    \item{n_t}{int. number of observed data points for each subject}
#'    \item{zeta}{float. between-subject variance for the means of the emission distributions}
#'    \item{Q}{float. between-subject variance for the intercepts of the transition-probability matrix}
#' }
#'
#' @seealso The simulation was conducted using the following two programs. The first program is the resource
#' manager found at <https://github.com/JasperHG90/sleepsimR-api/releases>. In particular, versions 1.3.1,
#' 1.3.2 and 1.3.3 were used. These versions contain different iterations of the simulation study. For more
#' information, see \link[sleepsimRdata]{scen}, \link[sleepsimRdata]{scen_subs} and
#' \link[sleepsimRdata]{scen_rerun}. The program used to run the simulations can be found here
#' <https://github.com/JasperHG90/sleepsimR-run/releases>. In particular, version 1.3 was used. The
#' simulations were executed on a cluster managed by SURF <https://www.surf.nl/en>. For more information
#' about the architectural design of the simulation study, visit
#' <https://github.com/JasperHG90/sleepsimR-documentation>.
#' For a definition of the simulation metrics, see: Morris, Tim P., Ian R. White, and Michael J. Crowther. "Using simulation studies to evaluate statistical methods." Statistics in medicine 38.11 (2019): 2074-2102.
"simulation_data_emiss_means"

#' Parsed simulation results for emission distribution between-subject variances.
#'
#' This dataset contains two separate datasets in a list. The first dataset is raw data that contains all
#' parameter estimates obtained from the simulation results (median MAP, SE MAP, lower 95\% CCI and upper
#' 95\% CCI). The second dataset summarizes the raw in terms of the simulation metrics (bias, MSE etc.) for
#' each simulation scenario. To reproduce these datasets, follow the instructions in the
#' "data-raw/3_preprocess_simulation_results" found in this R library.
#'
#' @format A list with two entries:
#' data_preprocessed: A data frame with 323.937 rows and 14 variables. :
#' \describe{
#'    \item{iteration_id}{string. unique id of the iteration}
#'    \item{emiss_var_short}{string. One of 'EEG_mean_beta', 'EOG_median_theta' or 'EOG_min_beta'.}
#'    \item{state}{string. Latent state. One of 'state1' (Awake), 'state2' (NREM) or 'state3' (REM)}
#'    \item{median}{numeric. MAP value (median) of the posterior distribution.}
#'    \item{SE}{numeric. MAP value (SE) of the posterior distribution.}
#'    \item{var}{string. Concatenation of 'emiss_var_short' plus the statistic (median/SE) plus the
#'    state name.}
#'    \item{true_val}{numeric. Value of the population parameter used in the simulation study. These values are
#'    the same as those used in the sleepsimR library. See
#'    <https://github.com/JasperHG90/sleepsimR/blob/master/R/zzz.R>}
#'    \item{lower}{numeric. Lower 95\% CCI for the median MAP value.}
#'    \item{upper}{numeric. Upper 95\% CCI for the median MAP value.}
#'    \item{scenario_id}{string .unique id of the simulation scenario}
#'    \item{n}{int. number of subjects}
#'    \item{n_t}{int. number of observed data points for each subject}
#'    \item{zeta}{float. between-subject variance for the means of the emission distributions}
#'    \item{Q}{float. between-subject variance for the intercepts of the transition-probability matrix}
#' }
#' summary_by_scenario: A data frame with 1.296 rows and 21 variables. In my study, I use 144 scenarios. Given
#' that I use three emission variables, each of which has 3 latent states, each scenario is represented 9 times
#' in this dataset.
#' \describe{
#'    \item{scenario_id}{string .unique id of the simulation scenario}
#'    \item{emiss_var_short}{string. One of 'EEG_mean_beta', 'EOG_median_theta' or 'EOG_min_beta'.}
#'    \item{bias}{bias of the parameter estimates in the scenario relative to the population parameter.}
#'    \item{bias_mcmc_se}{MCMC standard error of bias estimate.}
#'    \item{empirical_se}{Empirical standard error of the parameter estimates in the scenario.}
#'    \item{empirical_se_mcmc_se}{MCMC standard error of the empirical standard error.}
#'    \item{MSE}{MSE of the parameter estimates in the scenario relative to the population parameter.}
#'    \item{MSE_mcmc_se}{MCMC standard error of the MSE estimate.}
#'    \item{coverage}{Percentage of cases in a scenario where the population parameter is included in
#'    the 95\% CCI.}
#'    \item{coverage_mcmc_se}{MCMC standard error of the coverage estimate.}
#'    \item{bias_corr_coverage}{Percentage of cases in a scenario where the average estimated parameter
#'    value is included in the 95\% CCI.}
#'    \item{bias_corr_coverage_mcmc_se}{MCMC standard error of the bias corrected coverage estimate.}
#'    \item{multimodal}{P-value of a test for multimodality in the distribution of parameter estimates.
#'    See \link[multimode]{modetest}.}
#'    \item{modSE}{Model standard error of the parameter estimates in the scenario.}
#'    \item{modSE_mcmc_se}{MCMC standard error of the modSE estimate.}
#'    \item{n}{int. number of subjects}
#'    \item{n_t}{int. number of observed data points for each subject}
#'    \item{zeta}{float. between-subject variance for the means of the emission distributions}
#'    \item{Q}{float. between-subject variance for the intercepts of the transition-probability matrix}
#' }
#'
#' @seealso The simulation was conducted using the following two programs. The first program is the resource
#' manager found at <https://github.com/JasperHG90/sleepsimR-api/releases>. In particular, versions 1.3.1,
#' 1.3.2 and 1.3.3 were used. These versions contain different iterations of the simulation study. For more
#' information, see \link[sleepsimRdata]{scen}, \link[sleepsimRdata]{scen_subs} and
#' \link[sleepsimRdata]{scen_rerun}. The program used to run the simulations can be found here
#' <https://github.com/JasperHG90/sleepsimR-run/releases>. In particular, version 1.3 was used. The
#' simulations were executed on a cluster managed by SURF <https://www.surf.nl/en>. For more information
#' about the architectural design of the simulation study, visit
#' <https://github.com/JasperHG90/sleepsimR-documentation>.
#' For a definition of the simulation metrics, see: Morris, Tim P., Ian R. White, and Michael J. Crowther. "Using simulation studies to evaluate statistical methods." Statistics in medicine 38.11 (2019): 2074-2102.
"simulation_data_emiss_varmu"

#' Parsed simulation results for between-subject transition probabilities.
#'
#' This dataset contains two separate datasets in a list. The first dataset is raw data that contains all
#' parameter estimates obtained from the simulation results (median MAP, SE MAP, lower 95\% CCI and upper
#' 95\% CCI). The second dataset summarizes the raw in terms of the simulation metrics (bias, MSE etc.) for
#' each simulation scenario. To reproduce these datasets, follow the instructions in the
#' "data-raw/3_preprocess_simulation_results" found in this R library.
#'
#' @format A list with two entries:
#' data_preprocessed: A data frame with 323.937 rows and 14 variables. :
#' \describe{
#'    \item{iteration_id}{string. unique id of the iteration}
#'    \item{emiss_var_short}{string. One of 'EEG_mean_beta', 'EOG_median_theta' or 'EOG_min_beta'.}
#'    \item{state}{string. Latent state. One of 'state1' (Awake), 'state2' (NREM) or 'state3' (REM)}
#'    \item{median}{numeric. MAP value (median) of the posterior distribution.}
#'    \item{SE}{numeric. MAP value (SE) of the posterior distribution.}
#'    \item{var}{string. Concatenation of 'emiss_var_short' plus the statistic (median/SE) plus the
#'    state name.}
#'    \item{true_val}{numeric. Value of the population parameter used in the simulation study. These values are
#'    the same as those used in the sleepsimR library. See
#'    <https://github.com/JasperHG90/sleepsimR/blob/master/R/zzz.R>}
#'    \item{lower}{numeric. Lower 95\% CCI for the median MAP value.}
#'    \item{upper}{numeric. Upper 95\% CCI for the median MAP value.}
#'    \item{scenario_id}{string .unique id of the simulation scenario}
#'    \item{n}{int. number of subjects}
#'    \item{n_t}{int. number of observed data points for each subject}
#'    \item{zeta}{float. between-subject variance for the means of the emission distributions}
#'    \item{Q}{float. between-subject variance for the intercepts of the transition-probability matrix}
#' }
#' summary_by_scenario: A data frame with 1.296 rows and 21 variables. In my study, I use 144 scenarios. Given
#' that I use three emission variables, each of which has 3 latent states, each scenario is represented 9 times
#' in this dataset.
#' \describe{
#'    \item{scenario_id}{string .unique id of the simulation scenario}
#'    \item{emiss_var_short}{string. One of 'EEG_mean_beta', 'EOG_median_theta' or 'EOG_min_beta'.}
#'    \item{bias}{bias of the parameter estimates in the scenario relative to the population parameter.}
#'    \item{bias_mcmc_se}{MCMC standard error of bias estimate.}
#'    \item{empirical_se}{Empirical standard error of the parameter estimates in the scenario.}
#'    \item{empirical_se_mcmc_se}{MCMC standard error of the empirical standard error.}
#'    \item{MSE}{MSE of the parameter estimates in the scenario relative to the population parameter.}
#'    \item{MSE_mcmc_se}{MCMC standard error of the MSE estimate.}
#'    \item{coverage}{Percentage of cases in a scenario where the population parameter is included in
#'    the 95\% CCI.}
#'    \item{coverage_mcmc_se}{MCMC standard error of the coverage estimate.}
#'    \item{bias_corr_coverage}{Percentage of cases in a scenario where the average estimated parameter
#'    value is included in the 95\% CCI.}
#'    \item{bias_corr_coverage_mcmc_se}{MCMC standard error of the bias corrected coverage estimate.}
#'    \item{multimodal}{P-value of a test for multimodality in the distribution of parameter estimates.
#'    See \link[multimode]{modetest}.}
#'    \item{modSE}{Model standard error of the parameter estimates in the scenario.}
#'    \item{modSE_mcmc_se}{MCMC standard error of the modSE estimate.}
#'    \item{n}{int. number of subjects}
#'    \item{n_t}{int. number of observed data points for each subject}
#'    \item{zeta}{float. between-subject variance for the means of the emission distributions}
#'    \item{Q}{float. between-subject variance for the intercepts of the transition-probability matrix}
#' }
#'
#' @seealso The simulation was conducted using the following two programs. The first program is the resource
#' manager found at <https://github.com/JasperHG90/sleepsimR-api/releases>. In particular, versions 1.3.1,
#' 1.3.2 and 1.3.3 were used. These versions contain different iterations of the simulation study. For more
#' information, see \link[sleepsimRdata]{scen}, \link[sleepsimRdata]{scen_subs} and
#' \link[sleepsimRdata]{scen_rerun}. The program used to run the simulations can be found here
#' <https://github.com/JasperHG90/sleepsimR-run/releases>. In particular, version 1.3 was used. The
#' simulations were executed on a cluster managed by SURF <https://www.surf.nl/en>. For more information
#' about the architectural design of the simulation study, visit
#' <https://github.com/JasperHG90/sleepsimR-documentation>.
#' For a definition of the simulation metrics, see: Morris, Tim P., Ian R. White, and Michael J. Crowther. "Using simulation studies to evaluate statistical methods." Statistics in medicine 38.11 (2019): 2074-2102.
"simulation_data_gamma_prob"

#' Simulation results for the baseline scenarios, zeta=0.25 (scenarios 1-5A)
#'
#' This dataset contains the simulation results for baseline scenarios 1-5A. To reproduce these datasets, follow the instructions in the
#' "data-raw/5_preprocess_baseline_results" found in this R library.
#'
#' @format A data frame with 27 rows and 32 variables. :
#' \describe{
#'    \item{emiss_var_short}{string. One of 'EEG_mean_beta', 'EOG_median_theta' or 'EOG_min_beta'.}
#'    \item{state}{string. Latent state. One of 'state1' (Awake), 'state2' (NREM) or 'state3' (REM) or transitions from state to state.}
#'    \item{pbias_scenario_{1-5}}{String. Percent bias of the parameter estimates versus the population parameter for a scenario. MC error given in parentheses.}
#'    \item{empirical_se_scenario_{1-5}}{String. Empirical SE for a scenario. MC error given in parentheses.}
#'    \item{modSE_scenario_{1-5}}{String. Average model SE/CCI for a scenario. MC error given in parentheses.}
#'    \item{MSE_scenario_{1-5}}{String. MSE for a scenario. MC error given in parentheses.}
#'    \item{coverage_scenario_{1-5}}{String. Percent coverage for a scenario. MC error given in parentheses.}
#'    \item{bias_corr_coverage_scenario_{1-5}}{numeric. Percent bias-corrected coverage for a scenario. MC error given in parentheses.}
#' }
#'
#' @seealso The simulation was conducted using the following two programs. The first program is the resource
#' manager found at <https://github.com/JasperHG90/sleepsimR-api/releases>. In particular, versions 1.5.1 - 1.5.4
#' were used. These versions contain different baseline scenarios. For more
#' information, see \link[sleepsimRdata]{scen_baseline}. The program used to run the simulations can be found here
#' <https://github.com/JasperHG90/sleepsimR-run/releases>. In particular, versions 1.5.1 - 1.5.2 were used. The
#' simulations were executed on a cluster managed by SURF <https://www.surf.nl/en>. For more information
#' about the architectural design of the simulation study, visit
#' <https://github.com/JasperHG90/sleepsimR-documentation>.
#' For a definition of the simulation metrics, see: Morris, Tim P., Ian R. White, and Michael J. Crowther. "Using simulation studies to evaluate statistical methods." Statistics in medicine 38.11 (2019): 2074-2102.
"scen_out_zeta_025"

#' Simulation results for the baseline scenarios, zeta=0.5 (scenarios 1-5B)
#'
#' This dataset contains the simulation results for baseline scenarios 1-5B. To reproduce these datasets, follow the instructions in the
#' "data-raw/5_preprocess_baseline_results" found in this R library.
#'
#' @format A data frame with 27 rows and 32 variables. :
#' \describe{
#'    \item{emiss_var_short}{string. One of 'EEG_mean_beta', 'EOG_median_theta' or 'EOG_min_beta'.}
#'    \item{state}{string. Latent state. One of 'state1' (Awake), 'state2' (NREM) or 'state3' (REM) or transitions from state to state.}
#'    \item{pbias_scenario_{1-5}}{String. Percent bias of the parameter estimates versus the population parameter for a scenario. MC error given in parentheses.}
#'    \item{empirical_se_scenario_{1-5}}{String. Empirical SE for a scenario. MC error given in parentheses.}
#'    \item{modSE_scenario_{1-5}}{String. Average model SE/CCI for a scenario. MC error given in parentheses.}
#'    \item{MSE_scenario_{1-5}}{String. MSE for a scenario. MC error given in parentheses.}
#'    \item{coverage_scenario_{1-5}}{String. Percent coverage for a scenario. MC error given in parentheses.}
#'    \item{bias_corr_coverage_scenario_{1-5}}{numeric. Percent bias-corrected coverage for a scenario. MC error given in parentheses.}
#' }
#'
#' @seealso The simulation was conducted using the following two programs. The first program is the resource
#' manager found at <https://github.com/JasperHG90/sleepsimR-api/releases>. In particular, versions 1.5.1 - 1.5.4
#' were used. These versions contain different baseline scenarios. For more
#' information, see \link[sleepsimRdata]{scen_baseline}. The program used to run the simulations can be found here
#' <https://github.com/JasperHG90/sleepsimR-run/releases>. In particular, versions 1.5.1 - 1.5.2 were used. The
#' simulations were executed on a cluster managed by SURF <https://www.surf.nl/en>. For more information
#' about the architectural design of the simulation study, visit
#' <https://github.com/JasperHG90/sleepsimR-documentation>.
#' For a definition of the simulation metrics, see: Morris, Tim P., Ian R. White, and Michael J. Crowther. "Using simulation studies to evaluate statistical methods." Statistics in medicine 38.11 (2019): 2074-2102.
"scen_out_zeta_025"
