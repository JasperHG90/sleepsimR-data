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
#' <>. The docker app used to preprocess the data can be found here:
#' <>. The script used generate this dataset can be found here:
#' <>.
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
