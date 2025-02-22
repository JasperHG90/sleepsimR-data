% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{scen_subs}
\alias{scen_subs}
\title{Specification of the simulation scenarios (6.912 iterations)}
\format{
A data frame with 6.912 rows and 12 variables. :
\describe{
   \item{n}{int. number of subjects}
   \item{n_t}{int. number of observed data points for each subject}
   \item{zeta}{float. between-subject variance for the means of the emission distributions}
   \item{Q}{float. between-subject variance for the intercepts of the transition-probability matrix}
   \item{scenario_id}{string .unique id of the simulation scenario}
   \item{rank}{int. iteration number of the simulation scenario}
   \item{iteration_id}{string. unique id of the iteration}
   \item{dsim_seed}{int. random seed used to generate the data}
   \item{model_seed}{int. random seed used to run the mHMM}
   \item{start_gamma}{json. initial values for the between-subject transition probability matrix}
   \item{start_emiss}{json. intiial values for each of the 3 emission distributions}
}
}
\usage{
scen_subs
}
\description{
This dataset contains simulation settings for 6.912 iterations (approx. 48 iterations
per simulation scenario) of my simulation study. It was created using the function
\link[sleepsimR]{generate_scenarios} with the default seed. This subset was run first
to ensure that the simulation architecture worked as desired and to obtain a variance
estimate with which I could calculate the expected MCMC SE and required iterations for
my simulation study.
}
\seealso{
Version 1.3.1 <https://github.com/JasperHG90/sleepsimR-api/releases/tag/v1.3.1> of the
sleepsimR-api program contains this dataset to manage the simulation study. For more information,
visit <https://github.com/JasperHG90/sleepsimR-api>
}
\keyword{datasets}
