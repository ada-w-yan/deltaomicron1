% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/camostat_model_run.R
\name{run_exp_camostat_shared_pars}
\alias{run_exp_camostat_shared_pars}
\title{mutlstrain fit with some shared parameters}
\usage{
run_exp_camostat_shared_pars(
  save_dir,
  strain,
  PCR = FALSE,
  n_L = 1,
  n_I = 1,
  camostat_only = FALSE,
  different_eclipse = FALSE,
  reparameterise = FALSE,
  k1_endosomal_limit = FALSE,
  separate_pars = c("prob_infection_tmprss2", "prob_infection_endosomal",
    "log10_burst_size"),
  mvr = TRUE,
  length_run = 2,
  run_flag = TRUE
)
}
\arguments{
\item{save_dir}{string ending in "/": directory to save results in}

\item{strain}{character vector.  names of strains.}

\item{PCR}{logical.  whether to include PCR data (TRUE) or only plaque data (FALSE)}

\item{camostat_only}{logical.  Whether to fit camostat data only (TRUE) or but no drug and camostat data (FALSE)}

\item{different_eclipse}{logical.  if TRUE, fit different eclipse phases for the pathways}

\item{k1_endosomal_limit}{logical. if TRUE, 0 <= k1_endosomal/k1_tmprss2 <= 1}

\item{mvr}{logical. if TRUE, use multivariate proposal; if FALSE, use univariate proposal}

\item{length_run}{integer.  If length == 1, short run (for testing purposes); if length == 2, normal run; if length = 3, extra long run.
otherwise, run until convergence}

\item{run_flag}{logial.  if TRUE, run and postprocess; if FALSE, postprocess only}
}
\value{
NULL (results saved to file)
}
\description{
mutlstrain fit with some shared parameters
}
