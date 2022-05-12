#' estimate parameters for the different models in the study
#'
#' @param save_dir string ending in "/": directory to save results in
#' @param strain character vector.  names of strains.
#' @param PCR logical.  whether to include PCR data (TRUE) or only plaque data (FALSE)
#' @param camostat_only logical.  Whether to fit camostat data only (TRUE) or but no drug and camostat data (FALSE)
#' @param different_eclipse logical.  if TRUE, fit different eclipse phases for the pathways
#' @param reparameterise: logical.  if FALSE, use beta, p etc, if TRUE, use probability of infection (beta T_0) / c, burst size, delta, c
#' @param k1_endosomal_limit logical. if TRUE, 0 <= k1_endosomal/k1_tmprss2 <= 1
#' @param mvr logical.  If TRUE, use multivariate proposal, if FALSE, use univariate proposal.
#' @param length_run integer.  If length == 1, short run (for testing purposes); if length == 2, normal run; if length = 3, extra long run.
#' otherwise, run until convergence
#' @param run_flag logial.  if TRUE, run and postprocess; if FALSE, postprocess only
#' @return NULL (results saved to file)
#' @export
run_exp_camostat <- function(save_dir,
                             strain,
                             PCR = FALSE,
                             n_L = 1,
                             n_I = 1,
                             camostat_only = FALSE,
                             different_eclipse = FALSE,
                             reparameterise = FALSE,
                             k1_endosomal_limit = FALSE,
                             mvr = TRUE,
                             length_run = 2,
                             run_flag = TRUE) {

  # number of parallel chains to run
  n_replicates <- 3
  # number of temperatures for parallel tempering
  n_temperatures <- 5
  # set random number generator seed
  seeds <- set_default_seeds(n_replicates, n_temperatures)

  filenames = paste0(save_dir, vapply(seeds, function(x) x[1], double(1)))


  n_replicates <- 3

  CREATE_PRIOR_FUNC <- function(x) {
    function(y) 0
  }

  # define parameters for MCMC
  if(length_run == 1) {
    mcmcPars <- gen_mcmcPars_timing()
  } else if(length_run == 2){
    mcmcPars <- gen_mcmcPars(temperature = seq_len(n_temperatures), parallel_tempering_iter = 10)
  } else if (length_run == 3) {
    mcmcPars <- gen_mcmcPars_long(temperature = seq_len(n_temperatures), parallel_tempering_iter = 10)
  } else {
    stop("unknown length run")
  }

  # define function to generate random starting values
  generate_start_values <- function(x, y, z, f) generate_random_start(x, y, f)

  n_strains <- length(strain)

  assay <- vcapply(strain, function(x) ifelse(x %in% c("delta", "omicron"), "plaque", "tcid50"))
  assay <- unique(assay)

  if(length(assay) > 1) {
    stop("cannot simultaneously fit datasets with TCID50 and plaque assays")
  }

  if(camostat_only) {
    if(PCR) {
      stop("not yet implemented")
    }
    if(n_strains == 1) {
      get_data <- function() {
        read_camostat_data(strain, PCR = PCR) %>%
          select(-V_no_drug)
      }
    } else {
      get_data <- function() {
        data_df <- lapply(strain, read_camostat_data(x, PCR = PCR))
        data_df <- lapply(data_df, function(x) select(x, -V_no_drug))
        names(data_df) <- strain
        data_df
      }
    }
    specify_parameters <- function() {
      parTab <- specify_camostat_parameters_fn(n_L = n_L, n_I = n_I, assay = assay,
                                               different_eclipse = different_eclipse,
                                               reparameterise = reparameterise,
                                               k1_endosomal_limit = k1_endosomal_limit,
                                               strain_names = strain)
      if(reparameterise) {
        stop("not yet implemented")
      }
      parTab <- parTab[parTab$names != "log10_beta_tmprss2",]
      if(different_eclipse) {
        parTab <- parTab[parTab$names != "k1_tmprss2",]
      }
      parTab
    }

  } else {
    if(n_strains == 1) {
      get_data <- function() read_camostat_data(strain, PCR = PCR)
    } else {
      get_data <- function()  {
        data_df <- lapply(strain, read_camostat_data, PCR = PCR)
        names(data_df) <- strain
        data_df
      }

    }
    specify_parameters <- function() specify_camostat_parameters_fn(n_L = n_L, n_I = n_I, assay = assay,
                                                                    different_eclipse = different_eclipse,
                                                                    reparameterise = reparameterise,
                                                                    k1_endosomal_limit = k1_endosomal_limit,
                                                                    strain_names = strain, PCR = PCR)


  }



  if(camostat_only) {
    gen_summary_statistics_fn <- NULL
  } else {
    gen_summary_statistics_fn <- gen_summary_statistics_fn_camostat
  }

  # run MCMC
  run_all(starting_point_seed = seeds,
          filenames = filenames,
          run_parallel = TRUE,
          mvr = mvr,
          specify_parameters = specify_parameters,
          get_data = get_data,
          CREATE_POSTERIOR_FUNC = function(x, y, z) CREATE_POSTERIOR_FUNC_fn_camostat(x, y, z, camostat_only),
          CREATE_PRIOR_FUNC = CREATE_PRIOR_FUNC,
          mcmcPars = mcmcPars,
          generate_start_values = generate_start_values,
          gen_summary_statistics = gen_summary_statistics_fn,
          calc_residuals = NULL,
          plot_flags = c(trace = FALSE,
                         diagnostics= FALSE,
                         priors = FALSE,
                         posteriors = FALSE,
                         bivariate = FALSE),
          plot_dynamics = plot_dynamics_camostat,
          plot_residuals = NULL,
          run = c("run" = run_flag, "process" = TRUE))
}

#' estimate parameters for the different models in the study
#'
#' @param save_dir string ending in "/": directory to save results in
#' @param strain character vector.  names of strains.
#' @param PCR logical.  whether to include PCR data (TRUE) or only plaque data (FALSE)
#' @param camostat_vec logical vector of length n, where n is the number of conditions.  Each value
#' indicates whether camostat was used in the condition.
#' @param amphoB_vec logical vector of length n, where n is the number of conditions.  Each value
#' indicates whether amphoB was used in the condition.
#' @param Calu3 logical. if TRUE, fit model to Calu3 data, if FALSE, use HAE data
#' @param mvr logical.  If TRUE, use multivariate proposal, if FALSE, use univariate proposal.
#' @param length_run integer.  If length == 1, short run (for testing purposes); if length == 2, normal run; if length = 3, extra long run.
#' otherwise, run until convergence
#' @param run_flag logial.  if TRUE, run and postprocess; if FALSE, postprocess only
#' @return NULL (results saved to file)
#' @export
run_exp_camostat_amphoB <- function(save_dir,
                                    strain,
                                    camostat_vec,
                                    amphoB_vec,
                                    Calu3 = FALSE,
                                    PCR = FALSE,
                                    mvr = TRUE,
                                    length_run = 2,
                                    run_flag = TRUE) {

  # number of parallel chains to run
  n_replicates <- 3
  # number of temperatures for parallel tempering
  n_temperatures <- 5
  # set random number generator seed
  seeds <- set_default_seeds(n_replicates, n_temperatures)

  filenames = paste0(save_dir, vapply(seeds, function(x) x[1], double(1)))


  n_replicates <- 3

  CREATE_PRIOR_FUNC <- function(x) {
    function(y) 0
  }

  # define parameters for MCMC
  if(length_run == 1) {
    mcmcPars <- gen_mcmcPars_timing()
  } else if(length_run == 2){
    mcmcPars <- gen_mcmcPars(temperature = seq_len(n_temperatures), parallel_tempering_iter = 10)
  } else if (length_run == 3) {
    mcmcPars <- gen_mcmcPars_long(temperature = seq_len(n_temperatures), parallel_tempering_iter = 10)
  } else {
    stop("unknown length run")
  }

  # define function to generate random starting values
  generate_start_values <- function(x, y, z, f) generate_random_start(x, y, f)

  if(Calu3) {
    get_data <- function() read_Calu3_data(strain, PCR = PCR)
  } else {
    get_data <- function() read_camostat_data(strain, PCR = PCR)
  }

  specify_parameters <- function() specify_camostat_amphoB_parameters_fn(strain, Calu3 = Calu3, PCR = PCR)

  gen_summary_statistics_fn <- gen_summary_statistics_fn_camostat#_amphoB
  CREATE_POSTERIOR_FUNC_fn <- function(x, y, z) CREATE_POSTERIOR_FUNC_fn_camostat_amphoB(x, y, z, camostat_vec, amphoB_vec)

  # run MCMC
  run_all(starting_point_seed = seeds,
          filenames = filenames,
          run_parallel = TRUE,
          mvr = mvr,
          specify_parameters = specify_parameters,
          get_data = get_data,
          CREATE_POSTERIOR_FUNC = CREATE_POSTERIOR_FUNC_fn,
          CREATE_PRIOR_FUNC = CREATE_PRIOR_FUNC,
          mcmcPars = mcmcPars,
          generate_start_values = generate_start_values,
          gen_summary_statistics = gen_summary_statistics_fn,
          calc_residuals = NULL,
          plot_flags = c(trace = FALSE,
                         diagnostics= FALSE,
                         priors = FALSE,
                         posteriors = FALSE,
                         bivariate = FALSE),
          plot_dynamics = plot_dynamics_camostat,
          plot_residuals = NULL,
          run = c("run" = run_flag, "process" = TRUE))
}

#' estimate prob_infect_tmprss2 and prob_infect_camostat, fix other parameters
#'
#' @inheritParams specify_camostat_parameters_fn_fit_subset_pars
#' @inheritParams run_exp_camostat
#' @param mvr logical. if TRUE, use multivariate proposal; if FALSE, use univariate proposal
#' @return NULL (results saved to file)
#' @export
run_exp_camostat_fit_subset_pars <- function(save_dir,
                                             strain,
                                             old_values,
                                             old_parTab,
                                             subset_pars = c("prob_infection_tmprss2", "prob_infection_endosomal"),
                                             mvr = TRUE,
                                             length_run = 2,
                                             run_flag = TRUE) {

  # number of parallel chains to run
  n_replicates <- 3
  # number of temperatures for parallel tempering
  n_temperatures <- 5
  # set random number generator seed
  seeds <- set_default_seeds(n_replicates, n_temperatures)

  filenames = paste0(save_dir, vapply(seeds, function(x) x[1], double(1)))


  n_replicates <- 3

  CREATE_PRIOR_FUNC <- function(x) {
    function(y) 0
  }

  # define parameters for MCMC
  if(length_run == 1) {
    mcmcPars <- gen_mcmcPars_timing()
  } else if(length_run == 2){
    mcmcPars <- gen_mcmcPars(temperature = seq_len(n_temperatures), parallel_tempering_iter = 10)
  } else if (length_run == 3) {
    mcmcPars <- gen_mcmcPars_long(temperature = seq_len(n_temperatures), parallel_tempering_iter = 10)
  } else {
    stop("unknown length run")
  }

  # define function to generate random starting values
  generate_start_values <- function(x, y, z, f) generate_random_start(x, y, f)

  get_data <- function() read_camostat_data(strain)
  specify_parameters <- function() specify_camostat_parameters_fn_fit_subset_pars(old_values, old_parTab, subset_pars)

  gen_summary_statistics_fn <- gen_summary_statistics_fn_camostat

  # run MCMC
  run_all(starting_point_seed = seeds,
          filenames = filenames,
          run_parallel = TRUE,
          mvr = mvr,
          specify_parameters = specify_parameters,
          get_data = get_data,
          CREATE_POSTERIOR_FUNC = function(x, y, z) CREATE_POSTERIOR_FUNC_fn_camostat(x, y, z, FALSE),
          CREATE_PRIOR_FUNC = CREATE_PRIOR_FUNC,
          mcmcPars = mcmcPars,
          generate_start_values = generate_start_values,
          gen_summary_statistics = NULL,#gen_summary_statistics_fn,
          calc_residuals = NULL,
          plot_flags = c(trace = FALSE,
                         diagnostics= FALSE,
                         priors = FALSE,
                         posteriors = FALSE,
                         bivariate = FALSE),
          plot_dynamics = plot_dynamics_camostat,
          plot_residuals = NULL,
          run = c("run" = run_flag, "process" = TRUE))
}

#' mutlstrain fit with some shared parameters
#'
#' @inheritParams specify_camostat_parameters_fn_fit_subset_pars
#' @inheritParams run_exp_camostat
#' @param mvr logical. if TRUE, use multivariate proposal; if FALSE, use univariate proposal
#' @return NULL (results saved to file)
#' @export
run_exp_camostat_shared_pars <- function(save_dir,
                                         strain,
                                         PCR = FALSE,
                                         n_L = 1,
                                         n_I = 1,
                                         camostat_only = FALSE,
                                         different_eclipse = FALSE,
                                         reparameterise = FALSE,
                                         k1_endosomal_limit = FALSE,
                                         separate_pars = c("prob_infection_tmprss2", "prob_infection_endosomal", "log10_burst_size"),
                                         mvr = TRUE,
                                         length_run = 2,
                                         run_flag = TRUE) {

  if(!all(separate_pars %in% c("prob_infection_tmprss2", "prob_infection_endosomal", "log10_burst_size"))) {
    stop("not yet implemented")
  }

  # number of parallel chains to run
  n_replicates <- 3
  # number of temperatures for parallel tempering
  n_temperatures <- 5
  # set random number generator seed
  seeds <- set_default_seeds(n_replicates, n_temperatures)

  filenames = paste0(save_dir, vapply(seeds, function(x) x[1], double(1)))


  n_replicates <- 3

  CREATE_PRIOR_FUNC <- function(x) {
    function(y) 0
  }

  # define parameters for MCMC
  if(length_run == 1) {
    mcmcPars <- gen_mcmcPars_timing()
  } else if(length_run == 2){
    mcmcPars <- gen_mcmcPars(temperature = seq_len(n_temperatures), parallel_tempering_iter = 10)
  } else if (length_run == 3) {
    mcmcPars <- gen_mcmcPars_long(temperature = seq_len(n_temperatures), parallel_tempering_iter = 10)
  } else {
    stop("unknown length run")
  }

  # define function to generate random starting values
  generate_start_values <- function(x, y, z, f) generate_random_start(x, y, f)

  assay <- vcapply(strain, function(x) ifelse(x %in% c("delta", "omicron"), "plaque", "tcid50"))
  assay <- unique(assay)

  if(length(assay) > 1) {
    stop("cannot simultaneously fit datasets with TCID50 and plaque assays")
  }

  n_strains <- length(strain)

  if(n_strains == 1) {
    get_data <- function() read_camostat_data(strain, PCR = PCR)
  } else {
    get_data <- function()  {
      data_df <- lapply(strain, read_camostat_data, PCR = PCR)
      names(data_df) <- strain
      data_df
    }

  }

  specify_parameters <- function() specify_camostat_parameters_fn_share_subset_pars(n_L = n_L,
                                                                                    n_I = n_I,
                                                                                    different_eclipse = different_eclipse,
                                                                                    assay = "plaque",
                                                                                    PCR = PCR,
                                                                                    reparameterise = reparameterise,
                                                                                    k1_endosomal_limit = k1_endosomal_limit,
                                                                                    strain_names = strain,
                                                                                    separate_pars = separate_pars)

  gen_summary_statistics_fn <- gen_summary_statistics_fn_camostat

  # run MCMC
  run_all(starting_point_seed = seeds,
          filenames = filenames,
          run_parallel = TRUE,
          mvr = mvr,
          specify_parameters = specify_parameters,
          get_data = get_data,
          CREATE_POSTERIOR_FUNC = function(x, y, z) CREATE_POSTERIOR_FUNC_fn_camostat(x, y, z, FALSE),
          CREATE_PRIOR_FUNC = CREATE_PRIOR_FUNC,
          mcmcPars = mcmcPars,
          generate_start_values = generate_start_values,
          gen_summary_statistics = gen_summary_statistics_fn,
          calc_residuals = NULL,
          plot_flags = c(trace = FALSE,
                         diagnostics= FALSE,
                         priors = FALSE,
                         posteriors = FALSE,
                         bivariate = FALSE),
          plot_dynamics = plot_dynamics_camostat,
          plot_residuals = NULL,
          run = c("run" = run_flag, "process" = TRUE))
}
