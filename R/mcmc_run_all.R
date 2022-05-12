#' wrapper for running and postprocessing MCMC chain
#'
#' @param starting_pars NULL for randomly generated start values, or
#' n_replicates by n_temperatures list of vectors of length n_pars,
#' where n_pars is the number of parameters, n_replicates is the number of parallel
#' chains to run to assess convergence, and n_temperatures is the number of
#' temperatures used for parallel tempering.  If no parallel tempering, use
#' an n_replicates  list of vectors of length n_pars
#' @param starting_point_seed list of n_replicates vectors, each of length
#' n_temperatures, used for random number generator for generating starting
#' point for each chain.
#' @param filenames vector of n_replicates filenames in which to store results
#'(without extension)
#' @param run_flags logical vector with elements
#' run: run MCMC if TRUE
#' process: read in a process MCMC chain if TRUE
#' @param mvr logical.  if TRUE, use multivariate proposal, else use univariate proposal
#' @param mvr_init list of parameters for multivariate proposal.  if NULL, use defaults
#' @param specify_parameters closure to specify parameters by generating the parTab
#' data frame.
#' @param generate_data_and_posterior closure to load data and create the function
#' used to evaluate the likelihood.
#' @param mcmcPars list of tuning parameters for MCMC
#' @param mcmc_seed list of length n_replicates, each element being an integer.
#' Used to seed RNG for MCMC.  if NULL, use default
#' @param generate_start_values closure to generate starting parameter values
#' if starting_pars = NULL
#' @inheritParams postprocess_for_plotting
#' @export
run_all <- function(starting_pars = NULL,
                    starting_point_seed,
                    filenames,
                    run_flags = c(run = TRUE, process = TRUE),
                    run_parallel = TRUE,
                    mvr = TRUE,
                    mvr_init = NULL,
                    specify_parameters,
                    get_data,
                    CREATE_POSTERIOR_FUNC,
                    CREATE_PRIOR_FUNC,
                    mcmcPars = gen_mcmcPars(),
                    mcmc_seed = NULL,
                    generate_start_values = function(x, y, z) generate_random_start(x, y),
                    gen_summary_statistics = NULL,
                    calc_residuals = NULL,
                    get_prediction_ci_df = NULL,
                    plot_flags = c(trace = TRUE,
                                   diagnostics= TRUE,
                                   priors = TRUE,
                                   posteriors = TRUE,
                                   bivariate = FALSE),
                    plot_dynamics = NULL,
                    plot_residuals = NULL) {

  if(run_flags[["run"]]) {
    setup_and_run_MCMC(starting_pars,
                       starting_point_seed,
                       filenames,
                       run_parallel,
                       mvr,
                       mvr_init,
                       specify_parameters,
                       get_data,
                       CREATE_POSTERIOR_FUNC,
                       CREATE_PRIOR_FUNC,
                       mcmcPars,
                       mcmc_seed,
                       generate_start_values)
  }

  if(run_flags[["process"]]) {
    filename <- paste0(filenames[1], ".RData")
    postprocessed_list <- postprocess_for_plotting(filename = filename,
                                                   run_parallel = run_parallel,
                                                   gen_summary_statistics = gen_summary_statistics,
                                                   calc_residuals = calc_residuals,
                                                   get_prediction_ci_df = get_prediction_ci_df,
                                                   plot_flags = plot_flags)

    plot_postprocessed_chain(postprocessed_list = postprocessed_list,
                             plot_dynamics = plot_dynamics,
                             plot_residuals = plot_residuals,
                             plot_flags = plot_flags)

  }
}
