#' setup and run MCMC
#'
#' @inheritParams run_all
setup_and_run_MCMC <-function(starting_pars = NULL,
                              starting_point_seed,
                              filenames,
                              run_parallel = TRUE,
                              mvr = TRUE,
                              mvr_init = NULL,
                              specify_parameters,
                              get_data,
                              CREATE_POSTERIOR_FUNC,
                              CREATE_PRIOR_FUNC,
                              mcmcPars = gen_mcmcPars(),
                              mcmc_seed = NULL,
                              generate_start_values = function(x, y, z) generate_random_start(x, y)){

  lapply(filenames, mkdir_from_filename)
  ###############################################################################
  ## setting up random number generation
  ###############################################################################

  ## provide some seeds to generate starting points, if not given
  if(is.null(names(mcmcPars))) {
    n_temperatures <- length(mcmcPars[[1]][["temperature"]])
  } else {
    n_temperatures <- length(mcmcPars[["temperature"]])
  }

  n_replicates <- length(filenames)
  parallel_tempering_flag <- (n_temperatures > 1)

  if(missing(starting_point_seed)){
    starting_point_seed <- make_seed(n_replicates, n_temperatures)
  }

  ## if seeds manually given, check that the number of seeds is the same as the
  ## number of filenames
  if(parallel_tempering_flag){

    tryCatch(invisible(vapply(starting_point_seed, identity, double(n_temperatures))),
             error = function(c) "length of temperatures not equal to length of seeds"
    )

  }

  if(is.null(starting_pars) && n_replicates != length(starting_point_seed)){
    stop("length of filenames not equal to length of seeds")
  }

  ###############################################################################
  ## specifying model parameters
  ###############################################################################

  parTab <- specify_parameters()

  ################################################################################
  ## generating data and defining likelihood function
  ################################################################################

  data_df <- get_data()

  PRIOR_FUNC <- CREATE_PRIOR_FUNC(parTab)
  f <- CREATE_POSTERIOR_FUNC(parTab, data_df, PRIOR_FUNC)

  # check that likelihood function works

  check_likelihood <- function(f, values) {
    f_output <- f(parTab$values)
    lik <- f_output$lik

    stopifnot(is.numeric(lik) && length(lik) == 1 && !is.na(lik) && lik != -Inf)
  }

  check_likelihood(f, parTab$values)

  ################################################################################
  ## generating random starting values
  ################################################################################

  if(is.null(starting_pars)) {
    if(parallel_tempering_flag){
      startTab <- lapply(starting_point_seed,
                         function(y) lapply(y, function(x) generate_start_values(parTab, x, data_df, f)))
    } else {
      startTab <- lapply(starting_point_seed, function(x) generate_start_values(parTab, x, data_df, f))
    }
  } else {

    assign_starting_pars <- function(parTab, starting_pars) {
      startTab <- parTab
      startTab$values <- starting_pars
      startTab
    }

    if(parallel_tempering_flag){
      startTab <- lapply(starting_pars,
                         function(y) lapply(y, function(x) assign_starting_pars(parTab, x)))
    } else {
      startTab <- lapply(starting_pars, function(x) assign_starting_pars(parTab, x))
    }
  }

  if(mvr){

    if(is.null(mvr_init)) {
      define_initial_mvrPars <- function(parTab) {
        n_row_covMat <- sum(parTab$fixed == 0)
        covMat <- diag(nrow(parTab))
        list(covMat,2.38/sqrt(n_row_covMat),w=0.8)
      }

      mvrPars <- define_initial_mvrPars(parTab)

      if(parallel_tempering_flag){
        mvrPars <- rep(list(mvrPars), n_temperatures)
      }

      mvrPars <- rep(list(mvrPars), n_replicates)
    } else {
      mvrPars <- mvr_init
    }

  } else {
    mvrPars <- NULL
  }

  save_all(filenames[1])

  ################################################################################
  ## running MCMC
  ################################################################################
  message("setup complete")
  posterior_list <- lazymcmc::run_MCMC_loop(startTab, data_df, mcmcPars, filenames,
                                            CREATE_POSTERIOR_FUNC, mvrPars, PRIOR_FUNC,
                                            seed = mcmc_seed,
                                            run_parallel = run_parallel)

  list2here(posterior_list, var_names = c("diagnostics", "output"))
  rm(posterior_list)

  save_all(filenames[1])
  invisible(NULL)
}

#' continue running MCMC
#'
#'
cont_MCMC <-function(old_dir,
                     new_filenames,
                     length_run) {
  old_list <- readRDS(paste0(old_dir, "1_output.rds"))
  old_input <- load_into_list(paste0(old_dir, "1.RData"))
  starting_pars <- lapply(old_list, function(x) x$current_pars)
  n_pars <- length(starting_pars[[1]][[1]])
  unfixed_pars <- old_input$parTab$fixed == 0
  if(old_input$mvr) {
    mvr_init <- rep(list(vector("list", old_input$n_temperatures)), old_input$n_replicates)
    covMat <- matrix(0, nrow = n_pars, ncol = n_pars)
    for(i in seq_len(old_input$n_replicates)) {
      for(j in seq_len(old_input$n_temperatures)) {
        covMat[unfixed_pars, unfixed_pars] <- old_list[[i]]$covMat[[j]]
        mvr_init[[i]][[j]] <- list(covMat = covMat, scale = old_list[[i]]$scale[[j]], w = 0.8)
      }
    }
  } else {
    mvr_init <- NULL
  }

  # define parameters for MCMC
  if(length_run == 1) {
    mcmcPars <- gen_mcmcPars_timing()
  } else if(length_run == 2){
    mcmcPars <- gen_mcmcPars(temperature = seq_len(old_input$n_temperatures), parallel_tempering_iter = 10)
  } else if (length_run == 3) {
    mcmcPars <- gen_mcmcPars_long(temperature = seq_len(old_input$n_temperatures), parallel_tempering_iter = 10)
  } else {
    stop("unknown length run")
  }

  mcmcPars$adaptive_period <- mcmcPars$max_adaptive_period <- 0
  mcmcPars <- rep(list(mcmcPars), old_input$n_replicates)
  for(i in seq_len(old_input$n_replicates)) {
    mcmcPars[[i]]$temperature <- old_list[[i]]$temperatures
  }

  setup_and_run_MCMC(starting_pars = starting_pars,
                     filenames = new_filenames,
                     run_parallel = FALSE,#old_input$run_parallel,
                     mvr = old_input$mvr,
                     mvr_init = mvr_init,
                     specify_parameters = old_input$specify_parameters,
                     get_data = old_input$get_data,
                     CREATE_POSTERIOR_FUNC = old_input$CREATE_POSTERIOR_FUNC,
                     CREATE_PRIOR_FUNC = old_input$CREATE_PRIOR_FUNC,
                     mcmcPars = mcmcPars)
}
