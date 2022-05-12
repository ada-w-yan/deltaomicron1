calc_sol_camostat_preincub <- function(pars, solving_time, with_camostat) {

  model_filename <- "model_n_cell_types.R"

  gen <- compile_model(model_filename)
  par_names_mod <- get_mod_info(model_filename)$pars

  pars$beta <- calc_beta_pathway(pars, ifelse(with_camostat, "endosomal", "both"))

  ## function to extract the number of latent and infectious cells at t = 0
  extract_iv <- function(sol, col_name, tol = 1e-6) {
    cols <- (grepl(paste0(col_name, "\\["),colnames(sol)) | colnames(sol) == col_name)
    iv <- unname(sol[2, cols])
    iv[iv <=tol] <- 0
    iv
  }

  pars_incubation <- pars
  pars_incubation <- pars_incubation[par_names_mod]
  mod <- gen$new(user = pars_incubation)

  sol_incubation <- solve_ODE_error_handling(mod, solving_time)
  sol_incubation
}

#' calculate solution to ODEs for each pathway, and extract median, 95\% CI and max LL
#'
#' @param filename .RData file produced by lazymcmc
#' @param model string: "full", "endosomal_only", "tmprss2_only"
#' @param solving_time_max time to solve until.   Default = 72 hpi.
#' @return tibble with columns model, t, V, T1... and probs, which can take the values 0.025, 0.5, 0.975 or "max_LL"
#'
#' @import dplyr
#' @import odin
calc_sol_chain <- function(filename, model, solving_time_max = 72) {

  if(!(model %in% c("full", "endosomal_only", "tmprss2_only"))) {
    stop("unknown model")
  }

  # load posterior samples
  chain <- get_MCMC_chain(filename)
  n_samples <- 1000 # number of samples to use to calculate 95% CI
  max_LL_params <- get_max_LL_params(filename) # get maximum likelihood parameters
  # sample uniformly from rest of chain
  chain <- chain[seq(1, nrow(chain), length.out = n_samples),] %>%
    as_tibble %>%
    select(names(max_LL_params))# %>%
  # rbind(max_LL_params)

  # times for which to solve ODE
  solving_time <- seq(-max_LL_params[["incubation_period"]], solving_time_max)

  stages <- "n_L" %in% names(max_LL_params)
  different_eclipse <- any(grepl("k1_endosomal", names(max_LL_params)))
  reparameterise <- any(grepl("prob_infection", names(max_LL_params)))
  PCR <- any(grepl("c_tot", names(max_LL_params)))

  if(different_eclipse) {
    if(PCR) {
      if(stages) {
        model_filename <- "model_different_eclipse_stages_dual.R"
      } else {
        stop("not yet implemented")
      }
    } else if(stages) {
      model_filename <- "model_different_eclipse_stages.R"
    } else {
      stop("not yet implemented")
    }
  } else {
    if(PCR) {
      stop("not yet implemented")
    }
    if(stages) {
      model_filename <- ifelse(model == "tmprss2_only", "model_teiv_stages.R", "model_n_cell_types_stages.R")
    } else {
      model_filename <- ifelse(model == "tmprss2_only", "model_teiv.R", "model_n_cell_types.R")
    }
  }

  gen <- compile_model(model_filename)
  par_names_mod <- get_mod_info(model_filename)$pars


  parTab <- load_into_list(filename, "parTab")$parTab
  strain_names <- get_strain_names(parTab$names)
  n_strains <- length(strain_names)
  transform_pars_fn <- ifelse(model == "tmprss2_only" && (!different_eclipse), transform_pars_wrapper_tmprss2_only, transform_pars_wrapper_camostat)
  transform_pars <- transform_pars_fn(parTab)

  calc_sol <- function(pars) {

    pars <- transform_pars(pars)

    calc_sol_strain <- function(strain) {
      names(pars) <- gsub(paste0("_", strain), "", names(pars))

      if(different_eclipse) {
        if(model == "endosomal_only") {
          pars$beta_tmprss2 <- pars$k1_tmprss2 <- 0
        } else if(model == "tmprss2_only") {
          pars$beta_endosomal <- pars$k1_endosomal <- 0
        }
      } else if(model == "tmprss2_only") {
        pars$beta <- pars$beta_tmprss2
      } else {
        pars$beta <- rep(pars$beta_endosomal, 2)

        if(model == "full") {
          pars$beta[1] <- pars$beta[1] + pars$beta_tmprss2
        }
      }
      if(PCR) {
        pars$V_tot0 <- pars$V_0 * pars$RNA_ratio_inoculum
      }

      ## function to extract the number of latent and infectious cells at t = 0

      extract_iv <- function(sol, col_name, dims, tol = 1e-6) {

        cols <- (grepl(paste0(col_name, "\\["),colnames(sol)) | colnames(sol) == col_name)
        iv <- unname(sol[2, cols])

        if(!missing(dims)) {
          iv <- array(iv, dims)
        } else {
          iv <- as.numeric(iv)
        }
        iv[iv <=tol] <- 0
        iv
      }

      pars_incubation <- pars
      incubation_period <- pars$incubation_period

      pars_incubation <- pars_incubation[par_names_mod]
      mod <- gen$new(user = pars_incubation)

      sol_incubation <- solve_ODE_error_handling(mod, c(-incubation_period, 0))

      if((!stages) || (model == "tmprss2_only" && (!different_eclipse))) {
        pars$L_0 <- extract_iv(sol_incubation, "L")
        pars$I_0 <- extract_iv(sol_incubation, "I")
      } else  {
        if(different_eclipse) {
          pars$L_0 <- extract_iv(sol_incubation, "L", c(pars$n_C + 1, pars$n_L))
        } else {
          pars$L_0 <- extract_iv(sol_incubation, "L", c(pars$n_C, pars$n_L))
        }
        pars$I_0 <- extract_iv(sol_incubation, "I", c(pars$n_C, pars$n_I))
      }


      pars$T_0 <- extract_iv(sol_incubation, "T1")

      if(PCR) {
        pars$V_tot0 <- extract_iv(sol_incubation, "V_tot") * pars$prop_wash
      }
      pars$V_0 <- extract_iv(sol_incubation, "V") * pars$prop_wash

      pars <- pars[par_names_mod]

      do.call(mod$set_user,pars)

      solving_time_subset <- solving_time
      solving_time_subset <-
        solving_time_subset[solving_time_subset >= 0]
      sol <- solve_ODE_error_handling(mod, solving_time_subset) %>%
        rbind(sol_incubation[1, ], .)

      sol <- as.data.frame(sol)

      colnames(sol) <- sub("[", ".", colnames(sol), fixed = TRUE) %>%
        sub("]", "", ., fixed = TRUE)
      if(PCR) {
        out_compartments <- c("t", "V", "V_tot")
      } else {
        out_compartments <- c("t", "V")
      }
      sol <- sol[,out_compartments]
      sol$strain <- strain

      sol

    }
    lapply(strain_names, calc_sol_strain) %>%
      bind_rows()
  }

  probs <- c(0.025, .5, .975)

  prctiles <- apply(chain, 1, calc_sol) %>%
    bind_rows %>%
    group_by(t, strain) %>%
    summarise(across(everything(), quantile, probs = probs)) %>%
    mutate(probs = c("lower", "median", "upper"))

  max_LL_trajectory <- calc_sol(max_LL_params) %>%
    mutate(probs = "max_LL")

  prctiles <- prctiles %>%
    bind_rows(max_LL_trajectory) %>%
    mutate(model = model) %>%
    ungroup()

  if(n_strains == 1) {
    prctiles <- prctiles %>% select(-strain)
  }

  prctiles
}

#' for a set of parameters without stages, add stages and recalculate model solution
#'
#' @param model string:"full", "endosomal_only", or "tmprss2_only"
#' @param pars vector of parameter values
#' @param n_L number of latent stages
#' @param n_I number of infectious stages
#' @return data frame outputted by odin
#'
#' @import dplyr
#' @import odin
calc_sol_add_stages <- function(model, pars, n_L, n_I) {

  if(!(model %in% c("full", "endosomal_only", "tmprss2_only"))) {
    stop("unknown model")
  }

  # times for which to solve ODE
  solving_time <- seq(-1, 72)

  model_filename <- ifelse(model == "tmprss2_only", "model_teiv_stages.R", "model_n_cell_types_stages.R")

  gen <- compile_model(model_filename)
  par_names_mod <- get_mod_info(model_filename)$pars


  parTab <- specify_camostat_parameters_fn(n_L = n_L, n_I = n_I)
  parTab_without_stages <- specify_camostat_parameters_fn()
  names(pars) <- parTab_without_stages$names
  pars_w_stages <- double(nrow(parTab))
  names(pars_w_stages) <- parTab$names
  pars_w_stages[names(pars)] <- pars
  pars_w_stages[["n_L"]] <- n_L
  pars_w_stages[["n_I"]] <- n_I

  pars <- pars_w_stages
  transform_pars_fn <- ifelse(model == "tmprss2_only", transform_pars_wrapper_tmprss2_only, transform_pars_wrapper_camostat)
  transform_pars <- transform_pars_fn(parTab)


  pars <- transform_pars(pars)
  if(model == "tmprss2_only") {
    pars$beta <- pars$beta_tmprss2
  } else {
    pars$beta <- rep(pars$beta_endosomal, 2)

    if(model == "full") {
      pars$beta[1] <- pars$beta[1] + pars$beta_tmprss2
    }
  }

  ## function to extract the number of latent and infectious cells at t = 0

  extract_iv <- function(sol, col_name, dims, tol = 1e-6) {
    cols <- (grepl(paste0(col_name, "\\["),colnames(sol)) | colnames(sol) == col_name)
    iv <- unname(sol[2, cols])
    if(!missing(dims)) {
      iv <- array(iv, dims)
    }
    iv[iv <=tol] <- 0
    iv
  }

  pars_incubation <- pars
  incubation_period <- pars$incubation_period
  pars_incubation <- pars_incubation[par_names_mod]

  mod <- gen$new(user = pars_incubation)

  sol_incubation <- solve_ODE_error_handling(mod, c(-incubation_period, 0))

  if(model == "tmprss2_only") {
    pars$L_0 <- extract_iv(sol_incubation, "L")
    pars$I_0 <- extract_iv(sol_incubation, "I")
  } else {
    pars$L_0 <- extract_iv(sol_incubation, "L", c(pars$n_C, pars$n_L))
    pars$I_0 <- extract_iv(sol_incubation, "I", c(pars$n_C, pars$n_I))
  }

  pars$T_0 <- extract_iv(sol_incubation, "T1")

  pars$V_0 <- extract_iv(sol_incubation, "V") * pars$prop_wash

  pars <- pars[par_names_mod]
  do.call(mod$set_user,pars)

  solving_time_subset <- solving_time
  solving_time_subset <-
    solving_time_subset[solving_time_subset >= 0]
  sol <- solve_ODE_error_handling(mod, solving_time_subset) %>%
    rbind(sol_incubation[1, ], .)

  sol <- as.data.frame(sol)

  colnames(sol) <- sub("[", ".", colnames(sol), fixed = TRUE) %>%
    sub("]", "", ., fixed = TRUE)

  sol

}

# note erasing unused options
calc_sol_drugs <- function(pars, solving_time, camostat_vec, amphoB_vec) {

  model_filename <- "model_different_eclipse_stages_dual.R"

  gen <- compile_model(model_filename)
  par_names_mod <- get_mod_info(model_filename)$pars

  solving_time_subset <- solving_time
  solving_time_subset <-
    solving_time_subset[solving_time_subset >= 0]

  calc_sol_inner <- function(camostat, amphoB) {

    if(camostat) {
      pars$beta_tmprss2 <- 0
    }
    if(!amphoB) {
      pars$beta_endosomal <- pars$beta_endosomal * (1 - pars$IFITM)
    }

    ## function to extract the number of latent and infectious cells at t = 0

    extract_iv <- function(sol, col_name, dims, tol = 1e-6) {
      cols <- (grepl(paste0(col_name, "\\["),colnames(sol)) | colnames(sol) == col_name)
      iv <- unname(sol[2, cols])
      if(!missing(dims)) {
        iv <- array(iv, dims)
      }
      iv[iv <=tol] <- 0
      iv
    }

    pars_incubation <- pars
    incubation_period <- pars$incubation_period
    pars_incubation <- pars_incubation[par_names_mod]
    mod <- gen$new(user = pars_incubation)

    sol_incubation <- solve_ODE_error_handling(mod, c(-incubation_period, 0))

    pars$L_0 <- extract_iv(sol_incubation, "L", c(pars$n_C + 1, pars$n_L))
    pars$I_0 <- extract_iv(sol_incubation, "I", c(pars$n_C, pars$n_I))
    pars$T_0 <- extract_iv(sol_incubation, "T1")
    pars$V_tot0 <- extract_iv(sol_incubation, "V_tot") * pars$prop_wash
    pars$V_0 <- extract_iv(sol_incubation, "V") * pars$prop_wash

    pars <- pars[par_names_mod]
    do.call(mod$set_user,pars)


    sol <- solve_ODE_error_handling(mod, solving_time_subset) %>%
      rbind(sol_incubation[1, ], .)

    sol <- as.data.frame(sol)
    sol <- sol[,c("V", "V_tot")]
    colname_suffix <- paste0(ifelse(camostat, "_Camostat", ""), ifelse(amphoB, "_AmphoB", ""))
    colnames(sol) <- paste0(colnames(sol), colname_suffix)

    sol
  }

  sols <- Map(calc_sol_inner, camostat_vec, amphoB_vec)
  solving_time_df <- data.frame(t = c(-pars$incubation_period, solving_time_subset))
  sols <- c(list(solving_time_df), sols) %>%
    do.call(cbind, .)

  sols

}

#' calculate solution to ODEs for each pathway, and extract median, 95\% CI and max LL
#'
#' @param filename .RData file produced by lazymcmc
#' @param model string: "full", "endosomal_only", "tmprss2_only"
#' @param solving_time_max time to solve until.   Default = 72 hpi.
#' @return tibble with columns model, t, V, T1... and probs, which can take the values 0.025, 0.5, 0.975 or "max_LL"
#'
#' @import dplyr
#' @import odin
calc_sol_chain_camostat_amphoB <- function(filename, model, solving_time_max = 72) {

  if(!(model %in% c("full", "endosomal_only", "tmprss2_only", "full_no_IFITM", "endosomal_no_IFITM"))) {
    stop("unknown model")
  }

  # load posterior samples
  chain <- get_MCMC_chain(filename)
  n_samples <- 1000 # number of samples to use to calculate 95% CI
  max_LL_params <- get_max_LL_params(filename) # get maximum likelihood parameters
  # sample uniformly from rest of chain
  chain <- chain[seq(1, nrow(chain), length.out = n_samples),] %>%
    as_tibble %>%
    select(names(max_LL_params))# %>%
  # rbind(max_LL_params)

  # times for which to solve ODE
  solving_time <- seq(-max_LL_params[["incubation_period"]], solving_time_max)

  PCR <- any(grepl("c_tot", names(max_LL_params)))

  if(PCR) {
    model_filename <- "model_different_eclipse_stages_dual.R"
  } else {
    model_filename <- "model_different_eclipse_stages.R"
  }

  gen <- compile_model(model_filename)
  par_names_mod <- get_mod_info(model_filename)$pars

  parTab <- load_into_list(filename, "parTab")$parTab
  transform_pars <- transform_pars_wrapper_camostat(parTab)

  calc_sol <- function(pars) {

    pars <- transform_pars(pars)

    if(grepl("endosomal", model)) {
      pars$beta_tmprss2 <- pars$k1_tmprss2 <- 0
    } else if(model == "tmprss2_only") {
      pars$beta_endosomal <- pars$k1_endosomal <- 0
    }
    if(grepl("no_IFITM", model)) {
      pars$IFITM <- 0
    }
    pars$beta_endosomal <- pars$beta_endosomal * (1 - pars$IFITM)

    if(PCR) {
      pars$V_tot0 <- pars$V_0 * pars$RNA_ratio_inoculum
    }

    ## function to extract the number of latent and infectious cells at t = 0

    extract_iv <- function(sol, col_name, dims, tol = 1e-6) {

      cols <- (grepl(paste0(col_name, "\\["),colnames(sol)) | colnames(sol) == col_name)
      iv <- unname(sol[2, cols])

      if(!missing(dims)) {
        iv <- array(iv, dims)
      } else {
        iv <- as.numeric(iv)
      }
      iv[iv <=tol] <- 0
      iv
    }

    pars_incubation <- pars
    incubation_period <- pars$incubation_period

    pars_incubation <- pars_incubation[par_names_mod]
    mod <- gen$new(user = pars_incubation)

    sol_incubation <- solve_ODE_error_handling(mod, c(-incubation_period, 0))

    pars$L_0 <- extract_iv(sol_incubation, "L", c(pars$n_C + 1, pars$n_L))

    pars$I_0 <- extract_iv(sol_incubation, "I", c(pars$n_C, pars$n_I))

    pars$T_0 <- extract_iv(sol_incubation, "T1")

    if(PCR) {
      pars$V_tot0 <- extract_iv(sol_incubation, "V_tot") * pars$prop_wash
    }
    pars$V_0 <- extract_iv(sol_incubation, "V") * pars$prop_wash

    pars <- pars[par_names_mod]

    do.call(mod$set_user,pars)

    solving_time_subset <- solving_time
    solving_time_subset <-
      solving_time_subset[solving_time_subset >= 0]
    sol <- solve_ODE_error_handling(mod, solving_time_subset) %>%
      rbind(sol_incubation[1, ], .)

    sol <- as.data.frame(sol)

    colnames(sol) <- sub("[", ".", colnames(sol), fixed = TRUE) %>%
      sub("]", "", ., fixed = TRUE)
    if(PCR) {
      out_compartments <- c("t", "V", "V_tot")
    } else {
      out_compartments <- c("t", "V")
    }
    sol <- sol[,out_compartments]

    sol
  }

  probs <- c(0.025, .5, .975)

  prctiles <- apply(chain, 1, calc_sol) %>%
    bind_rows %>%
    group_by(t) %>%
    summarise(across(everything(), quantile, probs = probs)) %>%
    mutate(probs = c("lower", "median", "upper"))

  max_LL_trajectory <- calc_sol(max_LL_params) %>%
    mutate(probs = "max_LL")

  prctiles <- prctiles %>%
    bind_rows(max_LL_trajectory) %>%
    mutate(model = model) %>%
    ungroup()

  prctiles
}

calc_trajectories <- function(filename) {
  models <- c("full", "endosomal_only", "tmprss2_only")
  trajectories <- lapply(models, calc_sol_chain, filename = paste0(dir_name, "1.RData")) %>%
    bind_rows()
}

calc_trajectories_camostat_amphoB <- function(filename) {
  models <- c("full", "endosomal_only", "tmprss2_only", "full_no_IFITM", "endosomal_no_IFITM")
  trajectories <- lapply(models, calc_sol_chain_camostat_amphoB, filename = filename) %>%
    bind_rows()
}
