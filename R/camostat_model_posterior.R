
#'creates a closure to calculate the likelihood and model predictions for given
#' parameter values
#'
#' @return closure to calculate likelihood and model predictions for given
#' parameter values
CREATE_POSTERIOR_FUNC_fn_camostat <-
  function(parTab, data_df, PRIOR_FUNC, camostat_only) {

    transform_pars <- transform_pars_wrapper_camostat(parTab)

    if(camostat_only) {
      prediction_compartments <- "V_Camostat"
    } else {
      prediction_compartments <- c("V_no_drug", "V_Camostat")
    }
    # function to calculate log likelihood
    f <- function(pars) {
      if(PRIOR_FUNC(pars) == -Inf) {
        return(list(lik = -Inf, misc = NA))
      }

      transformed_pars <- transform_pars(pars)

      stages <- "n_L" %in% names(transformed_pars)
      PCR <- any(grepl("c_tot", names(transformed_pars)))
      different_eclipse <- any(grepl("k1_endosomal", names(transformed_pars)))
      grep_str <- "beta_tmprss2"
      beta_strs <- grep(grep_str, names(transformed_pars), value = TRUE)
      n_strains <- length(beta_strs)
      strain_names <- gsub(paste0(grep_str, "_"), "", beta_strs)

      if(camostat_only) {
        prediction_compartments <- "V_Camostat"
      } else {
        prediction_compartments <- c("V_no_drug", "V_Camostat")
      }

      if(PCR) {
        prediction_compartments <- c(prediction_compartments,
                                     gsub("V_", "V_tot_", prediction_compartments))
      }

      sigma <- transformed_pars[["sigma"]]
      obs_threshold <- transformed_pars[["obs_threshold"]]
      if(PCR) {
        sigma_tot <- transformed_pars[["sigma_tot"]]
      }

      if(PCR) {
        sigma_vec <- rep(c(sigma, sigma_tot), each = 2)
      } else {
        sigma_vec <- rep(sigma, 2)
      }
      obs_threshold_vec <- rep(obs_threshold, length(sigma_vec)) # no PCR data points below threshold so can set to whatever

      if(n_strains == 1) {
        solving_time <- make_solving_time(data_df$t)
        if(PCR) {
          transformed_pars$V_tot0 <- transformed_pars$V_0 * transformed_pars$RNA_ratio_inoculum
        }
        sol <- calc_sol_camostat(transformed_pars, solving_time, camostat_only, stages, different_eclipse)
        sol <- sol[,c("t", prediction_compartments)]

        individual_lik <- calc_LL_log10normal_threshold_wrapper(
          data_df,
          prediction_compartments,
          sol,
          prediction_compartments,
          sigma_vec,
          obs_threshold_vec)

        misc <- sol %>%
          as.list %>%
          lapply(name_predictions) %>%
          unlist


        lik <- sum(individual_lik)
      } else {
        calc_sol_inner <- function(strain_name) {
          # find times at which to solve ODEs
          solving_time <- make_solving_time(data_df[[strain_name]]$t)

          names(transformed_pars) <- gsub(paste0("_", strain_name), "", names(transformed_pars))
          names(transformed_pars) <- gsub(paste0(strain_name, "_"), "", names(transformed_pars))
          if(PCR) {
            transformed_pars$V_tot0 <- transformed_pars$V_0 * transformed_pars$RNA_ratio_inoculum
          }

          sol <- calc_sol_camostat(transformed_pars, solving_time, camostat_only, stages, different_eclipse)
          sol <- sol[,c("t", prediction_compartments)]

          individual_lik <- calc_LL_log10normal_threshold_wrapper(
            data_df[[strain_name]],
            prediction_compartments,
            sol,
            prediction_compartments,
            sigma_vec,
            obs_threshold_vec)

          sol <- sol[,prediction_compartments]
          misc <- sol %>%
            as.list %>%
            lapply(name_predictions) %>%
            unlist
          if(PCR) {
            names(misc) <- gsub("V_no_drug", paste0("V_", strain_name, "_no_drug"), names(misc))
            names(misc) <- gsub("V_Camostat", paste0("V_", strain_name, "_Camostat"), names(misc))
            names(misc) <- gsub("V_tot", paste0("V_tot_", strain_name), names(misc))
          } else {
            names(misc) <- gsub("V", paste0("V_", strain_name), names(misc))
          }


          lik <- sum(individual_lik)
          list(lik = lik, misc = misc)
        }

        sols <- lapply(strain_names, calc_sol_inner)

        lik <- lapply(sols, function(x) x$lik) %>% unlist %>% sum
        misc <- lapply(sols, function(x) x$misc) %>% unlist

      }

      list(lik = lik, misc = misc)

    }
    f
  }


#'creates a closure to calculate the likelihood and model predictions for given
#' parameter values
#'
#' @return closure to calculate likelihood and model predictions for given
#' parameter values
CREATE_POSTERIOR_FUNC_fn_camostat_amphoB <-
  function(parTab, data_df, PRIOR_FUNC, camostat_vec, amphoB_vec) {

    transform_pars <- transform_pars_wrapper_camostat(parTab)
    PCR <- any(grepl("c_tot", parTab$names))

    prediction_compartments <- paste0("V", ifelse(camostat_vec, "_Camostat", ""), ifelse(amphoB_vec, "_AmphoB", ""))
    if(PCR) {
      prediction_compartments <- c(prediction_compartments,
                                   gsub("V", "V_tot", prediction_compartments))
    }

    # function to calculate log likelihood
    f <- function(pars) {
      if(PRIOR_FUNC(pars) == -Inf) {
        return(list(lik = -Inf, misc = NA))
      }

      transformed_pars <- transform_pars(pars)

      sigma <- transformed_pars[["sigma"]]
      obs_threshold <- transformed_pars[["obs_threshold"]]
      if(PCR) {
        sigma_tot <- transformed_pars[["sigma_tot"]]
      }

      if(PCR) {
        sigma_vec <- rep(c(sigma, sigma_tot), each = 2 * length(camostat_vec))
      } else {
        sigma_vec <- rep(sigma, 2 * length(camostat_vec))
      }
      obs_threshold_vec <- rep(obs_threshold, length(sigma_vec)) # no PCR data points below threshold so can set to whatever

      solving_time <- make_solving_time(data_df$t)
      if(PCR) {
        transformed_pars$V_tot0 <- transformed_pars$V_0 * transformed_pars[[grep("RNA_ratio_inoculum", names(transformed_pars))]]
      }
      sol <- calc_sol_drugs(transformed_pars, solving_time, camostat_vec, amphoB_vec)

      sol <- sol[,c("t", prediction_compartments)]

      individual_lik <- calc_LL_log10normal_threshold_wrapper(
        data_df,
        prediction_compartments,
        sol,
        prediction_compartments,
        sigma_vec,
        obs_threshold_vec)

      misc <- sol %>%
        as.list %>%
        lapply(name_predictions) %>%
        unlist

      lik <- sum(individual_lik)

      list(lik = lik, misc = misc)

    }
    f
  }

# function to name a vector of predictions
name_predictions <- function(prediction) {
  names(prediction) <- as.character(seq_along(prediction))
  prediction
}
