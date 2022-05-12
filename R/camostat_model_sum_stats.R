calc_beta_pathway <- function(pars, pathway) {
  if("IFITM" %in% names(pars) && !grepl("no_IFITM", pathway)) {
    pars$beta_endosomal <- pars$beta_endosomal * (1 - pars$IFITM)
  }
  if(pathway == "tmprss2") {
    beta1 <- c(pars$beta_tmprss2, 0)
  } else if(grepl("endosomal", pathway)) {
    beta1 <- rep(pars$beta_endosomal, 2)
  } else {
    beta1 <- c(pars$beta_tmprss2, 0) + rep(pars$beta_endosomal, 2)
  }
  beta1
}

get_pathways <- function(pars) {
  if("IFITM" %in% names(pars)) {
    pathways <- c("tmprss2", "endosomal", "both", "endosomal_no_IFITM", "both_no_IFITM")
  } else {
    pathways <- c("tmprss2", "endosomal", "both")
  }
  pathways
}

calc_R_0_camostat <- function(pars) {

  calc_R_0_pathway <- function(pathway) {
    pars$beta <- calc_beta_pathway(pars, pathway)
    R_0 <- sum(pars$beta * pars$p_inf * pars$T_0 / pars$delta) /
      (pars$c_inf + sum(pars$beta * pars$T_0))

    R_0
  }

  pathways <- get_pathways(pars)
  R_0 <- vnapply(pathways, calc_R_0_pathway)

  names(R_0) <- paste0("R_0_", names(R_0))
  R_0
}

calc_r_camostat <- function(pars, warn_if_neg = FALSE) {
  calc_r_pathway <- function(pathway) {
    pars$beta <- calc_beta_pathway(pars, pathway)

    eigenvalue_mat <- matrix(c(-pars$k1, 0, pars$k1, 0, 0,
                               0, -pars$k1, 0, pars$k1, 0,
                               0, 0, -pars$delta, 0, pars$p_inf,
                               0, 0, 0, -pars$delta, pars$p_inf,
                               pars$beta[1] * pars$T_0[1], pars$beta[2] * pars$T_0[2], 0, 0,
                               -(pars$beta[1] * pars$T_0[1] + pars$beta[2] * pars$T_0[2] + pars$c_inf)),
                             ncol = 5)
    ev <- eigen(eigenvalue_mat)
    ev <- ev$values
    ev <- ev[Im(ev) == 0]
    stopifnot(length(ev) > 0)
    r <- max(as.double(ev))
    if(warn_if_neg && r <= 0) {
      warning("negative value of r")
    }
    r
  }

  pathways <- get_pathways(pars)
  r <- vnapply(pathways, calc_r_pathway)

  names(r) <- paste0("r_", names(r))
  r
}

calc_gen_time_camostat <- function(pars) {

  calc_gen_time_pathway <- function(pathway) {
    if(!("n_I" %in% names(pars))) {
      pars$n_I <- 1
    }
    pars$beta <- calc_beta_pathway(pars, pathway)
    virus_mean_time <- 1 / (pars$c_inf + sum(pars$beta * pars$T_0))
    mean_latent_period <- 1/pars$k1
    mean_virion_production_time <- 1 / pars$delta * (pars$n_I + 1) / 2 / pars$n_I

    virus_mean_time + mean_latent_period + mean_virion_production_time
  }

  pathways <- get_pathways(pars)
  gen_time <- vnapply(pathways, calc_gen_time_pathway)
  names(gen_time) <- paste0("gen_time_", names(gen_time))
  gen_time

}

#' create closure to calculate summary statistics and provide bounds for them
#'
#' @param parTab data frame
#' @return a list with elements
#' "calc_summary": closure to calculate
gen_summary_statistics_fn_camostat <-
  function(parTab) {
    different_eclipse <- "k1_endosomal" %in% parTab$names
    reparameterise <- "log10_burst_size" %in% parTab$names

    strain_names <- get_strain_names(parTab$names)
    n_strains <- length(strain_names)

    dummy_pars <- parTab$names
    names(dummy_pars) <- parTab$names
    pathways <- get_pathways(dummy_pars)

    if(different_eclipse) {
      sum_stat_names <- c("R_0", "r", "doubling_time")
    } else {
      sum_stat_names <- c("R_0", "r", "gen_time", "doubling_time")
    }

    sum_stat_names_pathway <- outer(sum_stat_names, pathways, paste0) %>% t %>%
      as.character
    # sum_stat_names_pathway <- c(sum_stat_names_pathway, "R_0_tmprss2_on_R_0_both", "R_0_endosomal_on_R_0_both")
    # if(reparameterise) {
      # sum_stat_names_reparameterise <- c("log10_beta_endosomal", "log10_beta_tmprss2", "log10_p_inf", "prob_infect")
      # sum_stat_names_pathway <- c(sum_stat_names_pathway, sum_stat_names_reparameterise)
    # }

    stages <- "n_L" %in% parTab$names

    # function to transform parameters to calculate summary statistics

    transform_pars <- transform_pars_wrapper_camostat(parTab)

    # functions to calculate summary statistics given parameter values

    # closure to generate function to calculate summary statistics given parameter values
    calc_summary <- function(values) {

      values <- transform_pars(values)
      if(different_eclipse) {
        calc_funcs <- list(calc_R_0_camostat, calc_r_camostat_stages_different_eclipse)
      } else {
        calc_funcs <- list(calc_R_0_camostat, calc_r_camostat, calc_gen_time_camostat)

        if(stages) {
          calc_funcs[[2]] <- calc_r_camostat_stages
        }
      }

      calc_summary_strain <- function(strain) {
        names(values) <- gsub(paste0("_", strain), "", names(values))
        sum_stats <-
          vapply(calc_funcs, function(f)
            f(values), double(length(pathways))) %>%
          as.numeric

        names(sum_stats) <- sum_stat_names_pathway[seq_along(sum_stats)]

        doubling_times <- log(2) / sum_stats[paste0("r", pathways)]
        names(doubling_times) <- paste0("doubling_time", pathways)
        sum_stats <- c(sum_stats, doubling_times)

        # sum_stats[["R_0_tmprss2_on_R_0_both"]] <- sum_stats[["R_0tmprss2"]] / sum_stats[["R_0both"]]
        # sum_stats[["R_0_endosomal_on_R_0_both"]] <- sum_stats[["R_0endosomal"]] / sum_stats[["R_0both"]]

        # if(reparameterise) {
        #   sum_stat_names_reparameterise <- gsub("log10_", "", sum_stat_names_reparameterise)
        #   sum_beta_T_0 <- sum(calc_beta_pathway(values, "both") * values[["T_0"]])
        #   values[["prob_infect"]] <- 10^(sum_beta_T_0 / (sum_beta_T_0 + values[["c_inf"]]))
        #   sum_stats <- c(sum_stats, log10(unlist(values[sum_stat_names_reparameterise])))
        # }

        if(n_strains > 1) {
          names(sum_stats) <- paste0(names(sum_stats), "_", strain)
        }
        sum_stats
      }

      sum_stats <- lapply(strain_names, calc_summary_strain) %>% unlist
      sum_stats
    }

    dummy_values <- calc_summary(parTab$values)
    par_names_plot <- names(dummy_values)
    lower_bound <- rep(0, length(par_names_plot))
    if(different_eclipse) {
      upper_bound <- rep(c(100, 1, 100), length(pathways))#, 1, 1)
    } else {
      upper_bound <- c(rep(c(100, 1, 100, 100), 3), 1, 1)
    }
    # if(reparameterise) {
    #   upper_bound <- c(upper_bound, 1, 1, 1e3, 1)
    # }

    upper_bound <- rep(upper_bound, n_strains)
    names(par_names_plot) <- names(lower_bound) <- names(upper_bound) <- par_names_plot

    # nice parameter names for plotting

    list(
      calc_summary = calc_summary,
      lower_bound = lower_bound,
      upper_bound = upper_bound,
      par_names_plot = par_names_plot
    )
  }


#' calculates r for the TEIV model with Erlang-distributed latent and infectious periods, and different cell types
#'
#' @param params a list or named numeric vector with the parameter values for
#' beta: infectivity in inverse units of V (in dT/dt)
#' T_0: initial number of target cells in units of T
#' p: production rate of infectious virus from cells
#' c_inf: decay rate of infectious virus
#' k1: rate from latent to infectious
#' delta: decay rate of infectious cells
#' n_L: shape parameter for distribution of latent period
#' n_I: shape parameter for distribution of infectious period
#' beta_inf: infectivity in inverse units of T (in dV/dt)
#' if loss of virus due to enry into target cells is ignored, beta_inf = 0
#' @return r: numeric vector of length 1
calc_r_camostat_stages <- function(params, warn_if_neg = FALSE) {
  calc_r_pathway <- function(pathway) {
    params$beta <- calc_beta_pathway(params, pathway)
    # different parameter naming systems
    beta_T_0 <- params[["beta"]] * params[["T_0"]]
    p <- params[["p_inf"]]
    c_inf <- params[["c_inf"]] + sum(params[["beta"]] * params[["T_0"]])
    k <- params[["k1"]]
    delta <- params[["delta"]]
    n_L <- params[["n_L"]]
    n_I <- params[["n_I"]]
    n_C <- length(params[["T_0"]])



    make_LL_II_mat <- function(n_row, value) {

      stopifnot(n_row > 0, isTRUE(all.equal(round(n_row), n_row)))

      if(n_row == 1) {
        matrix(-n_row * value)
      } else {
        vec <- double(n_row)
        vec[seq_len(2)] <- c(-n_row * value, n_row * value)
        mat <- toeplitz(vec)
        mat[upper.tri(mat)] <- 0
        mat
      }
    }

    make_mats <- function(k, delta) {
      LL_mat <- make_LL_II_mat(n_L, k)

      LI_mat <- matrix(0, nrow = n_L, ncol = n_I)

      IL_mat <- t(LI_mat)
      IL_mat[1,n_L] <- n_L * k

      II_mat <- make_LL_II_mat(n_I, delta)
      mat1 <- rbind(LL_mat, IL_mat)
      mat2 <- rbind(LI_mat, II_mat)
      cbind(mat1, mat2)
    }

    LI_mats <- Map(make_mats, k, delta)
    LI_mats <- Matrix::bdiag(LI_mats)

    # for more than one cell type
    LI_mats <- rep(list(LI_mats), n_C)
    LI_mats <- Matrix::bdiag(LI_mats)
    #
    p_vec <- vapply(p, function(x) c(double(n_L), rep(x, n_I)), numeric(n_L + n_I)) %>%
      as.numeric
    p_vec <- rep(p_vec, n_C)
    V_mat_row <- matrix(c(p_vec, -c_inf), nrow = 1)
    V_mat_col <- vapply(beta_T_0, function(x) c(x, double(n_L + n_I - 1)), numeric(n_L + n_I)) %>%
      as.numeric %>%
      matrix(., ncol = 1)
    mat3 <- cbind(LI_mats, V_mat_col)

    eigenvalue_mat <- rbind(mat3, V_mat_row)
    ev <- eigen(eigenvalue_mat)
    ev <- ev$values
    ev <- ev[Im(ev) == 0]
    stopifnot(length(ev) > 0)
    r <- max(as.double(ev))
    if(warn_if_neg && r <= 0) {
      warning("negative value of r")
    }
    r
  }
  pathways <- get_pathways(params)

  r <- vnapply(pathways, calc_r_pathway)
  names(r) <- paste0("r_", names(r))
  r
}

#' calculates r for the TEIV model with Erlang-distributed latent and infectious periods,
#' with TMPRSS2+ and TMPRSS2- cells, where the eclipse phase for the
#' TMPRSS2 and endosomal pathways are different
#'
#' @param params a list or named numeric vector with the parameter values for
#' beta: infectivity in inverse units of V (in dT/dt)
#' T_0: initial number of target cells in units of T
#' p: production rate of infectious virus from cells
#' c_inf: decay rate of infectious virus
#' k1: rate from latent to infectious
#' delta: decay rate of infectious cells
#' n_L: shape parameter for distribution of latent period
#' n_I: shape parameter for distribution of infectious period
#' beta_inf: infectivity in inverse units of T (in dV/dt)
#' if loss of virus due to enry into target cells is ignored, beta_inf = 0
#' @return r: numeric vector of length 1
calc_r_camostat_stages_different_eclipse <- function(params, warn_if_neg = FALSE) {

  # if endosomal pathway only, set k1 to k1_endosomal and use calc_r_camostat_stages
  # if tmprss2 pathway only, set k1 to k1_tmprss2 and use calc_r_camostat_stages

  params_endosomal <- params_tmprss2 <- params
  params_endosomal[["k1"]] <- params[["k1_endosomal"]]
  params_tmprss2[["k1"]] <- params[["k1_tmprss2"]]

  r_endosomal <- calc_r_camostat_stages(params_endosomal, warn_if_neg = warn_if_neg)
  r_endosomal_no_IFITM <- r_endosomal[["r_endosomal_no_IFITM"]]
  r_endosomal <- r_endosomal[["r_endosomal"]]

  r_tmprss2 <- calc_r_camostat_stages(params_tmprss2, warn_if_neg = warn_if_neg)
  r_tmprss2 <- r_tmprss2[["r_tmprss2"]]

  # now calc r for both pathways together

  calc_r_inner <- function(pathway) {
    params$beta <- calc_beta_pathway(params, pathway)
    p <- params[["p_inf"]]
    c_inf <- params[["c_inf"]] + sum(params[["beta"]] * params[["T_0"]])
    delta <- params[["delta"]]
    n_L <- params[["n_L"]]
    n_I <- params[["n_I"]]
    n_C <- length(params[["T_0"]])



    make_LL_II_mat <- function(n_row, value) {

      stopifnot(n_row > 0, isTRUE(all.equal(round(n_row), n_row)))

      if(n_row == 1) {
        matrix(-n_row * value)
      } else {
        vec <- double(n_row)
        vec[seq_len(2)] <- c(-n_row * value, n_row * value)
        mat <- toeplitz(vec)
        mat[upper.tri(mat)] <- 0
        mat
      }
    }

    LL_endosomal_mat <- make_LL_II_mat(n_L, params[["k1_endosomal"]])
    LL_tmprss2_mat <- make_LL_II_mat(n_L, params[["k1_tmprss2"]])
    II_mat <- make_LL_II_mat(n_I, delta)

    LI_mats <- Matrix::bdiag(list(LL_endosomal_mat, LL_tmprss2_mat, LL_endosomal_mat, II_mat, II_mat))
    LI_mats[3*n_L + 1, n_L] <- LI_mats[4 * n_L + 1, 3 * n_L] <- n_L* params[["k1_endosomal"]]
    LI_mats[3*n_L + 1, 2 * n_L] <- n_L* params[["k1_tmprss2"]]

    p_vec <- vapply(p, function(x) c(double(n_L), rep(x, n_I)), numeric(n_L + n_I)) %>%
      as.numeric
    p_vec <- rep(p_vec, n_C)
    V_mat_row <- matrix(c(double(3 * n_L), rep(p, 2 * n_I), -c_inf), nrow = 1)
    V_mat_col <- matrix(double(3 * n_L + 2 * n_I), ncol = 1)
    if(pathway == "both") {
      params$beta_endosomal <- params$beta_endosomal * (1 - params$IFITM)
    }

    V_mat_col[1,1] <- params[["beta_endosomal"]] * params[["T_0"]][1]
    V_mat_col[n_L + 1,1] <- params[["beta_tmprss2"]] * params[["T_0"]][1]
    V_mat_col[2 * n_L + 1,1] <- params[["beta_endosomal"]] * params[["T_0"]][2]
    mat3 <- cbind(LI_mats, V_mat_col)

    eigenvalue_mat <- rbind(mat3, V_mat_row)
    ev <- eigen(eigenvalue_mat)
    ev <- ev$values
    ev <- ev[Im(ev) == 0]
    stopifnot(length(ev) > 0)
    r_both <- max(as.double(ev))
    if(warn_if_neg && r_both <= 0) {
      warning("negative value of r")
    }
    r_both
  }

  pathways <- c("both", "both_no_IFITM")
  r_both <- vnapply(pathways, calc_r_inner)
  names(r_both) <- pathways
  r_both_no_IFITM <- r_both[["both_no_IFITM"]]
  r_both <- r_both[["both"]]

  r <- c(r_tmprss2 = r_tmprss2, r_endosomal = r_endosomal, r_both = r_both,
         r_endosomal_no_IFITM = r_endosomal_no_IFITM,
         r_both_no_IFITM = r_both_no_IFITM)
  r
}

#' test whether the ratio of a summary statistic between two strains
#' differs from 1
#'
#' @param data_dir character vector of length 2.  directories where results for the two strains live
#' @param par_name parameter to test
#' @return a list of length n, where n is the number of pair combinations of data_dir.
#' each element contains a list with the elements prctile_ratios and p_value.
#' prctile_ratio is a vector which gives the 2.5, 50 and 97.5th percentiles of the ratio of each
#' p_value is a scalar where each element is twice the proportion of the
#' sampled ratios greater than 1, or the proportion of the sampled ratios smaller
#' than 1 (whichever is the smaller)
#' @export
pairwise_stats_tests <- function(data_dir, par_name) {

  # read in parameters
  chains <- lapply(paste0(data_dir, "1_summary_chain_combined.csv"), readr::read_csv, col_select = par_name)
  chains <- lapply(chains, unlist)

  # undo log transforms of parameters
  if(grepl("log10", par_name)) {
    chains <- lapply(chains, function(x) 10^x)
  }

  n_samples <- 1e6

  # take n_samples samples of each summary statistic for each strain
  samples <- vapply(chains, sample, numeric(n_samples), size = n_samples, replace = TRUE)

  # calculate the ratios of the sampled summary statistics between the strains
  ratio_samples <- samples[,1]/ samples[,2]

  # calculate the percentiles and p values of the ratios
  prctile <- quantile(ratio_samples, c(.025, .50, .975))
  p_value <- min(sum(ratio_samples > 1), sum(ratio_samples < 1)) / length(ratio_samples) * 2

  list(prctile = prctile, p_value = p_value)
}
