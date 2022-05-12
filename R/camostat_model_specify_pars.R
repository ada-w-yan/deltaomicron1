#' make parameter table
#'
#' @param n_L: integer. number of latent stages
#' @param n_I: integer.  number of infectious stages
#' @param assay: string.  "plaque" or "tcid50" (only affects observation threshold)
#' @param reparameterise: logical.  if FALSE, use beta, p etc, if TRUE, use probability of infection (beta T_0) / c, burst size, delta, c
#' @param k1_endosomal_limit logical. if TRUE, 0 <= k1_endosomal/k1_tmprss2 <= 1
#' @param strain_names vector of strings.  Names of virus strains.
#' @return data_df frame with columns
#' values: default parameter value
#' names: internal name for parameter
#' fixed:  if 0, fit parameter, otherwise fix parameter
#' lower_bound: lower bound of parameter for fitting
#' upper_bound: upper bound of parameter for fitting
#' steps: width of proposal distribution in parameter space
#' names_plot: character vector which which to label parameter distribution plot
specify_camostat_parameters_fn <-
  function(n_L = 1,
           n_I = 1,
           different_eclipse = FALSE,
           assay = "plaque",
           PCR = FALSE,
           reparameterise = FALSE,
           k1_endosomal_limit = FALSE,
           strain_names) {

    if((!different_eclipse) && k1_endosomal_limit) {
      stop("cannot have k1_endosomal_limit = TRUE and different_eclipse = FALSE")
    }
    n_strains <- length(strain_names)
    if(n_strains > 1 && (!k1_endosomal_limit || !reparameterise || !different_eclipse)) {
      stop("not yet implemented")
    }
    if(PCR && any(strain_names %in% c("WT", "DeltaCS"))) {
      stop("no PCR data for WT and DeltaCS")
    }

    default_step <- 0.1

    # inefficient but easier to read
    # check lower and upper bounds for everything

    # initial number of target cells
    T_0 <- 5e5

    # supernatant volume differs between experiments
    supernatant_vol <- .2

    parTab <- data.frame(
      values = -0.5,
      names = "log10_c_inf",
      fixed = 0,
      lower_bound = -2,
      upper_bound = -.5,
      steps = default_step,
      names_plot = "$\\log_{10}c_{inf}$",
      stringsAsFactors = FALSE
    )

    if(PCR) {
      parTab <- rbind(parTab, data.frame(
        values = 0.5,
        names = "c_tot",
        fixed = 0,
        lower_bound = 0,
        upper_bound = 10^-.5,
        steps = default_step,
        names_plot = "$c_{tot}$",
        stringsAsFactors = FALSE
      ))
    }

    parTab <-
      rbind(
        parTab,
        data.frame(
          values = supernatant_vol,
          names = "supernatant_vol",
          fixed = 1,
          lower_bound = -Inf,
          upper_bound = Inf,
          steps = default_step,
          names_plot = "supernatant_vol"
        )
      )

    if(reparameterise) {
      pathways <- c("tmprss2", "endosomal")
      if(n_strains > 1) {
        make_parTab_strain <- function(strain_name) {
          data.frame(
            values = .4,
            names = paste0("prob_infection_", pathways, "_", strain_name),
            fixed = 0,
            lower_bound = 0,
            upper_bound = 1,
            steps = default_step,
            names_plot = paste0("prob infection ", pathways, " ", strain_name)
          )
        }
        parTab_strains <- lapply(strain_names, make_parTab_strain) %>%
          bind_rows
        parTab <- rbind(parTab, parTab_strains)
      } else {
        parTab <-
          rbind(
            parTab,
            data.frame(
              values = .4,
              names = paste0("prob_infection_", pathways),
              fixed = 0,
              lower_bound = 0,
              upper_bound = 1,
              steps = default_step,
              names_plot = paste0("prob infection ", pathways)
            )
          )
      }

    } else {
      parTab <-
        rbind(
          parTab,
          data.frame(
            values = -5,
            names = "log10_beta_tmprss2",
            fixed = 0,
            lower_bound = -10,
            upper_bound = -1,
            steps = default_step,
            names_plot = "$\\log_{10}\\beta_{TMPRSS2}$"
          )
        )

      parTab <-
        rbind(
          parTab,
          data.frame(
            values = -5,
            names = "log10_beta_endosomal",
            fixed = 0,
            lower_bound = -10,
            upper_bound = -1,
            steps = default_step,
            names_plot = "$\\log_{10}\\beta_endosomal$"
          )
        )
    }


    if(different_eclipse) {
      parTab <-
        rbind(
          parTab,
          data.frame(
            values = .25,
            names = "k1_endosomal",
            fixed = 0,
            lower_bound = 1/72,
            upper_bound = 1/2,
            steps = default_step,
            names_plot = "$k_endosomal$"
          )
        )

      if(k1_endosomal_limit) {
        parTab <-
          rbind(
            parTab,
            data.frame(
              values = .25,
              names = "k1_endosomal_on_k1_tmprss2",
              fixed = 0,
              lower_bound = 0,
              upper_bound = 1,
              steps = default_step,
              names_plot = "$k_endosomal/k_tmprss2$"
            )
          )
      } else {
        parTab <-
          rbind(
            parTab,
            data.frame(
              values = .25,
              names = "k1_tmprss2",
              fixed = 0,
              lower_bound = 1/72,
              upper_bound = 1/2,
              steps = default_step,
              names_plot = "$k_tmprss2$"
            )
          )
      }

    } else {
      parTab <-
        rbind(
          parTab,
          data.frame(
            values = .25,
            names = "k1",
            fixed = 0,
            lower_bound = 1/24,
            upper_bound = 1/2,
            steps = default_step,
            names_plot = "$k$"
          )
        )
    }


    parTab <-
      rbind(
        parTab,
        data.frame(
          values = -.5,
          names = "log10_delta",
          fixed = 0,
          lower_bound = -2.5,
          upper_bound = 0,
          steps = default_step,
          names_plot = "$\\log_{10}\\delta$"
        )
      )

    # note the priors are not equivalent
    if(reparameterise) {

      parTab <- rbind(parTab,
                      data.frame(values = 0,
                                 names = "log10_burst_size",
                                 fixed = 0,
                                 lower_bound = 0,
                                 upper_bound = 6,
                                 steps = default_step,
                                 names_plot = "log 10 burst size"))
    } else {
      parTab <- rbind(parTab,
                      data.frame(values = 0,
                                 names = "log10_p_inf",
                                 fixed = 0,
                                 lower_bound = -2,
                                 upper_bound = 3,
                                 steps = default_step,
                                 names_plot = "$\\log_{10}p_{inf}$"))
    }

    if(n_L > 1) {
      parTab <- rbind(
        parTab,
        data.frame(
          values = n_L,
          names = "n_L",
          fixed = 1,
          lower_bound = -Inf,
          upper_bound = Inf,
          steps = default_step,
          names_plot = "n_L"
        )
      )
    }

    if(n_I > 1) {
      parTab <- rbind(
        parTab,
        data.frame(
          values = n_I,
          names = "n_I",
          fixed = 1,
          lower_bound = -Inf,
          upper_bound = Inf,
          steps = default_step,
          names_plot = "n_I"
        )
      )
    }

    ## initial conditions

    parTab <- rbind(
      parTab,
      data.frame(
        values = log10(0.05),
        names = "log10_MOI",
        fixed = 0,
        lower_bound = log10(0.05) - .5,
        upper_bound = log10(0.05) + .5,
        steps = default_step,
        names_plot = "log10 MOI"
      )
    )

    # hardcoded from inoculum data, not great
    if(PCR) {
      parTab_RNA_ratio_inoculum <- data.frame(values = c(255000, 222000),
                                              names = paste0("RNA_ratio_inoculum_", c("omicron", "delta")),
                                              fixed = 1,
                                              lower_bound = -Inf,
                                              upper_bound = Inf,
                                              steps = default_step,
                                              names_plot = "RNA_ratio_inoculum")
      parTab_RNA_ratio_inoculum <- parTab_RNA_ratio_inoculum[parTab_RNA_ratio_inoculum$names %in% paste0("RNA_ratio_inoculum_", strain_names),]
      parTab <- rbind(parTab, parTab_RNA_ratio_inoculum)
    }

    parTab <- rbind(
      parTab,
      data.frame(
        values = 0.5,
        names = "prop_wash",
        fixed = 0,
        lower_bound = 0,
        upper_bound = .1,
        steps = default_step,
        names_plot = "prop wash"
      )
    )

    if(PCR) {
      parTab <- rbind(
        parTab,
        data.frame(
          values = 0.5,
          names = "log10_p_tot_on_p_inf",
          fixed = 0,
          lower_bound = 0,
          upper_bound = 6,
          steps = default_step,
          names_plot = "log10_p_tot_on_p_inf"
        )
      )
    }


    parTab <-
      rbind(
        parTab,
        data.frame(
          values = 1,
          names = "incubation_period",
          fixed = 1,
          lower_bound = -Inf,
          upper_bound = Inf,
          steps = default_step,
          names_plot = "incubation period"
        )
      )

    parTab <-
      rbind(
        parTab,
        data.frame(
          values = T_0,
          names = "T_0",
          fixed = 1,
          lower_bound = -Inf,
          upper_bound = Inf,
          steps = default_step,
          names_plot = "$T_0$"
        )
      )

    parTab <-
      rbind(
        parTab,
        data.frame(
          # Table 1, weight all studies equally
          values = (0.0376 + 0.0158 + 0.0043 + 0.0471) / 4,
          names = "prop_ace2_pos",
          fixed = 1,
          lower_bound = -Inf,
          upper_bound = Inf,
          steps = default_step,
          names_plot = "$ACE2^+$"
        )
      )

    parTab <-
      rbind(
        parTab,
        data.frame(
          # Table 1, weight all studies equally
          values = (1.09 / 3.76 + 0.73 / 1.58 + 0.06 / 0.43 + 2.52 / 4.71) / 4,
          # proportion tmprss2+ out of ace2+
          names = "prop_tmprss2_pos",
          fixed = 1,
          lower_bound = -Inf,
          upper_bound = Inf,
          steps = default_step,
          names_plot = "$TMPRSS2^+$"
        )
      )
    parTab <-
      rbind(
        parTab,
        data.frame(
          values = .4,
          names = "sigma",
          fixed = 0,
          lower_bound = 0,
          upper_bound = 1,
          steps = default_step,
          names_plot = "$\\sigma$"
        )
      )

    if(PCR) {
      parTab <-
        rbind(
          parTab,
          data.frame(
            values = .4,
            names = "sigma_tot",
            fixed = 0,
            lower_bound = 0,
            upper_bound = 1,
            steps = default_step,
            names_plot = "$\\sigma_{tot}$"
          )
        )
    }

    obs_threshold <- switch(assay,
                            plaque = 10,
                            tcid50 = 3,
                            stop("unknown assay"))
    obs_threshold <- obs_threshold * supernatant_vol

    parTab <-
      rbind(
        parTab,
        data.frame(
          values = obs_threshold,
          names = "obs_threshold",
          fixed = 1,
          lower_bound = -Inf,
          upper_bound = Inf,
          steps = default_step,
          names_plot = "obs threshold"
        )
      )


    # observation threshold of 10 pfu.  treat values below threshold as censored

    # convert logicals to numerics
    parTab$fixed <- as.numeric(parTab$fixed)
    # set bounds for fixed parameters to c(-Inf Inf)
    parTab[parTab$fixed == 1, "lower_bound"] <- -Inf
    parTab[parTab$fixed == 1, "upper_bound"] <- Inf

    parTab
  }

#' make parameter table where only some parameters are fitted
#' specify fixed parameter values
#' @param old_values named numeric vector of fixed parameter values
#' @param old_parTab old parameter table
#' @param subset_pars character vector.  Subset of parameters to be fitted.
#' @return data frame: new parameter table
specify_camostat_parameters_fn_fit_subset_pars <- function(old_values, old_parTab, subset_pars) {
  # old_max_LL_values <- get_max_LL_params(paste0(old_dir_name, "1.RData"))
  parTab <- old_parTab
  parTab$fixed <- 1
  parTab[parTab$names %in% subset_pars, "fixed"] <- 0
  parTab[parTab$fixed == 1, "lower_bound"] <- -Inf
  parTab[parTab$fixed == 1, "upper_bound"] <- Inf
  parTab$values <- old_values
  parTab
}

#'
#'
specify_camostat_parameters_fn_share_subset_pars <- function(n_L = 1,
                                                             n_I = 1,
                                                             different_eclipse = FALSE,
                                                             assay = "plaque",
                                                             PCR = FALSE,
                                                             reparameterise = FALSE,
                                                             k1_endosomal_limit = FALSE,
                                                             strain_names,
                                                             separate_pars) {
  # to do: implement other parameters

  stopifnot(all(separate_pars %in% c("prob_infection_endosomal", "prob_infection_tmprss2", "log10_burst_size")))
  parTab <- specify_camostat_parameters_fn (n_L = n_L,
             n_I = n_I,
             different_eclipse = different_eclipse,
             assay = assay,
             PCR = PCR,
             reparameterise = reparameterise,
             k1_endosomal_limit = k1_endosomal_limit,
             strain_names = strain_names)
  # temp fix
  parTab <- parTab[!grepl("prob_infection", parTab$names),]
  pathways <- c("tmprss2", "endosomal")
  parTab <- rbind(
    parTab,
    data.frame(
      values = .4,
      names = paste0("prob_infection_", pathways),
      fixed = 0,
      lower_bound = 0,
      upper_bound = 1,
      steps = .1,
      names_plot = paste0("prob infection ", pathways)
    )
  )

  parTab_separate_pars <- parTab[parTab$names %in% separate_pars,]
  parTab_together_pars <- parTab[!(parTab$names %in% separate_pars),]
  make_strain_parTab <- function(strain_name) {
    parTab_separate_pars$names <- paste0(parTab_separate_pars$names, "_", strain_name)
    parTab_separate_pars
  }
  parTab_separate <- lapply(strain_names, make_strain_parTab) %>%
    bind_rows()
  parTab <- bind_rows(parTab_together_pars, parTab_separate)
  parTab
}

specify_camostat_amphoB_parameters_fn <- function(strain, Calu3, PCR) {
  parTab <- specify_camostat_parameters_fn(n_L = 10,
             n_I = 10,
             different_eclipse = TRUE,
             assay = "plaque",
             PCR = PCR,
             reparameterise = TRUE,
             k1_endosomal_limit = TRUE,
             strain_names = strain)
  # susceptibility to IFITM -- 0 is 0 susceptibility, 1 is complete abrogation of infection
  parTab <- rbind(parTab, data.frame(
    values = 0,
    names = "IFITM",
    fixed = 1,
    lower_bound = 0,
    upper_bound = 1,
    steps = .1,
    names_plot = "IFITM",
    stringsAsFactors = FALSE
  ))

  if(Calu3) {
    parTab[parTab$names %in% c("prop_ace2_pos", "prop_tmprss2_pos"), "values"] <- 1
    parTab[parTab$names == "obs_threshold", "values"] <- 25
    parTab[parTab$names == "supernatant_vol", "values"] <- .5
    parTab[parTab$names == "log10_MOI", c("values", "lower_bound", "upper_bound")] <-
      log10(0.01) + c(0, -.5, .5)
  }
  parTab
}
