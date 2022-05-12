#' determine strain names from parameter names
#'
#' @param par_names character vector: parameter names
#' @return character vector string names
get_strain_names <- function(par_names) {
  par_names <- gsub("log10_", "", par_names)
  reparameterise <- any(grepl("burst_size", par_names))
  # note only works if prob_infection differs between strains!
  grep_str <- ifelse(reparameterise, "prob_infection_tmprss2", "beta_tmprss2")
  strain_strs <- grep(grep_str, par_names, value = TRUE)
  strain_names <- gsub(paste0(grep_str, "_"), "", strain_strs)
  strain_names
}


#' create a closure which transforms an unnamed numeric vector into
#' lists of named parameters to calculate the likelihood
#'
#' @param parTab data frame containing priors
#' @return a function which transforms an unnamed numeric vector into
#' lists of named parameters to calculate the likelihood
transform_pars_wrapper_camostat <- function(parTab) {

  par_names <- parTab$names
  log10_ind <- grep("log10", par_names)
  par_names <- sub("log10_", "", par_names)

  transform_pars <- function(pars) {
    pars[log10_ind] <- 10 ^(pars[log10_ind])
    names(pars) <- par_names
    n_C <- 2 # number of cell types
    pars <- as.list(pars)

    different_eclipse <- any(grepl("k1_endosomal",  par_names))
    reparameterise <- any(grepl("burst_size", par_names))
    k1_endosomal_limit <- any(grepl("k1_endosomal_on_k1_tmprss2", par_names))

    strain_names <- get_strain_names(par_names)
    n_strains <- length(strain_names)

    calc_beta_from_prob <- function(pars, pathway, strain_name = NA) {

      prob_name <- paste0("prob_infection_", pathway)
      if(!is.na(strain_name)) {
        prob_name <- paste0(prob_name, "_", strain_name)
      }

      T_0 <- pars$T_0 * pars$prop_ace2_pos
      if(pathway == "tmprss2") {
        T_0 <- T_0 * pars$prop_tmprss2_pos
      }
      pars$c_inf * pars[[prob_name]] / (1 - pars[[prob_name]]) / T_0
    }

    pathways <- c("endosomal", "tmprss2")
    if(reparameterise) {
      # only works if beta is different for each strain
      if(n_strains > 1) {
        beta_grid <- expand_grid(pathway = pathways, strain_name = strain_names)
        betas <- Map(function(x, y) calc_beta_from_prob(pars, x, y), beta_grid$pathway, beta_grid$strain_name)

        names(betas) <- paste0("beta_", beta_grid$pathway, "_", beta_grid$strain_name)
      } else {
        betas <- lapply(pathways, calc_beta_from_prob, pars = pars)
        names(betas) <- paste0("beta_", pathways)
      }
      pars <- c(pars, betas)
      # TBC
      if(n_strains > 1 && paste0("burst_size_", strain_names[1]) %in% names(pars)) {
        make_pars_list <- function(strain_name) {
          p_inf <- pars[[paste0("burst_size_", strain_name)]] * pars$delta
          p_tot <- pars$p_tot_on_p_inf * p_inf
          values <- list(p_inf = p_inf, p_tot = p_tot)
          names(values) <- paste0(c("p_inf_", "p_tot_"), strain_name)
          values
        }
        values <- lapply(strain_names, make_pars_list) %>% unlist
        pars <- c(pars, values)
      } else {
        pars$p_inf <- pars$burst_size * pars$delta
        pars$p_tot <- pars$p_tot_on_p_inf * pars$p_inf
      }

    }

    if(k1_endosomal_limit) {
      pars <- c(pars, list(k1_tmprss2 = pars$k1_endosomal / pars$k1_endosomal_on_k1_tmprss2))
    }

    if("n_L" %in% par_names) {
      pars <- c(pars,
                list(n_C = n_C, L_0 = matrix(0, n_C, pars$n_L), I_0 = matrix(0, n_C, pars$n_I),
                     V_0 = pars$MOI * pars$T_0))
      if(different_eclipse) {
        pars$L_0 <- matrix(0, n_C + 1, pars$n_L)
      }
    } else {
      pars <- c(pars,
                list(n_C = n_C, L_0 = double(n_C), I_0 = double(n_C),
                     V_0 = pars$MOI * pars$T_0))
    }

    pars$T_0 <- pars$T_0 * pars$prop_ace2_pos * c(pars$prop_tmprss2_pos, 1 - pars$prop_tmprss2_pos)
    # cell type 1 is tmprss2+

    pars
  }
  transform_pars
}





#' plot viral dynamics
#'
#' @param parTab data frame
#' @param data_df
#' @param chain data table of samples from MCMC chain
#' @return list of ggplot objects: plots of data_df, 95% credible intervals and
#' and sampled trajectories for viral loads
#' @import dplyr
#' @export
plot_dynamics_camostat <- function(parTab, data_df, chain){

  # extract the credible intervals for the experiment being plotted
  chain_names <- names(chain)

  strain_names <- get_strain_names(parTab$names)
  n_strains <- length(strain_names)

  PCR <- any(grepl("tot", parTab$names))

  if(n_strains == 1) {
    data_df <- list(data_df)
    names(data_df) <- strain_names
  }
  plot_dynamics_strain <- function(strain, plot_PCR) {

    # retrieve prediction times from data_df frame
    prediction_times <- sort(unique(c(0, data_df[[strain]]$t)))

    camostat_only <- !any(grepl("no_drug", colnames(data_df[[strain]])))
    # prediction_compartments <- c("T1.1", "T1.2", "L1.1", "L1.2", "I1.1", "I1.2", "V") %>%
    #   outer(c("_no_drug", "_Camostat"), paste0) %>% as.character
    # prediction_compartments <- paste0("V_", c("no_drug", "Camostat"))

    plot_prediction_drug <- function(drug) {
      if(n_strains == 1) {
        prediction_names <- paste0("V_", drug, ".", seq_along(prediction_times))
      } else {
        prediction_names <- paste0("V_", strain, "_", drug, ".", seq_along(prediction_times))
      }

      data_name <- paste0("V_", drug)
      if(plot_PCR) {
        prediction_names <- gsub("V_", "V_tot_", prediction_names)
        data_name <- gsub("V_", "V_tot_", data_name)
        y_label <- "RNA copy number"
        y_max <- 1e12
      } else {
        y_label <- "Viral load (pfu)"
        y_max <- 1e10
      }

      ci_samples <- lapply(c(FALSE, TRUE),
                           function(x) gen_prediction_ci(chain, prediction_names,
                                                         prediction_times, x))

      names(ci_samples) <- c("ci", "samples")

      ci_plot <- ggplot(ci_samples$ci, aes(x = t)) +
        geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), fill = "grey") +
        geom_point(data = data_df[[strain]], aes_string(y = data_name)) +
        geom_line(aes(y = `50%`)) +
        theme_bw() +
        scale_y_log10(y_label) +
        coord_cartesian(ylim = c(1, y_max), expand = FALSE) +
        xlab("Time")

      samples_plot <- ggplot(ci_samples$samples, aes(x = t)) +
        geom_line(aes(y = value, color = variable, group = variable)) +
        geom_point(data = data_df[[strain]], aes_string(y = data_name)) +
        theme_bw() +
        scale_y_log10(y_label) +
        coord_cartesian(ylim = c(1, y_max), expand = FALSE) +
        xlab("Time") +
        theme(legend.position = "none")

      out_list <- list(ci_plot, samples_plot)
      if(n_strains == 1) {
        names(out_list) <- paste0("V_", drug, "_", c("ci", "samples"))
      } else {
        names(out_list) <- paste0("V_", strain, "_", drug, "_", c("ci", "samples"))
      }
      if(plot_PCR) {
        names(out_list) <- gsub("V_", "V_tot_", names(out_list))
      }

      out_list
    }

    drugs <- c("no_drug", "Camostat")
    if(camostat_only) {
      plot_prediction_drug("Camostat")
    } else {
      lapply(drugs, plot_prediction_drug) %>%
        do.call(c, .)
    }
  }

  if(PCR) {
    plot_PCR <- c(FALSE, TRUE)
  } else {
    plot_PCR <- FALSE
  }
  strain_PCR_grid <- expand_grid(strain = strain_names, plot_PCR  = plot_PCR)

  g <- apply_named_args(strain_PCR_grid, 1, plot_dynamics_strain) %>%
    do.call(c, .)
  g
}

#' solve preincubation viral dynamics for one iteration of chain
#'
#' @param subchain data table of samples from MCMC chain, should have 1 row only
#' @return data frame
solve_from_posterior_preincub <- function(subchain) {
  parTab <- specify_camostat_parameters_fn()
  if("data.table" %in% class(subchain)) {
    subchain <- unlist(subchain)
  }
  fitted_pars <- parTab[parTab$fixed == 0, "names"]
  parTab[parTab$names %in% fitted_pars, "values"] <- subchain[fitted_pars]
  transform_pars <- transform_pars_wrapper_camostat(parTab)
  pars <- transform_pars(parTab$values)
  solving_time <- seq(-1, 0, by = .01)
  sol <- calc_sol_camostat_preincub(pars, solving_time, with_camostat = TRUE)
  sol
}

#' create a closure which transforms an unnamed numeric vector into
#' lists of named parameters to calculate the likelihood
#'
#' @param parTab data frame containing priors
#' @return a function which transforms an unnamed numeric vector into
#' lists of named parameters to calculate the likelihood
transform_pars_wrapper_tmprss2_only <- function(parTab) {

  par_names <- parTab$names
  log10_ind <- grep("log10", par_names)
  par_names <- sub("log10_", "", par_names)

  transform_pars <- function(pars) {

    pars[log10_ind] <- 10 ^(pars[log10_ind])
    names(pars) <- par_names
    pars <- as.list(pars)
    if("n_L" %in% names(pars)) {
      pars <- c(pars,
                list(L_0 = double(pars$n_L), I_0 = double(pars$n_I),
                     V_0 = pars$MOI * pars$T_0))
    } else {
      pars <- c(pars,
                list(L_0 = 0, I_0 = 0,
                     V_0 = pars$MOI * pars$T_0))
    }

    pars$T_0 <- pars$T_0 * pars$prop_ace2_pos * pars$prop_tmprss2_pos
    # cell type 1 is tmprss2+

    pars
  }
  transform_pars
}

plot_trajectories <- function(strain1, Calu3 = FALSE, trajectory_filename) {

  dir_name <- paste0(strsplit(trajectory_filename, "/", fixed = TRUE)[[1]][1], "/")

  trajectories <- readRDS(trajectory_filename)
  PCR <- "V_tot" %in% colnames(trajectories)

  if(Calu3) {
    data_df <- read_Calu3_data(strain = strain1)
  } else {
    data_df <- read_camostat_data(strain = strain1)
  }

  if(PCR) {
    data_df <- data_df %>%
      pivot_longer(starts_with("V"), names_to = "model", values_to = "V") %>%
      mutate(PCR = grepl("tot", model),
             model = gsub("_tot", "", model),
             model = recode(model, V = "full",
                            V_Camostat = "endosomal_only",
                            V_AmphoB = "full_no_IFITM",
                            V_Camostat_AmphoB = "endosomal_no_IFITM"))
  } else {
    data_df <- data_df %>%
      rename(full = "V_no_drug", endosomal_only = "V_Camostat") %>%
      pivot_longer(c("full", "endosomal_only"), names_to = "model", values_to = "V")
  }

  trajectories <- readRDS(trajectory_filename)

  data_df <- data_df %>%
    filter(t >= 0)
  plot_inner <- function(data_df, trajectories, PCR, facet1 = FALSE) {

    if(PCR) {
      trajectories <- trajectories %>%
        select(t, V_tot, probs, model) %>%
        pivot_wider(names_from = probs, values_from = V_tot) %>%
        filter(t >= 0)
      y_label = "RNA copy number"
      y_max <- 1e12
      data_df <- data_df %>%
        filter(PCR)
    } else {
      trajectories <- trajectories %>%
        select(t, V, probs, model) %>%
        pivot_wider(names_from = probs, values_from = V) %>%
        filter(t >= 0)
      if("PCR" %in% colnames(data_df)) {
        data_df <- data_df %>% filter(PCR == FALSE)
      }
      y_label = "Viral load (pfu)"
      y_max <- 1e8
    }

    levels_vec <- c("full", "full_no_IFITM", "endosomal_only", "endosomal_no_IFITM", "tmprss2_only")

    data_df <- data_df %>%
      mutate(model = factor(model, levels = levels_vec))
    trajectories <- trajectories %>%
      mutate(model = factor(model, levels = levels_vec))
    facet_labels <- c("Both pathways", "Both pathways no IFITM", "Endosomal", "Endosomal no IFITM", "TMPRSS2")
    names(facet_labels) <- levels_vec
    g <- ggplot(data_df) +
      geom_ribbon(data = trajectories, aes(x = t, ymin = lower, ymax = upper, fill = model, group = model),
                  alpha = .3) +
      geom_point(aes(x = t, y = V,
                     color = model, group = model)) +
      geom_line(data = trajectories, aes(x = t, y = max_LL, color = model, group = model)) +
      scale_y_log10(y_label) +
      coord_cartesian(ylim = c(1, y_max)) +
      theme_bw() +
      xlab("Time (hpi)") +
      scale_color_discrete(drop = FALSE,
                           labels = c("No drug", "Amphotericin B", "Camostat", "Camostat + \n Amphotericin B", ""),
                           name = "Drug") +
      scale_fill_discrete(guide = "none")
    if(facet1) {
      g <- g +
        facet_wrap(~model, nrow = 3, labeller = labeller(model = facet_labels))
    }
    if(!PCR) {
      g <- g + geom_hline(yintercept = 50, linetype = "dashed")
    }
    g
  }

  if(PCR) {
    plot_grid <- expand_grid(PCR = c(FALSE, TRUE), facet1 = c(FALSE, TRUE)) %>%
      mutate(filename = paste0(dir_name, "model_predictions_", strain1,
                               ifelse(PCR, "_PCR", ""),
                               ifelse(facet1, "_facet", ""), ".png"))
    g <- Map(function(x, y) plot_inner(data_df, trajectories, x, y), plot_grid$PCR, plot_grid$facet1)
    Map(function(x, y) ggsave(x, y, width = 5, height = 3),
        plot_grid$filename, g)
  } else {
    g <- plot_inner(data_df, trajectories, FALSE)
    ggsave(paste0(dir_name, "model_predictions_", strain1, ".png"), g, width = 5, height = 5)
  }

}

plot_trajectories_camostat <- function(strain1, Calu3 = FALSE, trajectory_filename) {
  
  dir_name <- paste0(strsplit(trajectory_filename, "/", fixed = TRUE)[[1]][1], "/")
  
  trajectories <- readRDS(trajectory_filename)
  PCR <- "V_tot" %in% colnames(trajectories)
  
  if(Calu3) {
    data_df <- read_Calu3_data(strain = strain1)
  } else {
    data_df <- read_camostat_data(strain = strain1)
  }
  
  if(PCR) {
    data_df <- data_df %>%
      pivot_longer(starts_with("V"), names_to = "model", values_to = "V") %>%
      mutate(PCR = grepl("tot", model),
             model = gsub("_tot", "", model),
             model = recode(model, V = "full",
                            V_Camostat = "endosomal_only",
                            V_AmphoB = "full_no_IFITM",
                            V_Camostat_AmphoB = "endosomal_no_IFITM"))
  } else {
    data_df <- data_df %>%
      rename(full = "V_no_drug", endosomal_only = "V_Camostat") %>%
      pivot_longer(c("full", "endosomal_only"), names_to = "model", values_to = "V")
  }
  
  trajectories <- readRDS(trajectory_filename)
  
  data_df <- data_df %>%
    filter(t >= 0)
  plot_inner <- function(data_df, trajectories, PCR, facet1 = FALSE) {
    
    if(PCR) {
      trajectories <- trajectories %>%
        select(t, V_tot, probs, model) %>%
        pivot_wider(names_from = probs, values_from = V_tot) %>%
        filter(t >= 0)
      y_label = "RNA copy number"
      y_max <- 1e12
      data_df <- data_df %>%
        filter(PCR)
    } else {
      trajectories <- trajectories %>%
        select(t, V, probs, model) %>%
        pivot_wider(names_from = probs, values_from = V) %>%
        filter(t >= 0)
      if("PCR" %in% colnames(data_df)) {
        data_df <- data_df %>% filter(PCR == FALSE)
      }
      y_label = "Viral load (pfu)"
      y_max <- 1e8
    }
    
    levels_vec <- c("full", "endosomal_only", "tmprss2_only")
    
    data_df <- data_df %>%
      filter(model %in% levels_vec) %>%
      mutate(model = factor(model, levels = levels_vec))
    trajectories <- trajectories %>%
      filter(model %in% levels_vec) %>%
      mutate(model = factor(model, levels = levels_vec))
    facet_labels <- c("Both pathways", "Endosomal", "TMPRSS2")
    names(facet_labels) <- levels_vec
    g <- ggplot(data_df) +
      geom_ribbon(data = trajectories, aes(x = t, ymin = lower, ymax = upper, fill = model, group = model),
                  alpha = .3) +
      geom_point(aes(x = t, y = V,
                     color = model, group = model)) +
      geom_line(data = trajectories, aes(x = t, y = max_LL, color = model, group = model)) +
      scale_y_log10(y_label) +
      coord_cartesian(ylim = c(1, y_max)) +
      theme_bw() +
      xlab("Time (hpi)") +
      scale_color_discrete(drop = FALSE,
                           labels = c("No drug", "Camostat", ""),
                           name = "Drug") +
      scale_fill_discrete(guide = "none")
    if(facet1) {
      g <- g +
        facet_wrap(~model, nrow = 3, labeller = labeller(model = facet_labels))
    }
    if(!PCR) {
      g <- g + geom_hline(yintercept = 50, linetype = "dashed")
    }
    g
  }
  
  if(PCR) {
    plot_grid <- expand_grid(PCR = c(FALSE, TRUE), facet1 = c(FALSE, TRUE)) %>%
      mutate(filename = paste0(dir_name, "model_predictions_", strain1,
                               ifelse(PCR, "_PCR", ""),
                               ifelse(facet1, "_facet", ""), ".png"))
    g <- Map(function(x, y) plot_inner(data_df, trajectories, x, y), plot_grid$PCR, plot_grid$facet1)
    Map(function(x, y) ggsave(x, y, width = 3, height = 3),
        plot_grid$filename, g)
  } else {
    g <- plot_inner(data_df, trajectories, FALSE)
    ggsave(paste0(dir_name, "model_predictions_", strain1, ".png"), g, width = 5, height = 5)
  }
  
}
