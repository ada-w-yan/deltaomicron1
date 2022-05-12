plot_postprocessed_chain <-function(postprocessed_list,
                                    plot_dynamics = NULL,
                                    plot_residuals = NULL,
                                    plot_flags = c(diagnostics = TRUE,
                                                   priors = TRUE,
                                                   posteriors = TRUE, 
                                                   bivariate = FALSE)){
  
  list2here(postprocessed_list, 
            var_names = c("chain", "parTab", "startTab", "prior_chain", 
                          "mcmcPars", "filenames", "data_df", 
                          "residuals"))
  if("diagnostics" %in% names(postprocessed_list) && plot_flags[["diagnostics"]]) {
    diagnostics <- postprocessed_list$diagnostics
      g <- plot_diagnostics(diagnostics$max_psrf, mcmcPars[["iterations"]]/10)
      
      try(ggsave_wch(paste0(filenames[1],"_diagnostics.pdf"), g))
  }
  
  remove(postprocessed_list)
  

  
  if(plot_flags[["priors"]]) {
    unfixed_pars <- parTab$fixed == 0 & !is.na(parTab$values)
    par_names_unfixed <- parTab$names[unfixed_pars]
    g <- plot_marginal_histogram_all(prior_chain, parTab, ncol = 0)
    
    ggsave_wrapper(g,
                   filename_fn = function(y) make_filename(filenames[1], 
                                                           y, 
                                                           "marginal_prior"),
                   width = 10,
                   height = 10,
                   filename_args = list(par_names_unfixed))
  }

  plot_posterior_partial <- function(startTab, chain, filenames, residuals, combined_flag) {
    plot_posterior(parTab, startTab, data_df, chain,
                   filenames, combined_flag, plot_dynamics,
                   residuals, plot_residuals, plot_flags)
  }
  
  if(is.data.frame(chain)) {
    # converged
    h <- plot_posterior_partial(startTab, 
                                chain, 
                                filenames[1], 
                                residuals, 
                                combined_flag = TRUE)
  } else {
    h <- Map(plot_posterior_partial,
             startTab, chain, filenames, residuals, combined_flag = FALSE)
  }
  invisible(NULL)
}

plot_posterior <- function(parTab, startTab, data_df, chain,
                           filename, combined_flag, plot_dynamics,
                           calc_residuals, plot_residuals, plot_flags){
  
  if(combined_flag){
    combined_str <- "_combined"
  } else {
    combined_str <- ""
  }
  
  if(!is.null(plot_dynamics)) {
    
    g <- plot_dynamics(parTab, data_df, chain)

    if(is.null(names(g))) {
      filename_args <- list(seq_along(g))
    } else {
      filename_args <- list(names(g))
    }
    
    ggsave_wrapper(g,
                   filename_fn = function(x) make_filename(filename, x,
                                             paste0("dynamics",combined_str)),
                   width  = 10, height = 10,
                   filename_args = filename_args)
  }
  
  if(!is.null(plot_residuals) && !is.null(residuals)){
    g <- plot_residuals(residuals, data_df)
    
    ggsave_wrapper(g,
                   filename_fn = function(x) make_filename(filename, x,
                                             paste0("residuals",combined_str)),
                   width  = 10, height = 10,
                   filename_args = list(seq_along(g) - 1))
  }
  
  unfixed_pars <- parTab$fixed == 0 & !is.na(parTab$values)
  par_names_unfixed <- parTab$names[unfixed_pars]
  
  if(plot_flags[["posteriors"]]) {
    g <- plot_marginal_histogram_all(chain, parTab, ncol = 0)
    
    ggsave_wrapper(g,
                   filename_fn = function(x) make_filename(filename, x,
                                             paste0("marginal",combined_str)),
                   width  = 10, height = 10,
                   filename_args = list(par_names_unfixed))
  }
  
  
  ## plot two-parameter correlations
  if(plot_flags[["bivariate"]] && length(par_names_unfixed) > 1){
    
    plot_bivariate_scatter_partial <- 
      function(zoom) plot_bivariate_scatter_all(chain,
                                                parTab,
                                                startTab,
                                                ncol = 0,
                                                n_samples = 1e3,
                                                zoom = zoom,
                                                real_data = TRUE)
    g <- lapply(c(FALSE, TRUE), plot_bivariate_scatter_partial)
    names(g) <- paste0(c("bivariate", "bivariate_zoom"), combined_str)
    
    par_name_mat <- expand.grid(par_names_unfixed, par_names_unfixed)
    
    filenames_bivariate <- apply(par_name_mat, 1, paste0, collapse = "")
    
    ggsave_partial <- function(element_name, g) {
      idx_mat <- expand.grid(seq_along(par_names_unfixed), 
                             seq_along(par_names_unfixed))
      duplicate <- idx_mat[,2] - idx_mat[,1] >= 0
      ggsave_wrapper(g[!duplicate], 
                     filename_fn = function(x) make_filename(filename, 
                                                             x, 
                                                             element_name),
                     width = 10, 
                     height = 10, 
                     filename_args = list(filenames_bivariate[!duplicate]))
    }
    
    lapply_w_name(g, ggsave_partial)
  }
  invisible(NULL)
}