#' Emulates the ggplot colour palette.
#'
#' \code{gg_color_hue} emulates the ggplot colour palette.
#'
#' @param n numeric vector of length 1. The number of colours required.
#' @return a character vector of length \code{n}. Each element is a hex colour code.
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#' Specifies the \code{y} scale of a \code{ggplot} object.
#'
#' \code{add_y_scale} specifies the \code{y} scale of a plot.
#'
#' Called by \code{\link{plot.results.pyro}} to unify plots.
#' @param in_g_obj \code{ggplot} object to specify scale for.
#' @param log_scale logical vector of length 1. If true, use log scale on
#' \code{y}-axis; use linear scale otherwise.
#' @param y_label character vector of length 1. Y-axis label.
#' @param y_max numeric vector of length 1. Upper limit of y-axis.
#' @return a \code{ggplot} object: \code{in_g_obj} with specified \code{y}-scale.
add_y_scale <- function(in_g_obj, log_scale = TRUE, y_label, y_max, y_min){
  if(log_scale) {
    scale_fn <- "scale_y_log10"
    if(missing(y_min)) {
      y_min <- 1
    }
  } else {
    scale_fn <- "scale_y_continuous"
    if(missing(y_min)) {
      y_min <- 0
    }
  }
  in_g_obj + do.call(scale_fn, list(y_label, expand = c(0,0))) +
    ggplot2::coord_cartesian(ylim = c(y_min,y_max))
}

#' plot the marginal histogram of a single fitted variable
#'
#' @param chain a single row of a data table containing the parameter values
#' @param x_lab parameter name to display on x-axis
#' @param xlim numeric vector of length 2.  limits of x-axis
#' @return ggplot object
plot_marginal_histogram <- function(chain, x_lab, xlim){
  vec <- unname(unlist(chain))
  plot_histogram_general(vec = vec, limits = xlim, plot_xlab = x_lab)
}

#' plotting the marginal histogram of all fitted variables
#'
#' @param chain data table containing iterations of MCMC chain
#' @param parTab data frame containing parameters
#' @param ncol if ncol = 0, output each plot as a separate ggplot object.
#' if ncol > 0, arrange into ncol columns in a single ggplot object.
#' @return ggplot object or a list of ggplot objects
plot_marginal_histogram_all <- function(chain, parTab, ncol = 1){
  unfixed_pars <- which(parTab$fixed == 0 & !is.na(parTab$values))

  g_list <- lapply(unfixed_pars, function(x) plot_marginal_histogram(
    chain = chain[,parTab[x,"names"],drop = FALSE,with = FALSE],
    x_lab = parTab[x,"names_plot"],
    xlim = c(parTab[x,"lower_bound"],parTab[x,"upper_bound"])
  ))

  arrange_grid(g_list, ncol)

}

#' make a trace plot of a single fitted variable
#'
#' @param chain data table containing iterations of MCMC chain after discarding
#' burn-in
#' @param true_value numeric vector of length 1. The true parameter value (if
#' unknown because we're using real data, set this to any value)
#' @param start_value Parameter value at start of MCMC chain (before discarding burn-in)
#' @param var_name name of parameter to plot in colnames(chain)
#' @param y_lab parameter name to display on plot
#' @param ylim if zoom = FALSE, numeric vector of length 2 giving lower and
#' upper limits of y axis. If either of these is NA, the range of the parameter
#' values in the chain + the start value is used. if zoom = TRUE, this parameter
#' is not used.
#' @param zoom logical vector of length 1. if TRUE, the y-axis limits are the
#' range of the parameter values in the chain (the start value is disregarded).
#' if FALSE, the limits are as descibed in ylim.
#' @param n_samples numeric vector of length 1.  number of samples to plot, taken
#' at evenly spaced intervals along the chain.  if 0, plot all samples.
#' @import ggplot2
#' @return ggplot object
plot_trace <- function(chain, true_value, start_value, var_name, y_lab, ylim, zoom,
                       n_samples = 0, real_data){
  chain <- thin_chain(chain,n_samples)
  space <- FALSE

  if(is.na(true_value)){
    true_value <- numeric(0)
  }
  if(any(is.na(ylim))){
    all_values <- c(as.matrix(chain[,var_name,with=FALSE]), true_value, start_value)
    ylim <- range(all_values)
    space <- TRUE
  }
  if(zoom){
    ylim <- range(c(as.matrix(chain[,var_name,with=FALSE]),true_value))
    space <- TRUE
  }

  g <- ggplot(chain, aes_string("sampno", var_name)) +
    geom_line() +
    scale_x_continuous("Iteration", expand = c(0,0),
                       labels = scales::scientific) +
    theme_bw()
  if(length(true_value) == 0){
    g <- g + geom_hline(yintercept = start_value, colour = gg_color_hue(2)[1]) +
      theme(text = element_text(size = 32))
  } else {
    hlines <- data.frame(yint = c(true_value,start_value),Value = c("true", "start"))
    g <- g + geom_hline(data = hlines,aes(yintercept = yint,colour = Value, linetype = Value)) +
      theme(legend.justification=c(1,1),
            legend.position=c(1,1),
            text = element_text(size = 32))
  }
  if(space){
    g <- g + scale_y_continuous(latex2exp::TeX(y_lab), limits = ylim,
                                breaks=ylim)
  } else {
    g <- g + scale_y_continuous(latex2exp::TeX(y_lab), limits = ylim,
                                breaks=ylim,
                                expand = c(0,0))
  }
  g
}

#' plot the trace of all fitted variables
#'
#' @param chain data table containing iterations of MCMC chain after discarding
#' burn-in
#' @param parTab data frame containing parameter information, created by
#' specify_parameters()
#' @param startTab data frame of the same form as parTab, containing the
#' parameter values at the start of the chain (before discrading burnin)
#' @param zoom logical vector of length 1. if TRUE, the y-axis limits are the
#' range of the parameter values in the chain (the start value is disregarded).
#' if FALSE, the limits are as given in parTab, unless at least one of these is NA,
#' in which case the range of the chain values + start value is used.
#' @param ncol if ncol = 0, output each plot as a separate ggplot object.
#' if ncol > 0, arrange into ncol columns in a single ggplot object.
#' @param n_samples numeric vector of length 1.  number of samples to plot, taken
#' at evenly spaced intervals along the chain.  if 0, plot all samples.
#' @param real_data logical vector ot length 1.  If TRUE, because we're using
#' real data, we don't know the true parameter values, so the values in parTab
#' are just example values.  If FALSE, we are using simulated data, and the
#' values in parTab are the true values, so plot these
#' @return ggplot object or a list of ggplot objects
plot_trace_all <- function(chain, parTab, startTab, zoom, ncol = 1,
                           n_samples = 0, real_data){

  unfixed_pars <- which(parTab$fixed == 0)
  if(real_data) {
    true_value <- rep(NA, nrow(parTab))
  } else {
    true_value <- parTab[,"values"]
  }

  g_list <- lapply(unfixed_pars, function(x) plot_trace(
    chain = chain[,c("sampno",as.character(parTab[[x,"names"]])), with = FALSE],
    true_value = true_value[x],
    start_value = startTab[x,"values"],
    var_name = parTab[x,"names"],
    y_lab = parTab[x,"names_plot"],
    ylim = c(parTab[x,"lower_bound"],parTab[x,"upper_bound"]),
    zoom = zoom,
    n_samples = n_samples))

  arrange_grid(g_list, ncol)
}

#' arrange a list of ggplot objects into a grid in a single ggplot object
#'
#' @param g_list list of ggplot objects
#' @param ncol if ncol = 0, this function returns the original g_list.
#' if ncol > 0, arrange g_list into ncol columns in a single ggplot object.
#' @return either a single ggplot object or a list of ggplot objects
arrange_grid <- function(g_list, ncol){
  if(ncol == 0){
    invisible(g_list)
  } else {
    g_list$ncol <- ncol
    g2 <- do.call(grid.arrange,g_list)
    invisible(g2)
  }
}

#' plot the potential scale reduction factor vs the number of MCMC iterations
#'
#' @param iterations numeric voector of length 1. the potential scale reduction
#' factor is calculated once every this many iterations
#' @param max_psrf numeric vector.  The maximum PSRF among the parameters each
#' time the PSRF is calculated
#' @return ggplot object
plot_diagnostics <- function(max_psrf,iterations){
  plot_df <- data.frame(iterations = iterations*seq_along(max_psrf),max_psrf = max_psrf)
  g <- ggplot(data = plot_df, aes(x = iterations, y = max_psrf)) +
    geom_point() +
    scale_y_continuous("Max PSRF", expand = c(0,0), limits = c(0,max(plot_df$max_psrf))) +
    scale_x_continuous("Iterations", expand = c(0,0), label=scientific_format(digits=2)) +
    geom_hline(yintercept = 1.1, linetype = "dashed", size = 1) +
    theme_bw() +
    theme(text = element_text(size = 24))
  g
}

#' make scatter plot of samples from joint posterior distribution across two parameters
#'
#' @param chain data table containing parameter values
#' @param true_value If simulated data is used, numeric vector of length 2
#' containing the true parameter values, otherwise c(NA, NA)
#' @param start_value nx2 matrix of starting parameter values, where n is the
#' number f parallel MCMC chains
#' @param var_name character vector of length 2.  Parameter names in colnames(chain)
#' @param lower_bounds numeric vector of length 2. Lower axis limits.
#' @param upper_bounds numeric vector of length 2.  Upper axis limits.
#' @param labels character vector of length 2. Parameter names to display on plot
#' @param n_samples numeric vector of length 1.  number of samples to plot, taken
#' at evenly spaced intervals along the chain.  if 0, plot all samples.
#' @param zoom logical vector of length 1. if TRUE, the axis limits are the
#' range of the parameter values in the chain (the start value is disregarded).
#' if FALSE, the limits are as given in lower_bounds and upper_bounds,
#' unless at least one of these is NA,
#' in which case the range of the chain values + start value is used.
#' @return a ggplot object
#' @import ggplot2
plot_bivariate_scatter <- function(chain, true_value, start_value, var_name,
                                   lower_bounds, upper_bounds,
                                   labels, n_samples = 0, zoom = FALSE){
  chain <- thin_chain(chain,n_samples)

  real_data <- any(is.na(true_value))
  if(!real_data){
    chain <- rbind(chain,as.list(true_value))
  }

  start_data_table <- data.table::as.data.table(start_value)
  # if start_value is a vector rather than a matrix, it gets transposed
  # automatically -- transpose back
  if(ncol(start_data_table) == 1){
    start_data_table <- t(start_data_table)
  }
  colnames(start_data_table) <- colnames(chain)

  if(real_data){
    colours <- gg_color_hue(3)
    if(zoom){
      xlim <- range(chain[,1])
      ylim <- range(chain[,2])
    } else {
      xlim <- c(lower_bounds[1], upper_bounds[1])
      ylim <- c(lower_bounds[2], upper_bounds[2])
    }

    g <- ggplot(data = chain, aes_string(x = var_name[1],y = var_name[2])) +
      geom_point() +
      scale_x_continuous(latex2exp::TeX(labels[1]), limits = xlim,
                         breaks=xlim) +
      scale_y_continuous(latex2exp::TeX(labels[2]), limits = ylim,
                         breaks=ylim) +
      theme_bw() +
      theme(legend.justification=c(1,1),
            legend.position = "none",
            text = element_text(size = 24))

  } else {
    chain <- rbind(chain,start_data_table)
    chain <- cbind(chain, data.frame(Value = c(rep("chain",nrow(chain)-nrow(start_data_table)-1),
                                               "true",paste0("start",1:nrow(start_data_table)))))
    levels(chain$Value) <- unique(c(levels(chain$Value),"start"))
    chain$Value[(nrow(chain)-nrow(start_value)+1)] <- "start"
    chain$Value <- factor(chain$Value, levels = sort(levels(chain$Value)))
    colours <- gg_color_hue(3)
    if(zoom){
      xlim <- range(chain[1:(nrow(chain)-nrow(start_data_table)),1])
      ylim <- range(chain[1:(nrow(chain)-nrow(start_data_table)),2])
    } else {
      xlim <- c(lower_bounds[1], upper_bounds[1])
      ylim <- c(lower_bounds[2], upper_bounds[2])
    }

    g <- ggplot(data = chain, aes_string(x = var_name[1],y = var_name[2])) +
      geom_point(aes(colour = Value, size = Value)) +
      scale_colour_manual(breaks = c("chain","start","true"),
                          values = c(colours[1],rep(colours[2],(nrow(start_data_table))),colours[3])) +
      scale_size_manual(breaks = c("chain","start","true"),
                        values = c(1,rep(3,(nrow(start_data_table)+1)))) +
      scale_x_continuous(latex2exp::TeX(labels[1]), limits = xlim,
                         breaks=xlim) +
      scale_y_continuous(latex2exp::TeX(labels[2]), limits = ylim,
                         breaks=ylim) +
      theme_bw() +
      theme(legend.justification=c(1,1),
            legend.position = "none",
            text = element_text(size = 24))
  }
  g
}

#' plot bivariate scatter plots of all fitted variables
#'
#' @param chain data table containing MCMC iterations after discarding burn-in
#' @param parTab data frame containing parameter information constructed by
#' specify_parameters()
#' @param startTab either a data frame of the same form as parTab, containing
#' starting parameter values (before discarding burn-in), or a list of such data
#' frames, one for each parallel chain
#' @param ncol if ncol = 0, output each plot as a separate ggplot object.
#' if ncol > 0, arrange into ncol columns in a single ggplot object.
#' @param n_samples numeric vector of length 1.  number of samples to plot, taken
#' at evenly spaced intervals along the chain.  if 0, plot all samples.
#' @param zoom logical vector of length 1. if TRUE, the axis limits are the
#' range of the parameter values in the chain (the start value is disregarded).
#' if FALSE, the limits are as given in parTab,
#' unless at least one of these is NA,
#' in which case the range of the chain values + start value is used.
#' @param real_data logical vector ot length 1.  If TRUE, because we're using
#' real data, we don't know the true parameter values, so the values in parTab
#' are just example values.  If FALSE, we are using simulated data, and the
#' values in parTab are the true values, so plot these
#' @return ggplot object or a list of ggplot objects
plot_bivariate_scatter_all <- function(chain, parTab, startTab,
                                       ncol = 1, n_samples = 0, zoom,
                                       real_data = real_data){
  unfixed_pars <- which(parTab$fixed == 0)
  parTab <- parTab[unfixed_pars,]
  # if I've passed a list of data frames instead of a single data frame
  if(!is.data.frame(startTab)){
    startTab <- lapply(startTab, function(x) x[unfixed_pars,])
  } else {
    startTab <- startTab[unfixed_pars,]
  }

  idx <- expand.grid(seq_along(unfixed_pars), seq_along(unfixed_pars))

  if(!is.data.frame(startTab)){
    start_value <- t(vapply(startTab, function(x) x$values, double(nrow(startTab[[1]]))))
  } else {
    start_value = t(as.matrix(startTab$values))
  }

  if(real_data) {
    true_value <- rep(NA, nrow(parTab))
  } else {
    true_value <- parTab[,"values"]
  }
  plot_wrapper <- function(x) {
    plot_bivariate_scatter(chain = chain[,parTab[x,"names"], with = FALSE],
                           true_value = true_value[x],
                           start_value = start_value[,x],
                           var_name = parTab[x,"names"],
                           lower_bounds = parTab[x,"lower_bound"],
                           upper_bounds = parTab[x,"upper_bound"],
                           labels = parTab[x,"names_plot"],
                           n_samples = n_samples,
                           zoom = zoom)
  }

  g_list <- apply(idx, 1, plot_wrapper)

  arrange_grid(g_list, ncol)
}

#' select evenly spaced iterations from MCMC chain
#'
#' @param chain data table containing MCMC iterations
#' @param n_samples numeric vector of length 1. number of samples to select
#' @return data table of selected iterations
thin_chain <- function(chain,n_samples){
  if(n_samples == 0 || n_samples >= nrow(chain)){
    chain
  } else {
    chain[round(seq(1,nrow(chain),length.out = n_samples)),]
  }
}

#' plot a nicer-looking histogram
#'
#' @param vec numeric vector of valuse from which to make histogram
#' @param special_value optional parameter.  If specified, should be numeric
#' vector of length 1. Plot the bar for this value in a different colour. Used,
#' for example, to draw attention to the bar at x = 0.
#' @param breaks optional parameter.  If specified, should be numeric vector.
#' Specifies breaks on horizontal axis.  If not given, use horizontal axis
#' limits as breaks.
#' @param vline_values optional parameter.  If specified, should be numeric
#' vector. Plot vertical lines at these values.
#' @param limits optional parameter. If specified, should be numeric vector of
#' length 2.  Use as horizontal axis limits. If not given use the range of vec
#' (if there are values below 0) or [0, max(vec)] (if there are no values below
#' zero).
#' @param plot_title optional parameter. If specified, should be character
#' vector of length 1. Title of plot.
#' @param plot_xlab character vector of length 1. Horizontal axis label.
#' @import ggplot2
#' @return ggplot object
plot_histogram_general <- function(vec, special_value, breaks, vline_values, limits, plot_title, plot_xlab){

  if(missing(breaks)){
    breaks <- limits
  }

  if(missing(limits)){
    if(any(vec < 0)){
      limits <- range(vec)
    } else {
      limits <- c(0,max(vec))
    }
  }

  breaks_list <- gen_nice_breaks(breaks, limits)

  if(breaks_list$order_magnitude == 1){
    cat_xlab_string <- " ($\\times 10$)"
  } else if (breaks_list$order_magnitude == 0) {
    cat_xlab_string <- ""
  } else {
    cat_xlab_string <- paste0(" ($\\times 10^{",breaks_list$order_magnitude,"}$)")
  }

  if(missing(special_value)){
    plot_df <- data.frame(vec = vec)
  } else {
    plot_df <- data.frame(vec = vec, is_special = vec == special_value)
  }

  g <- ggplot(plot_df, aes(x = vec)) +
    theme_bw() +
    scale_x_continuous(latex2exp::TeX(paste0(plot_xlab, cat_xlab_string)),
                       limits = limits,
                       breaks = breaks_list$breaks,
                       labels = breaks_list$break_names,
                       expand = c(0,0)) +
    scale_y_continuous("",breaks = NULL)

  if(!missing(plot_title)){
    g <- g + ggtitle(plot_title)
  }

  if(missing(special_value)){
    g <- g + geom_histogram(bins = 30) +
      theme(text = element_text(size = 32), aspect.ratio = 1)
  } else {
    g <- g + geom_histogram(aes(fill = is_special), bins = 30) +
      theme(text = element_text(size = 32), aspect.ratio = 1, legend.position = "none")
  }

  if(!missing(vline_values)) {
    g <- g + geom_vline(xintercept = vline_values, linetype = "dotted")
  }

  g
}

#' construct plot filename
#'
#' @param filename character vector of length 1. experiment/chain identifier
#' @param par_name character vector of length 1. parameter name
#' @param description character vector of length 1. description of plot
#' @return character vector of length 1.
make_filename <- function(filename, par_name, description) {
  paste0(filename,"_", description, par_name,".pdf")
}

#' darken a vector of hex codes for colours
#'
#' @param colour_in vector of hex codes for colours
#' @return vector of hex codes for darkened colours
darken <- function(colour_in) {
  offset_const <- -.3
  grDevices::adjustcolor(colour_in, offset = c(rep(offset_const, 3), 0))
}

#' lighten a vector of hex codes for colours
#'
#' @param colour_in vector of hex codes for colours
#' @return vector of hex codes for lighened colours
lighten <- function(colour_in) {
  offset_const <- .3
  grDevices::adjustcolor(colour_in, offset = c(rep(offset_const, 3), 0))
}

#' plot trace of log likelihood for each chain
#'
#' @param filename string containing filename of .RData file outputted by lazymcmc
#' @return ggplot object
#' @import ggplot2
plot_LL <- function(filename) {
  chain <- get_MCMC_chain(filename, raw = TRUE)

  process_chain <- function(chain_idx) {
    chain <- as.data.frame(chain[[chain_idx]]) %>%
      extract(c("sampno", "lnlike")) %>%
      cbind(., data.frame(idx = chain_idx))
  }

  chain <- lapply(seq_along(chain), process_chain) %>%
    do.call(rbind, .)
  chain$idx <- factor(chain$idx)
  ggplot(chain, aes(x = sampno, y = lnlike, color = idx, group = idx)) +
    geom_line() +
    theme_bw() +
    theme(legend.justification=c(1,1),
          legend.position=c(.9,.4),
          text = element_text(size = 16)) +
    xlab("Iteration") +
    ylab("Log likelihood") +
    guides(color = guide_legend(title="Chain"))
}

#' plot virus dynamics
#'
#' @return ggplot object
#' @import ggplot2
plot_dynamics_inner <- function(data_df, ci_df, samples_df, noise_var, true_var, y_label, error_string,
                                y_min, y_max, obs_threshold, x_label = "Time (hpi)"){
  epsilon <- 1e-10
  grey <- "#D8D8D8"
  point_shape <- 24
  point_size <- 3

  if(!missing(obs_threshold) && !("below_threshold" %in% colnames(data_df))) {
    data_df$below_threshold <- data_df[,noise_var] < obs_threshold
    data_df[data_df$below_threshold & !is.na(data_df$below_threshold), noise_var] <- 1
  }

  g <- ggplot(data=data_df) +
    xlab(x_label) + ggtitle("") +
    scale_x_continuous(expand = c(0,0)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))

  if(!missing(error_string) && !is.null(ci_df)){

    ci_df <- cbind(ci_df, data.frame("below_1" = ci_df[,error_string[1]] < 1,
                                     "above_max" = ci_df[,error_string[2]] > y_max))
    ci_df <- cbind(ci_df, data.frame("out_range" = ci_df[,"below_1"] | ci_df[,"above_max"]))
    ci_df[which(ci_df$below_1),error_string[1]] <- rep((1 + epsilon),sum(ci_df$below_1, na.rm = TRUE))
    ci_df[which(ci_df$above_max),error_string[2]] <- rep((y_max - epsilon),sum(ci_df$above_max, na.rm = TRUE))

    g <- g + geom_ribbon(data = ci_df[!is.na(ci_df[,error_string[1]]),],
                         mapping = aes_string(x = "t", ymin = error_string[1],
                                              ymax = error_string[2]), fill = grey)
  }

  if(!missing(samples_df) && !is.null(samples_df)){

    # samples_df <- melt(samples_df, id.vars = "t", na.rm = TRUE)
    g <- g + geom_line(data = samples_df, aes(x = t, y = value, color = variable),
                       show.legend = FALSE)
  }

  if(true_var %in% colnames(data_df)){
    g <- g + geom_line(data = data_df[!is.na(data_df[,true_var]),],
                       mapping = aes_string(x = "t", y = true_var))
  }

  if(!missing(noise_var) && noise_var %in% colnames(data_df)){
    if("below_threshold" %in% colnames(data_df)){
      g <- g + geom_point(mapping = aes_string(x = "t", y = noise_var,
                                               fill = "below_threshold"),
                          colour = "black",
                          shape = point_shape, size = point_size)
    } else {
      g <- g + geom_point(mapping = aes_string(x = "t", y = noise_var),
                          fill = "black",
                          colour = "transparent",
                          shape = point_shape, size = point_size)
    }
  }

  if(!missing(obs_threshold)){
    g <- g + geom_hline(yintercept = obs_threshold, linetype = "dotted") +
      scale_fill_manual(breaks=c(FALSE,TRUE),values=c("black", "white"))

    if(sum(data_df$below_threshold, na.rm = TRUE) == 0) { # if no values below threshold
      g <- g + theme(legend.position="none")
    } else {
      g <- g + theme(legend.justification=c(1,1),
                     legend.position=c(1,1))
    }
  } else {
    g <- g + theme(legend.position="none")
  }
  g <- add_y_scale(g, log_scale = TRUE,  y_label = y_label, y_min = y_min, y_max = y_max)
  g
}

#' plot marginal histogram, either coloured or faceted by chain
#'
#' @param dir_name string: directoy where results are saved
#' @param par_name string: parameter to plot
#' @return list of 2 ggplot objects.  First is marginal histogram coloured by chain,
#' second is marginal histogram faceted by chain
plot_marginal_by_chain <- function(dir_name, par_name) {
  chains <- get_MCMC_chain(paste0(dir_name, "1.RData"), bind = FALSE)

  make_par_tibble <- function(chain_no, par_name) {
    tibble(chain_no = chain_no, values = unlist(chains[[chain_no]][,par_name]))
  }

  par_tibble <- lapply(seq_along(chains), make_par_tibble, par_name = par_name) %>%
    bind_rows() %>%
    mutate(chain_no = factor(chain_no))

  g <- ggplot(par_tibble, aes(x = values, fill = chain_no, group = chain_no)) +
    geom_histogram(position = "dodge") +
    theme_bw() +
    xlab(par_name)
  ggsave(paste0(dir_name, par_name, ".png"), width = 3, height = 3)

  h <- ggplot(par_tibble, aes(x = values)) +
    facet_wrap(~chain_no, nrow = 1) +
    geom_histogram() +
    theme_bw() +
    xlab(par_name)
  ggsave(paste0(dir_name, par_name, "_individual.png"), width = 8, height = 3)
  list(g, h)
}
