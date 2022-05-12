# make a vector of prediction times from a vector of sampling times
make_prediction_time <-
  function(sampling_time,
           min_prediction_period,
           t_end) {
    sort(unique(c(
      sampling_time, seq(0, t_end, by = min_prediction_period)
    )))
  }

#' make a vector of solving times from a vector of times
#'
#' @param times_in numeric vector of times
#' @return numeric vector: sorted unique times, with t = 0
make_solving_time <- function(times_in) {
  sort(unique(c(0, times_in)))
}
