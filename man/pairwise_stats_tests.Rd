% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/camostat_model_sum_stats.R
\name{pairwise_stats_tests}
\alias{pairwise_stats_tests}
\title{test whether the ratio of a summary statistic between two strains
differs from 1}
\usage{
pairwise_stats_tests(data_dir, par_name)
}
\arguments{
\item{data_dir}{character vector of length 2.  directories where results for the two strains live}

\item{par_name}{parameter to test}
}
\value{
a list of length n, where n is the number of pair combinations of data_dir.
each element contains a list with the elements prctile_ratios and p_value.
prctile_ratio is a vector which gives the 2.5, 50 and 97.5th percentiles of the ratio of each
p_value is a scalar where each element is twice the proportion of the
sampled ratios greater than 1, or the proportion of the sampled ratios smaller
than 1 (whichever is the smaller)
}
\description{
test whether the ratio of a summary statistic between two strains
differs from 1
}
