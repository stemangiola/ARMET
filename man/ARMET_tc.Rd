% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ARMET_tc.R
\name{ARMET_tc}
\alias{ARMET_tc}
\title{ARMET-tc main}
\usage{
ARMET_tc(
  .data,
  .formula = ~1,
  .sample = NULL,
  .transcript = NULL,
  .abundance = NULL,
  reference = NULL,
  approximate_posterior = F,
  iterations = 250,
  sampling_iterations = 100,
  .n_markers = n_markers,
  do_regression = T,
  prior_survival_time = c(),
  model = stanmodels$ARMET_tc_fix_hierarchical
)
}
\arguments{
\item{.data}{A tibble}

\item{.formula}{A formula}

\item{.sample}{A column symbol}

\item{.transcript}{A column symbol}

\item{.abundance}{A column symbol}

\item{reference}{A tibble}

\item{approximate_posterior}{A boolean for variational Bayes}

\item{iterations}{An integer total iterations}

\item{sampling_iterations}{An integer. Sampling iteractions}

\item{.n_markers}{A tibble}

\item{do_regression}{A boolean}

\item{prior_survival_time}{An array}

\item{model}{A stan model}
}
\value{
An ARMET object
}
\description{
This function calls the stan model.
}
