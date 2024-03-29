% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/forwardBackward.R
\name{forward_backward}
\alias{forward_backward}
\title{Forward and Backward Probabilities for Hidden Markov Model}
\usage{
forward_backward(model, forward_only = FALSE, log_space = FALSE, threads = 1)
}
\arguments{
\item{model}{Object of class \code{hmm} or \code{mhmm}.}

\item{forward_only}{If \code{TRUE}, only forward probabilities are computed. The default is \code{FALSE}.}

\item{log_space}{Compute forward and backward probabilities in logarithmic scale instead of scaling.
The default is \code{FALSE}.}

\item{threads}{Number of threads used in parallel computing. The default is 1.}
}
\value{
List with components
  \item{forward_probs}{If \code{log_space = FALSE}, scaled forward probabilities, i.e. probability of state given
  observations up to that time point. If \code{log_space = TRUE},
  logarithms of non-scaled forward probabilities. }
  \item{backward_probs}{Scaled backward probabilities (\code{log_space = FALSE}),
  or logarithms of non-scaled backward probabilities(\code{log_space = TRUE}). }
  \item{scaling_factors}{Sum of non-scaled forward probabilities at each time point.
  Only computed if \code{log_space = FALSE}.}
  In case of multiple observations, these are computed independently for each sequence.
}
\description{
The \code{forward_backward} function computes scaled forward and backward probabilities of a hidden Markov model.
}
\examples{
# Load a pre-defined MHMM
data("mhmm_biofam")

# Compute forward and backward probabilities
fb <- forward_backward(mhmm_biofam)

# The most probable hidden state at time t
# given the observations up to time t for the first subject:
apply(fb$forward_probs[, , 1], 2, which.max)

}
