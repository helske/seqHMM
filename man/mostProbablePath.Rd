\name{mostProbablePath}
\alias{mostProbablePath}
\title{Most Probable Path of Hidden States of Hidden Markov Model given the sequence.}
\usage{
mostProbablePath(model, tol = 1e-15)
}
\arguments{
  \item{model}{Hidden Markov model of class
  \code{HMModel}.}

  \item{tol}{Tolerance parameter defining the smallest
  possible nonzero probability. Used for ensuring that
  there won't be infinite log-probabilities. Default is
  1e-15.}
}
\value{
Most probable path of states given the observations. In
case of multiple observations, most probable path is
computed independently for each sequence.
}
\description{
Function \code{mostProbablePath} computes the most probable
path of the hidden states of the hidden Markov model given
the observed sequence.
}

