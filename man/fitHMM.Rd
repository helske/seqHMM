\name{fitHMM}
\alias{fitHMM}
\title{Estimate parameters of Hidden Markov Model}
\usage{
fitHMM(model, method = "nlm", ...)
}
\arguments{
  \item{model}{Hidden Markov model of class HMModel.}

  \item{...}{Further arguments to nlm, such as arguments
  controlling the amount printing and number of
  iterations.}
}
\value{
List with two components: First component is output from
nlm and second is the estimated model.
}
\description{
Function \code{fitHMM} estimates the initial state,
transition and emission probabilities of hidden Markov
model using numerical maximization of log-likelihood.
Initial values for estimation are taken from the
corresponding components of model with preservation of
original zero probabilities.
}

