\name{fitHMM}
\alias{fitHMM}
\title{Estimate parameters of Hidden Markov Model}
\usage{
fitHMM(model, ...)
}
\arguments{
  \item{model}{Hidden Markov model of class
  \code{HMModel}.}

  \item{...}{Further arguments to \code{nlm}.}
}
\value{
List with two components: First component is output from
\code{nlm} and second is the estimated model.
}
\description{
Function \code{fitHMM} estimates the initial state,
transition and emission probabilities of hidden Markov
model using numerical maximization of log-likelihood.
Initial values for estimation are taken from the
corresponding components of model with preservation of
original zero probabilities.
}

