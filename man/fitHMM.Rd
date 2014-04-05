\name{fitHMM}
\alias{fitHMM}
\title{Estimate parameters of Hidden Markov Model}
\usage{
fitHMM(model, method = c("bobyqa", "nlm"), ...)
}
\arguments{
  \item{model}{Hidden Markov model of class
  \link{\code{HMModel}}.}

  \item{method}{Function used in optimization. Default is
  \link{\code{bobyqa}}, a derivative free method. Another
  option is \code{nlm} which uses Newton-type algorithm.}

  \item{...}{Further arguments to \link{\code{nlm}} or
  \link{\code{bobyqa}}, such as arguments controlling the
  amount printing and number of iterations. See
  corresponding functions for details.}
}
\value{
List with two components: First component is output from
\link{\code{nlm}} and second is the estimated model.
}
\description{
Function \code{fitHMM} estimates the initial state,
transition and emission probabilities of hidden Markov
model using numerical maximization of log-likelihood.
Initial values for estimation are taken from the
corresponding components of model with preservation of
original zero probabilities.
}

