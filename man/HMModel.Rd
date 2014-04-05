\name{HMModel}
\alias{HMModel}
\title{Build a Hidden Markov Model}
\usage{
HMModel(observations, transitionMatrix, emissionMatrix, initialProbs,
  stateNames = NULL)
}
\arguments{
  \item{observations}{Observations. Can be either vector
  (single sequence ), data.frame (multiple sequences, each
  row corresponds to one sequence), or list of
  vectors/data.frames where each component corresponds to
  one channel.}

  \item{transitionMatrix}{Matrix of transition
  probabilities.}

  \item{emissionMatrix}{Matrix of emission probabilities or
  in case multiple channels list of such matrices.}

  \item{initialProbs}{Vector of initial state
  probabilities.}

  \item{stateNames}{Optional labels for the observations.}
}
\value{
Object of class \code{HMModel}
}
\description{
Function HHModel constructs an object of class
\code{HMModel}.
}

