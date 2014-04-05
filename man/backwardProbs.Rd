\name{backwardProbs}
\alias{backwardProbs}
\title{Backward Probabilities for Hidden Markov Model}
\usage{
backwardProbs(model)
}
\arguments{
  \item{model}{Hidden Markov model of class
  \code{HMModel}.}
}
\value{
Backward probabilities in logarithm scale. In case of
multiple observations, these are computed independently for
each sequence.
}
\description{
Function \code{backwardProbs} computes the backward
probabilities of hidden Markov model in logarithm scale.
}

