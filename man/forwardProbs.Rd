\name{forwardProbs}
\alias{forwardProbs}
\title{Forward Probabilities for Hidden Markov Model}
\usage{
forwardProbs(model)
}
\arguments{
  \item{model}{Hidden Markov model of class
  \code{HMModel}.}
}
\value{
Forward probabilities in logarithm scale. In case of
multiple observations, these are computed independently for
each sequence.
}
\description{
Function \code{forwardProbs} computes the forward
probabilities of hidden Markov model in logarithm scale.
}

