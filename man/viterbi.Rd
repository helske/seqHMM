\name{viterbi}
\alias{viterbi}
\title{Viterbi Algorithm for Hidden Markov Model}
\usage{
viterbi(model)
}
\arguments{
  \item{model}{Hidden Markov model of class
  \code{HMModel}.}
}
\value{
Most probable path of states given the observations. In
case of multiple sequences, most probable path is computed
independently for each sequence.
}
\description{
Function \code{viterbi} computes the most probable path of
the hidden states of the hidden Markov model.
}

