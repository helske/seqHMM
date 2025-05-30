% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/posterior_probs.R
\name{posterior_probs}
\alias{posterior_probs}
\alias{posterior_probs.hmm}
\alias{posterior_probs.mhmm}
\alias{posterior_probs.nhmm}
\alias{posterior_probs.mnhmm}
\title{Posterior Probabilities for Hidden Markov Models}
\usage{
posterior_probs(model, ...)

\method{posterior_probs}{hmm}(model, ...)

\method{posterior_probs}{mhmm}(model, ...)

\method{posterior_probs}{nhmm}(model, ...)

\method{posterior_probs}{mnhmm}(model, ...)
}
\arguments{
\item{model}{A hidden Markov model object.}

\item{...}{Ignored.}
}
\value{
A data frame of posterior probabilities for each state and sequence.
}
\description{
Function \code{posterior_probs} computes the posterior probabilities of hidden
states of a (mixture) hidden Markov model.
}
\examples{
# Load a pre-defined MHMM
data("mhmm_biofam")

# Compute posterior probabilities
pb <- posterior_probs(mhmm_biofam)

}
