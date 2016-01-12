#' Reorganize a mixture hidden Markov model to a list of separate hidden Markov models
#' (covariates ignored)
#'
#' The \code{separate_mhmm} function reorganizes the parameters of a \code{mhmm} object
#' into a list where each list component is an object of class \code{hmm} consisting of the
#' parameters of the corresponding cluster.
#'
#' @export
#' @param model Mixture hidden Markov model of class \code{mhmm}.
#'
#' @return List with components of class \code{hmm}.
#'
#' @seealso \code{\link{build_mhmm}} and \code{\link{fit_model}}
#' for building and fitting MHMMs; and \code{\link{mhmm_biofam}} for
#' more information on the model used in examples.
#'
#' @examples
#' # Loading mixture hidden Markov model (mhmm object)
#' # of the biofam data
#' data("mhmm_biofam")
#'
#' # Separate models for clusters
#' sep_hmm <- separate_mhmm(mhmm_biofam)
#' 
#' # Plotting the model for the first cluster
#' plot(sep_hmm[[1]])


separate_mhmm <- function(model){

  divmodels <- replicate(model$n_clusters, list())

  for(i in 1:model$n_clusters){
    divmodels[[i]] <- build_hmm(observations=model$observations,
                                transition_probs=model$transition_probs[[i]],
                                emission_probs=model$emission_probs[[i]],
                                initial_probs=model$initial_probs[[i]])
  }
  divmodels
}
