#' Print Method for a Hidden Markov Model
#'
#' Prints the parameters of a (mixture) hidden Markov model.
#'
#' @export
#' @rdname print
#' @method print hmm
#' @param x Hidden Markov model of class \code{hmm} or \code{mhmm}.
#' @param ... Ignored.
#' @seealso \code{\link{build_hmm}} and \code{\link{fit_hmm}} for building and 
#'   fitting hidden Markov models and \code{\link{build_mhmm}} and 
#'   \code{\link{fit_mhmm}} for building and 
#'   fitting mixture hidden Markov models.
print.hmm <- function(x, ...){
  print.listof(list("transition matrix" = x$transition_matrix, 
    "emission matrix" = x$emission_matrix, 
    "initial probabilities" = x$initial_probs))
}