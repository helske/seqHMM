#' Print Method for a Hidden Markov Model
#'
#' Prints the parameters of a (mixture) hidden Markov model.
#'
#' @export
#' @rdname print
#' @method print hmm
#' @param x Hidden Markov model of class \code{hmm} or \code{mhmm}.
#' @param digits Minimum number of significant digits to print.
#' @param ... Further arguments to \code{print.default}.
#' @seealso \code{\link{build_hmm}} and \code{\link{fit_hmm}} for building and 
#'   fitting hidden Markov models and \code{\link{build_mhmm}} and 
#'   \code{\link{fit_mhmm}} for building and 
#'   fitting mixture hidden Markov models.
print.hmm <- function(x, digits = 3, ...){
  
  if (x$n_channels == 1) {
    print.listof(list(
      "Initial probabilities" = x$initial_probs,
      "Transition probabilities" = x$transition_probs, 
      "Emission probabilities" = x$emission_probs), digits = digits, ...)
  } else {
    print.listof(list("Initial probabilities" = x$initial_probs), digits = digits, ...)
    cat("\n")
    print.listof(list("Transition probabilities" = x$transition_probs), digits = digits, ...)
    cat("\n")
    cat("Emission probabilities :\n")
    print.listof(x$emission_probs, digits = digits, ...)
    cat("\n")
  }
}