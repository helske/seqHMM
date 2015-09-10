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
  if(x$n_channels == 1){
    print.listof(list("Transition matrix" = round(x$transition_matrix, digits = 3), 
                      "Emission matrix" = round(x$emission_matrix, digits = 3), 
                      "Initial probabilities" = round(x$initial_probs, digits = 3)))
  }else{
    print.listof(list("Transition matrix" = round(x$transition_matrix, digits = 3)))
    cat("\n")
    cat("Emission matrix :\n\n")
    for(i in 1:length(x$emission_matrix)){
      cat(names(x$emission_matrix)[i], ":\n")
      print(round(x$emission_matrix[[i]], digits = 3))
    }
    cat("\n\n")
    print.listof(list("Initial probabilities" = round(x$initial_probs, digits = 3)))
  }
}