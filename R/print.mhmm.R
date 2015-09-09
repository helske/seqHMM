#' @export
#' @method print mhmm
#' @rdname print
print.mhmm <- function(x, ...){
  print.listof(list("transition matrix" = x$transition_matrix, 
    "emission matrix" = x$emission_matrix, 
    "initial probabilities" = x$initial_probs, 
    "beta" = x$beta))
}