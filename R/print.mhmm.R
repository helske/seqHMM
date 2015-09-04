#' @export
#' @rdname print
#' @method print mhmm

print.mhmm <- function(x, ...){
  print(list(transition_matrix=x$transition_matrix, emission_matrix=x$emission_matrix, 
             initial_probs=x$initial_probs, beta=x$beta))
}