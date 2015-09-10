#' @export
#' @method print mhmm
#' @rdname print
print.mhmm <- function(x, ...){
  if(x$n_channels == 1){
    cat("Transition matrix :\n\n")
    print.listof(lapply(x$transition_matrix, round, digits = 3))
    cat("Emission matrix :\n\n")
    print.listof(lapply(x$emission_matrix, round, digits = 3))
    cat("Initial probabilities :\n\n")
    print.listof(lapply(x$initial_probs, round, digits = 3))
    cat("\n\n\n")
  }else{
    cat("Transition matrix :\n\n")
    print.listof(lapply(x$transition_matrix, round, digits = 3))
    cat("Emission matrix :\n\n")
    for(i in 1:length(x$emission_matrix)){
      cat(names(x$emission_matrix)[i], ":\n\n")
      print.listof(lapply(x$emission_matrix[[i]], round, digits = 3))
    }
    cat("Initial probabilities :\n\n")
    print.listof(lapply(x$initial_probs, round, digits = 3))
    cat("\n\n\n")
  }
}