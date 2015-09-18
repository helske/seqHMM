#' @export
#' @method print mhmm
#' @rdname print
print.mhmm <- function(x, ...){
  
  cat("Coefficients :\n\n")
  print(x$coefficients, digits = 3)
  
  cat("\nInitial probabilities :\n\n")
  print.listof(x$initial_probs, digits = 3)
  
  cat("Transition matrix :\n\n")
  print.listof(x$transition_matrix, digits = 3)
  
  cat("Emission matrix :\n\n")
  if (x$n_channels == 1) {
    print.listof(x$emission_matrix, digits = 3)
  } else {
    for(i in 1:length(x$emission_matrix)){
      cat(names(x$emission_matrix)[i], ":\n\n")
      print.listof(x$emission_matrix[[i]], digits = 3)
    }
  }
  cat("\n")
}