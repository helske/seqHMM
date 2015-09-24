#' @export
#' @method print mhmm
#' @rdname print
print.mhmm <- function(x, digits = 3, ...){
  
  cat("Coefficients :\n")
  print(x$coefficients, digits = digits, ...)
  
  cat("\nInitial probabilities :\n")
  print.listof(x$initial_probs, digits = digits, ...)
  
  cat("Transition probabilities :\n")
  print.listof(x$transition_matrix, digits = digits, ...)
  
  cat("Emission probabilities :\n")
  if (x$n_channels == 1) {
    print.listof(x$emission_matrix, digits = digits, ...)
  } else {
    for(i in 1:length(x$emission_matrix)){
      cat(names(x$emission_matrix)[i], ":\n\n")
      print.listof(x$emission_matrix[[i]], digits = digits, ...)
    }
  }
  cat("\n")
}