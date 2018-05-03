#' @export
#' @method print mhmm
#' @rdname print
print.mhmm <- function(x, digits = 3, ...){
  
  cat("Coefficients :\n")
  print(x$coefficients, digits = digits, ...)
  
  if (attr(x, "type") != "lcm") {
  cat("\nInitial probabilities :\n")
  print.listof(x$initial_probs, digits = digits, ...)
  
  cat("Transition probabilities :\n")
  print.listof(x$transition_probs, digits = digits, ...)
  } else cat("\n")
  
  if (attr(x, "type") != "mmm") {
    cat("Emission probabilities :\n")
    
    if (x$n_channels == 1) {
      if (attr(x, "type") == "lcm") { 
        tmp <- do.call(rbind, x$emission_probs)
        rownames(tmp) <- x$state_names
        print(tmp, digits = digits, ...)
      } else print.listof(x$emission_probs, digits = digits, ...)
    } else {
      if (attr(x, "type") == "lcm") {
        for (i in 1:x$n_channels) {
          cat(x$channel_names[i], ":\n")
          print(do.call(rbind, sapply(x$emission_probs,"[",i)), digits = digits, ...)
          cat("\n")
        }
      } else {
        for(i in 1:x$n_clusters){
          cat(names(x$emission_probs)[i], ":\n\n")
          print.listof(x$emission_probs[[i]], digits = digits, ...)
        }
      }
    }
  }
  
  cat("\n")
}