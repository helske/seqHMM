#' @export
#' @rdname print
#' @method print HMModel

print.mixHMModel <- function(x, ...){
  print(list(transitionMatrix=x$transitionMatrix, emissionMatrix=x$emissionMatrix, 
             initialProbs=x$initialProbs, beta=x$beta))
}