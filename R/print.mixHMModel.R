#' Print Method for a Mixture Hidden Markov Model
#'
#' Function \code{print.mixHMModel} prints the parameters of a mixture hidden Markov model.
#'
#'
#' @export
#' @param object Hidden Markov model of class \code{mixHMModel}.
#' @param ... Ignored.
#' @return Parameters of a hidden Markov model.
#' @seealso \code{\link{buildHMM}} and \code{\link{fitHMM}} for building and 
#'   fitting Hidden Markov models.


print.mixHMModel <- function(x){
  print(list(transitionMatrix=x$transitionMatrix, emissionMatrix=x$emissionMatrix, 
             initialProbs=x$initialProbs, beta=x$beta))
}