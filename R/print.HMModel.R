#' Print Method for a Hidden Markov Model
#'
#' Function \code{print.HMModel} prints the parameters of a hidden Markov model.
#'
#'
#' @export
#' @param object Hidden Markov model of class \code{HMModel}.
#' @param ... Ignored.
#' @return Parameters of a hidden Markov model.
#' @seealso \code{\link{buildHMM}} and \code{\link{fitHMM}} for building and 
#'   fitting Hidden Markov models.


print.HMModel <- function(x){
  print(list(transitionMatrix=x$transitionMatrix, emissionMatrix=x$emissionMatrix, 
       initialProbs=x$initialProbs))
}