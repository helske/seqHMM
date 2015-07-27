#' Print Method for a Hidden Markov Model
#'
#' Function \code{print.HMModel} prints the parameters of a hidden Markov model.
#'
#'
#' @export
#' @rdname print
#' @method print HMModel
#' @param x Hidden Markov model of class \code{HMModel} or \code{mixHMModel}.
#' @param ... Ignored.
#' @seealso \code{\link{buildHMM}} and \code{\link{fitHMM}} for building and 
#'   fitting hidden Markov models and \code{\link{buildMixHMM}} and 
#'   \code{\link{fitMixHMM}} for building and 
#'   fitting mixture hidden Markov models.


print.HMModel <- function(x, ...){
  print(list(transitionMatrix=x$transitionMatrix, emissionMatrix=x$emissionMatrix, 
       initialProbs=x$initialProbs))
}