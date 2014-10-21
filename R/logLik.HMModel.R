#' Log-likelihood of the Hidden Markov Model
#'
#' Function \code{logLik.HMModel} computes the log-likelihood value of a hidden Markov model.
#'
#'
#' @export
#' @aliases logLik logLik.HMModel
#' @param object Hidden Markov model of class \code{HMModel}.
#' @return Log-likelihood of hidden Markov model.
logLik.HMModel<-function(object,...){
  
  obsArray<-data.matrix(object$observations)-1
  obsArray[is.na(obsArray)]<-object$numberOfSymbols
  storage.mode(obsArray)<-"integer"
  
  logLikHMM(object$transitionMatrix, cbind(object$emissionMatrix,1), 
            object$initialProbs, obsArray)
}