

#' Compute BIC for Hidden Markov Model
#' 
#' @import stats
#' @export
#' @rdname BIC
#' @param object Object of class \code{HMModel} or \code{MCHMModel}.
#' @param trim Should trimming be performed before computing BIC? See \code{trimHMM} for details. Default is FALSE.
#' @param ... further parameters to \code{trimHMM}.
BIC.HMModel<-function(object,trim=FALSE,...){
  
  #if(expand)
  #  object<-MCtoSC(object)
  
  if(trim){
    trimmed<-trimHMM(object,return.loglik=TRUE,...)
    object<-trimmed$model
    loglik<-trimmed$loglik
  } else loglik<-logLik(object)
  
  
  if(object$numberOfChannels==1){
    -2*loglik + log(sum(!is.na(object$observations)))*
      (sum(object$initialProbs>0)+sum(object$transitionMatrix>0)+
         sum(object$emissionMatrix>0))
  } else{  
    -2*loglik + log(sum(!sapply(object$observations,is.na)))*
      (sum(object$initialProbs>0)+sum(object$transitionMatrix>0)+
         sum(sapply(object$emissionMatrix,function(x) sum(x>0))))
  }
}