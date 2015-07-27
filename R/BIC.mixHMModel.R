#' @export
#' @rdname BIC
BIC.mixHMModel<-function(object,trim=FALSE,...){

  
  if(trim){
    trimmed<-trimHMM(object,return.loglik=TRUE,...)
    object<-trimmed$model
    loglik<-trimmed$loglik
  } else loglik<-logLik(object)
  
  
  -2*loglik + log(object$numberOfSequences*object$lengthOfSequences)*
    (sum(unlist(object$initialProbs)>0)+sum(unlist(object$transitionMatrix)>0)+
       sum(unlist(object$emissionMatrix)>0))
}