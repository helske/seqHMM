#' @export
BIC.HMModel<-function(object,expand=FALSE, trim=FALSE,maxit=0,zerotol=1e-8,lltol=.Machine$double.eps ^ 0.5,...){
  
  if(expand)
    object<-MCtoSC(object)
  
  if(trim){
    trimmed<-trimHMM(object,maxit=maxit,return.loglik=TRUE,zerotol=zerotol,...)
    object<-trimmed$model
    loglik<-trimmed$loglik
  } else loglik<-logLik(object)
  
  if(object$numberOfChannels==1){
    b <- -2*loglik + log(sum(!is.na(object$observations)))*
      (sum(object$initialProbs>0)+sum(object$transitionMatrix>0)+sum(object$emissionMatrix>0))
  } else {
    
    b <- -2*loglik + log(sum(!sapply(object$observations,is.na)))*
      (sum(object$initialProbs>0)+sum(object$transitionMatrix>0)+
         sum(sapply(object$emissionMatrix,function(x) sum(x>0))))
  }
  
  b
}