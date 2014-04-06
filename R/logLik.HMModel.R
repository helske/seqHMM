#' Log-likelihood of the Hidden Markov Model
#'
#' Function \code{logLik.HMModel} computes the log-likelihood value of a hidden Markov model.
#'
#'
#' @export
#' @S3method logLik HMModel
#' @method logLik HMModel
#' @aliases logLik logLik.HMModel
#' @param object Hidden Markov model of class \code{HMModel}.
#' @return Log-likelihood of hidden Markov model.
logLik.HMModel<-function(object,...){
  if(object$numberOfChannels==1){
    if(object$numberOfSequences==1){
      obs<-as.integer(factor(object$observations,labels=1:object$numberOfSymbols,levels=object$symbolNames))
      miss<-is.na(obs)
      storage.mode(miss)<-"integer"
      ll<-.Fortran("hmmloglik",PACKAGE="seqHMM",NAOK = TRUE,
                   object$transitionMatrix,object$emissionMatrix,object$initialProbs,
                   obs,object$numberOfStates,object$numberOfSymbols,object$lengthOfSequences,
                   miss,logLik=double(1))$logLik
      
    } else {
      obs<-apply(object$observations,2,factor,labels=1:dim(object$emissionMatrix)[2],
                 levels=object$symbolNames)
      storage.mode(obs)<-"integer"
      miss<-is.na(obs)
      storage.mode(miss)<-"integer"
      ll<-.Fortran("mvhmmloglik",PACKAGE="seqHMM",NAOK = TRUE,
                   object$transitionMatrix,object$emissionMatrix,object$initialProbs,
                   obs,object$numberOfStates,object$numberOfSymbols,object$lengthOfSequences,
                   miss,object$numberOfSequences,logLik=double(1))$logLik
      
    } 
  } else {
    obsArray<-array(0,c(object$numberOfSequences,object$lengthOfSequences,object$numberOfChannels))
    maxNumberOfSymbols<-max(object$numberOfSymbols)
    emissionArray<-array(NA,c(object$numberOfStates,maxNumberOfSymbols,object$numberOfChannels))
    for(i in 1:object$numberOfChannels)
      emissionArray[,1:object$numberOfSymbols[i],i]<-object$emissionMatrix[[i]]
    
    if(object$numberOfSequences==1){
      for(i in 1:object$numberOfChannels)
        obsArray[,,i]<-as.integer(as.factor(object$observations[[i]]))
      miss<-is.na(obsArray)
      storage.mode(miss)<-"integer"
      ll<-.Fortran("mchmmloglik",PACKAGE="seqHMM",NAOK = TRUE,
                   object$transitionMatrix,emissionArray,object$initialProbs,
                   obsArray,object$numberOfStates,maxNumberOfSymbols,
                   object$lengthOfSequences,miss,logLik=double(1),object$numberOfChannels)$logLik
    } else {
      for(i in 1:object$numberOfChannels)
        obsArray[,,i]<-apply(object$observations[[i]],2,factor,labels=1:object$numberOfSymbols[[i]],
                             levels=object$symbolNames[[i]])
      storage.mode(obsArray)<-"integer"
      miss<-is.na(obsArray)
      storage.mode(miss)<-"integer"
      
      ll<-.Fortran("mvmchmmloglik",PACKAGE="seqHMM",NAOK = TRUE,
                   object$transitionMatrix,emissionArray,object$initialProbs,
                   obsArray,object$numberOfStates,maxNumberOfSymbols,
                   object$lengthOfSequences,miss,object$numberOfSequences,logLik=double(1),
                   object$numberOfChannels)$logLik
    }
  }    
  ll
}