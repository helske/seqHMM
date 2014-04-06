#' Forward Probabilities for Hidden Markov Model
#'
#' Function \code{forwardProbs} computes the forward probabilities of hidden Markov model in logarithm scale.
#'
#' @export 
#' @param model Hidden Markov model of class \code{HMModel}.
#' @return Forward probabilities in logarithm scale. In case of multiple observations,
#' these are computed independently for each sequence.
forwardProbs<-function(model){
  if(model$numberOfChannels==1){
    if(model$numberOfSequences==1){
      obs<-as.integer(factor(model$observations,labels=1:model$numberOfSymbols,levels=model$symbolNames))  
      miss<-is.na(obs)
      #if(miss[1]) stop("First observation cannot be missing!")
      storage.mode(miss)<-"integer"
      out<-.Fortran("forward",PACKAGE="seqHMM",NAOK = TRUE,
                    model$transitionMatrix,model$emissionMatrix,model$initialProbs,
                    obs,model$numberOfStates,
                    model$numberOfSymbols,model$lengthOfSequences,miss,
                    alpha=matrix(0,model$numberOfStates,model$lengthOfSequences))$alpha
    } else {
      obs<-apply(model$observations,2,factor,labels=1:model$numberOfSymbols,
                 levels=model$symbolNames)
      storage.mode(obs)<-"integer"
      miss<-is.na(obs)
      #if(any(miss[,1])) stop("First observation cannot be missing!")
      storage.mode(miss)<-"integer"
      out<-.Fortran("mvforward",PACKAGE="seqHMM",NAOK = TRUE,
                    object$transitionMatrix,object$emissionMatrix,object$initialProbs,
                    obs,model$numberOfStates,
                    model$numberOfSymbols,
                    model$lengthOfSequences,miss,model$numberOfSequences,
                    alpha=array(0,dim=c(model$numberOfStates,model$lengthOfSequences,model$numberOfSequences)))$alpha
    }
  } else {
    obsArray<-array(0,c(model$numberOfSequences,model$lengthOfSequences,model$numberOfChannels))
    maxNumberOfSymbols<-max(model$numberOfSymbols)
    emissionArray<-array(NA,c(model$numberOfStates,maxNumberOfSymbols,model$numberOfChannels))
    for(i in 1:model$numberOfChannels)
      emissionArray[,1:model$numberOfSymbols[i],i]<-model$emissionMatrix[[i]]
    
    if(model$numberOfSequences==1){
      for(i in 1:model$numberOfChannels)
        obsArray[,,i]<-as.integer(as.factor(model$observations[[i]]))
      miss<-is.na(obsArray)
      #if(any(miss[1,1,])) stop("First observation cannot be missing!")
      storage.mode(miss)<-"integer"
      out<-.Fortran("mcforward",PACKAGE="seqHMM",NAOK = TRUE,
                    model$transitionMatrix,emissionArray,model$initialProbs,
                    obsArray,model$numberOfStates,
                    maxNumberOfSymbols,model$lengthOfSequences,miss,
                    alpha=matrix(0,model$numberOfStates,model$lengthOfSequences),
                    model$numberOfChannels)$alpha
    } else {
      for(i in 1:model$numberOfChannels)
        obsArray[,,i]<-apply(model$observations[[i]],2,factor,labels=1:model$numberOfSymbols[[i]],
                             levels=model$symbolNames[[i]])
      storage.mode(obsArray)<-"integer"
      miss<-is.na(obsArray)
      #if(any(miss[,1,])) stop("First observation cannot be missing!")
      storage.mode(miss)<-"integer"
      out<-.Fortran("mcforward",PACKAGE="seqHMM",NAOK = TRUE,
                    model$transitionMatrix,emissionArray,model$initialProbs,
                    obsArray,model$numberOfStates,
                    maxNumberOfSymbols, model$lengthOfSequences,miss,model$numberOfSequences,
                    alpha=array(0,dim=c(model$numberOfStates,model$lengthOfSequences,model$numberOfSequences)),
                    model$numberOfChannels)$alpha
    }
    
  }
  out
}


#' Backward Probabilities for Hidden Markov Model
#'
#' Function \code{backwardProbs} computes the backward probabilities of hidden Markov model in logarithm scale.
#'
#' @export 
#' @param model Hidden Markov model of class \code{HMModel}.
#' @return Backward probabilities in logarithm scale. In case of multiple observations,
#' these are computed independently for each sequence.
backwardProbs<-function(model){
  if(model$numberOfChannels==1){
    if(model$numberOfSequences==1){
      obs<-as.integer(factor(model$observations,labels=1:model$numberOfSymbols,levels=model$symbolNames))  
      miss<-is.na(obs)
      #if(miss[1]) stop("First observation cannot be missing!")
      storage.mode(miss)<-"integer"
      out<-.Fortran("backward",PACKAGE="seqHMM",NAOK = TRUE,
                    model$transitionMatrix,model$emissionMatrix,model$initialProbs,
                    obs,model$numberOfStates,
                    model$numberOfSymbols,model$lengthOfSequences,miss,
                    beta=matrix(0,model$numberOfStates,model$lengthOfSequences))$beta
    } else {
      obs<-apply(model$observations,2,factor,labels=1:model$numberOfSymbols,
                 levels=model$symbolNames)
      storage.mode(obs)<-"integer"
      miss<-is.na(obs)
      #if(any(miss[,1])) stop("First observation cannot be missing!")
      storage.mode(miss)<-"integer"
      out<-.Fortran("mvbackward",PACKAGE="seqHMM",NAOK = TRUE,
                    object$transitionMatrix,object$emissionMatrix,object$initialProbs,
                    obs,model$numberOfStates,
                    model$numberOfSymbols,
                    model$lengthOfSequences,miss,model$numberOfSequences,
                    beta=array(0,dim=c(model$numberOfStates,model$lengthOfSequences,model$numberOfSequences)))$beta
    }
  } else {
    obsArray<-array(0,c(model$numberOfSequences,model$lengthOfSequences,model$numberOfChannels))
    maxNumberOfSymbols<-max(model$numberOfSymbols)
    emissionArray<-array(NA,c(model$numberOfStates,maxNumberOfSymbols,model$numberOfChannels))
    for(i in 1:model$numberOfChannels)
      emissionArray[,1:model$numberOfSymbols[i],i]<-model$emissionMatrix[[i]]
    
    if(model$numberOfSequences==1){
      for(i in 1:model$numberOfChannels)
        obsArray[,,i]<-as.integer(as.factor(model$observations[[i]]))
      miss<-is.na(obsArray)
      #if(any(miss[1,1,])) stop("First observation cannot be missing!")
      storage.mode(miss)<-"integer"
      out<-.Fortran("mcbackward",PACKAGE="seqHMM",NAOK = TRUE,
                    model$transitionMatrix,emissionArray,model$initialProbs,
                    obsArray,model$numberOfStates,
                    maxNumberOfSymbols,model$lengthOfSequences,miss,
                    beta=matrix(0,model$numberOfStates,model$lengthOfSequences),
                    model$numberOfChannels)$beta
    } else {
      for(i in 1:model$numberOfChannels)
        obsArray[,,i]<-apply(model$observations[[i]],2,factor,labels=1:model$numberOfSymbols[[i]],
                             levels=model$symbolNames[[i]])
      storage.mode(obsArray)<-"integer"
      miss<-is.na(obsArray)
      #if(any(miss[,1,])) stop("First observation cannot be missing!")
      storage.mode(miss)<-"integer"
      out<-.Fortran("mcbackward",PACKAGE="seqHMM",NAOK = TRUE,
                    model$transitionMatrix,emissionArray,model$initialProbs,
                    obsArray,model$numberOfStates,
                    maxNumberOfSymbols, model$lengthOfSequences,miss,model$numberOfSequences,
                    beta=array(0,dim=c(model$numberOfStates,model$lengthOfSequences,model$numberOfSequences)),
                    model$numberOfChannels)$beta
    }
    
  }
  out
  
}
