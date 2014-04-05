#' Posterior Probabilities for Hidden Markov Model
#'
#' Function \code{posteriorProbs} computes the posterior probabilities of hidden states of
#' hidden Markov model given the observations in logarithm scale.
#'
#' @export 
#' @param model Hidden Markov model of class \code{HMModel}.
#' @return Posterior probabilities in logarithm scale. In case of multiple observations,
#' these are computed independently for each sequence.
posteriorProbs<-function(model){
  if(model$numberOfChannels==1){
    if(model$numberOfSequences==1){
      obs<-as.integer(factor(model$observations,labels=1:model$numberOfSymbols,levels=model$symbolNames))  
      miss<-is.na(obs)
      #if(miss[1]) stop("First observation cannot be missing!")
      storage.mode(miss)<-"integer"
      out<-.Fortran("posterior",PACKAGE="LifeSequenceHMM",NAOK = TRUE,
                    model$transitionMatrix,model$emissionMatrix,model$initialProbs,
                    obs,model$numberOfStates,
                    model$numberOfSymbols,model$lengthOfSequences,miss,
                    posterior=matrix(0,model$numberOfStates,model$lengthOfSequences))$posterior
    } else {
      obs<-apply(model$observations,2,factor,labels=1:model$numberOfSymbols,
                 levels=model$symbolNames)
      storage.mode(obs)<-"integer"
      miss<-is.na(obs)
      #if(any(miss[,1])) stop("First observation cannot be missing!")
      storage.mode(miss)<-"integer"
      out<-.Fortran("mvposterior",PACKAGE="LifeSequenceHMM",NAOK = TRUE,
                    object$transitionMatrix,object$emissionMatrix,object$initialProbs,
                    obs,model$numberOfStates,
                    model$numberOfSymbols,
                    model$numberOfSequences,miss,model$lengthOfSequences,
                    posterior=array(0,dim=c(model$numberOfStates,model$lengthOfSequences,model$numberOfSequences)))$posterior
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
      out<-.Fortran("mcposterior",PACKAGE="LifeSequenceHMM",NAOK = TRUE,
                    model$transitionMatrix,emissionArray,model$initialProbs,
                    obsArray,model$numberOfStates,
                    maxNumberOfSymbols,model$lengthOfSequences,miss,
                    posterior=matrix(0,model$numberOfStates,model$lengthOfSequences),
                    model$numberOfChannels)$posterior
    } else {
      for(i in 1:model$numberOfChannels)
        obsArray[,,i]<-apply(model$observations[[i]],2,factor,labels=1:model$numberOfSymbols[[i]],
                             levels=model$symbolNames[[i]])
      storage.mode(obsArray)<-"integer"
      miss<-is.na(obsArray)
      #if(any(miss[,1,])) stop("First observation cannot be missing!")
      storage.mode(miss)<-"integer"
      out<-.Fortran("mcposterior",PACKAGE="LifeSequenceHMM",NAOK = TRUE,
                    model$transitionMatrix,emissionArray,model$initialProbs,
                    obsArray,model$numberOfStates,
                    maxNumberOfSymbols,model$lengthOfSequences,miss,model$numberOfSequences,
                    posterior=array(0,dim=c(model$numberOfStates,model$lengthOfSequences,model$numberOfSequences)),
                    model$numberOfChannels)$posterior
    }
    
  }
  out
}
