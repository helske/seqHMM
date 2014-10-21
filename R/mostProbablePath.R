#' Most Probable Path of Hidden States of Hidden Markov Model given the sequence.
#'
#' Function \code{mostProbablePath} computes the most probable path of the hidden states of the 
#' hidden Markov model given the observed sequence.
#'
#' @export 
#' @param model Hidden Markov model of class \code{HMModel}.
#' @return List which contains the most probable path of states given the observations and its log-probability. In case of multiple observations,
#' most probable path is computed independently for each sequence.
mostProbablePath<-function(model){
  
  model$transitionMatrix <- log(model$transitionMatrix)
  model$transitionMatrix[!is.finite(model$transitionMatrix)]<- -0.1*.Machine$double.xmax
  if(is.list(model$emissionMatrix)){
    model$emissionMatrix<-lapply(model$emissionMatrix,function(x){
      x<-log(x)
      x[!is.finite(x)]<- -0.1*.Machine$double.xmax
      x
    })
    
  } else {
    model$emissionMatrix <- log(model$emissionMatrix)
    model$emissionMatrix[!is.finite(model$emissionMatrix)]<- -0.1*.Machine$double.xmax
  }
  
  model$initialProbs <- log(model$initialProbs)
  model$initialProbs[!is.finite(model$initialProbs)]<- -0.1*.Machine$double.xmax
  
  if(model$numberOfChannels==1){
    
    if(model$numberOfSequences==1){
      obs<-as.integer(as.factor(unlist(model$observations)))
      miss<-is.na(obs)
      #if(miss[1]) stop("First observation cannot be missing!")
      storage.mode(miss)<-"integer"
      out<-.Fortran("viterbi",PACKAGE="seqHMM",NAOK = TRUE,
                    model$transitionMatrix,model$emissionMatrix,model$initialProbs,
                    obs,model$numberOfStates,
                    model$numberOfSymbols,model$lengthOfSequences,miss
                    ,q=integer(model$lengthOfSequences),logp=0)
      mpp<-rownames(model$transitionMatrix)[out$q]
    } else {
      obs<-apply(model$observations,2,factor,labels=1:model$numberOfSymbols,
                 levels=model$symbolNames)
      miss<-is.na(obs)
      q<-matrix(0,model$numberOfSequences,model$lengthOfSequences)
      storage.mode(q)<-storage.mode(obs)<-storage.mode(miss)<-"integer"
      out<-.Fortran("mvviterbi",PACKAGE="seqHMM",NAOK = TRUE,
                    model$transitionMatrix,model$emissionMatrix,model$initialProbs,
                    obs,model$numberOfStates,
                    model$numberOfSymbols,
                    model$lengthOfSequences,miss,model$numberOfSequences,
                    q=q,logp=numeric(model$numberOfSequences))
      mpp<-apply(out$q,2,function(x) rownames(model$transitionMatrix)[x])
    }
  } else {
    maxNumberOfSymbols<-max(model$numberOfSymbols)
    emissionArray<-array(NA,c(model$numberOfStates,maxNumberOfSymbols,model$numberOfChannels))
    for(i in 1:model$numberOfChannels)
      emissionArray[,1:model$numberOfSymbols[i],i]<-model$emissionMatrix[[i]]
    
    if(model$numberOfSequences==1){
      obsArray<-array(0,c(model$lengthOfSequences,model$numberOfChannels))
      
      for(i in 1:model$numberOfChannels)
        obsArray[,i]<-as.integer(as.factor(unlist(model$observations[[i]])))
      miss<-is.na(obsArray)
      #if(any(miss[1,1,])) stop("First observation cannot be missing!")
      storage.mode(miss)<-storage.mode(obsArray)<-"integer"
      out<-.Fortran("mcviterbi",PACKAGE="seqHMM",NAOK = TRUE,
                    model$transitionMatrix,emissionArray,model$initialProbs,
                    obsArray,model$numberOfStates,
                    maxNumberOfSymbols,model$lengthOfSequences,miss
                    ,q=integer(model$lengthOfSequences),model$numberOfChannels,logp=logp)
      mpp<-rownames(model$transitionMatrix)[out$q]
    } else {
      obsArray<-array(0,c(model$numberOfSequences,model$lengthOfSequences,model$numberOfChannels))
      
      for(i in 1:model$numberOfChannels)
        obsArray[,,i]<-apply(model$observations[[i]],2,factor,labels=1:model$numberOfSymbols[[i]],
                             levels=model$symbolNames[[i]])
      miss<-is.na(obsArray)
      q<-matrix(0,model$numberOfSequences,model$lengthOfSequences)
      storage.mode(q)<-storage.mode(miss)<-storage.mode(obsArray)<-"integer"
      out<-.Fortran("mvmcviterbi",PACKAGE="seqHMM",NAOK = TRUE,
                    model$transitionMatrix,emissionArray,model$initialProbs,
                    obsArray,model$numberOfStates,
                    maxNumberOfSymbols,
                    model$lengthOfSequences,miss,model$numberOfSequences,
                    q=q,model$numberOfChannels,logp=numeric(model$numberOfSequences))
      mpp<-apply(out$q,2,function(x) rownames(model$transitionMatrix)[x])
    }
    
  }
  list(mpp=mpp,logP=out$logp)
}
