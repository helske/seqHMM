#' Build a Hidden Markov Model
#' 
#' Function buildHMM constructs an object of class \code{HMModel}.
#' 
#' @import TraMineR
#' @importFrom Rcpp evalCpp
#' @export
#' @useDynLib seqHMM
#' @param observations TraMineR stslist containing the sequences, or a list of such objects (one for each channel).
#' @param transitionMatrix Matrix of transition probabilities.
#' @param emissionMatrix Matrix of emission probabilities or a list of such objects (one for each channel).
#' @param initialProbs Vector of initial state probabilities.
#' @param X NOT YET IMPLEMENTED. Time-constant covariates for initial state probabilities. 
#' Must be a matrix with dimensions p x k where p is the number of sequences and k is the number of covariates.
#' @param stateNames Optional labels for the hidden states
#' @return Object of class \code{HMModel}
#' 
buildHMM<-function(observations,transitionMatrix,emissionMatrix,initialProbs,X,stateNames=NULL){
  
  # determine number of states
  numberOfStates<-length(initialProbs)
  if(is.null(stateNames))
    stateNames<-as.character(1:numberOfStates)
  
  if(dim(transitionMatrix)[1]!=dim(transitionMatrix)[2] ||
       numberOfStates!=dim(transitionMatrix)[1])
    stop("Dimensions of transitionMatrix and length of initialProbs do not match.")
  
  if(!isTRUE(all.equal(rowSums(transitionMatrix),rep(1,dim(transitionMatrix)[1]),check.attributes=FALSE)))
    stop("Transition probabilities in transitionMatrix do not sum to one.")
  if(!isTRUE(all.equal(sum(initialProbs),1,check.attributes=FALSE)))
    stop("Initial state probabilities do not sum to one.")
  dimnames(transitionMatrix)<-list(from=stateNames,to=stateNames)
  
  if(is.list(emissionMatrix) && length(emissionMatrix)==1){
    emissionMatrix <- emissionMatrix[[1]]
    observations <- observations[[1]]
  }
  
  if(is.list(emissionMatrix)){
    if(length(observations)!=length(emissionMatrix)){
      stop("Number of channels defined by emissionMatrix differs from one defined by observations.")
    }
    numberOfChannels <- length(emissionMatrix)
    
    numberOfSequences<-nrow(observations[[1]])
    lengthOfSequences<-ncol(observations[[1]])
    
    symbolNames<-lapply(observations,alphabet)
    numberOfSymbols<-sapply(symbolNames,length)
    
    if(!isTRUE(all.equal(rep(numberOfStates,numberOfChannels),sapply(emissionMatrix,nrow),check.attributes=FALSE)))
      stop("Number of rows in emissionMatrix is not equal to the number of states.")
    if(!isTRUE(all.equal(numberOfSymbols,sapply(emissionMatrix,ncol),check.attributes=FALSE)))
      stop("Number of columns in emissionMatrix is not equal to the number of symbols.")
    if(!isTRUE(all.equal(c(sapply(emissionMatrix,rowSums)),rep(1,numberOfChannels*numberOfStates),check.attributes=FALSE)))
      stop("Emission probabilities in emissionMatrix do not sum to one.")
    
    channelNames<-names(observations)  
    if(is.null(channelNames))
      channelNames<-as.character(1:numberOfChannels)
    for(i in 1:numberOfChannels)
      dimnames(emissionMatrix[[i]])<-list(stateNames=stateNames,symbolNames=symbolNames[[i]])
    names(emissionMatrix)<-channelNames
  } else {
    numberOfChannels <- 1
    channelNames<-NULL
    numberOfSequences<-nrow(observations)
    lengthOfSequences<-ncol(observations)
    symbolNames<-alphabet(observations)
    numberOfSymbols<-length(symbolNames)
    
    if(numberOfStates!=dim(emissionMatrix)[1])
      stop("Number of rows in emissionMatrix is not equal to the number of states.")
    if(numberOfSymbols!=dim(emissionMatrix)[2])
      stop("Number of columns in emissionMatrix is not equal to the number of symbols.")
    if(!isTRUE(all.equal(rep(1,numberOfStates),rowSums(emissionMatrix),check.attributes=FALSE)))
      stop("Emission probabilities in emissionMatrix do not sum to one.")
    dimnames(emissionMatrix)<-list(stateNames=stateNames,symbolNames=symbolNames)
    
  }
  
  if(!missing(X)){
    warning("Covariates not yet implemented.")
    if(nrow(X)!=numberOfSequences)
      stop("Wrong dimensions in X.")
    numberOfCovariates<-ncol(X)
  }
  
  model<-list(observations=observations,transitionMatrix=transitionMatrix,
              emissionMatrix=emissionMatrix,initialProbs=initialProbs,stateNames=stateNames,
              symbolNames=symbolNames,channelNames=channelNames,lengthOfSequences=lengthOfSequences,
              numberOfSequences=numberOfSequences,
              numberOfSymbols=numberOfSymbols,numberOfStates=numberOfStates,
              numberOfChannels=numberOfChannels)
  class(model)<-"HMModel"
  model
}