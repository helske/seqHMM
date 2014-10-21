#' Build a Multichannel Hidden Markov Model
#' 
#' Function buildMCHMM constructs an object of class \code{MCHMModel}.
#' 
#' The set of unique observed symbols and it's ordering (which is used in emissionMatrix) is taken 
#' from internal representation of 
#' factor levels after removing unused levels. More precisely following code is used for each channel: 
#' \code{intersect(levels(observations[,1]),as.character(unlist(observations))))}.
#' 
#' 
#' @export
#' @useDynLib seqHMM
#' @param observations Observations. Must be a list of data frames where each data frame corresponds to one channel.
#' Each row corresponds to one sequence and each column one time point which must be a factor.
#' @param transitionMatrix Matrix of transition probabilities.
#' @param emissionMatrix List of emission probability matrices.
#' @param initialProbs Vector of initial state probabilities.
#' @param X NOT YET IMPLEMENTED. Time-constant covariates for initial state probabilities. 
#' Must be a matrix with dimensions p x k where p is the number of sequences and k is the number of covariates.
#' @param stateNames Optional labels for the observations.
#' @return Object of class \code{HMModel}
#' 
buildMCHMM<-function(observations,transitionMatrix,emissionMatrix,initialProbs,X,
                  stateNames=NULL){
  
  # determine number of states
  numberOfStates<-length(initialProbs)
  if(is.null(stateNames))
    stateNames<-as.character(1:numberOfStates)
  channelNames<-NULL
  
  if(dim(transitionMatrix)[1]!=dim(transitionMatrix)[2] ||
       numberOfStates!=dim(transitionMatrix)[1])
    stop("Dimensions of transitionMatrix and length of initialProbs do not match.")
  
  if(!isTRUE(all.equal(rowSums(transitionMatrix),rep(1,dim(transitionMatrix)[1]),check.attributes=FALSE)))
    stop("Transition probabilities in transitionMatrix do not sum to one.")
  if(!isTRUE(all.equal(sum(initialProbs),1,check.attributes=FALSE)))
    stop("Initial state probabilities do not sum to one.")
  dimnames(transitionMatrix)<-list(from=stateNames,to=stateNames)
  
  numberOfChannels<-length(observations)
  
  numberOfSequences<-nrow(observations[[1]])
  lengthOfSequences<-ncol(observations[[1]])

  symbolNames<-lapply(observations,function(x) intersect(levels(x[,1]),as.character(unlist(x))))
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
  
  
  if(!missing(X)){
    warning("Covariates not yet implemented.")
    if(nrow(X)!=p)
      stop("Wrong dimensions in X.")
    numberOfCovariates<-ncol(X)
  }
  
  model<-list(observations=observations,transitionMatrix=transitionMatrix,
              emissionMatrix=emissionMatrix,initialProbs=initialProbs,stateNames=stateNames,
              symbolNames=symbolNames,channelNames=channelNames,lengthOfSequences=lengthOfSequences,
              numberOfSequences=numberOfSequences,
              numberOfSymbols=numberOfSymbols,numberOfStates=numberOfStates,
              numberOfChannels=numberOfChannels)
  class(model)<-"MCHMModel"
  model
}