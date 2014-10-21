#' Build a Hidden Markov Model
#' 
#' Function buildHMM constructs an object of class \code{HMModel}.
#' 
#'
#' The set of unique observed symbols and it's ordering (which is used in emissionMatrix) is taken 
#' from internal representation of 
#' factor levels after removing unused levels. More precisely following code is used: 
#' \code{intersect(levels(observations[,1]),as.character(unlist(observations))))}.
#' 
#' @import TraMineR
#' @importFrom Rcpp evalCpp
#' @export
#' @useDynLib seqHMM
#' @param observations Data frame or TraMineR stslist containing observations where each row corresponds to one sequence and 
#' each column one time point which must be a factor.
#' @param transitionMatrix Matrix of transition probabilities.
#' @param emissionMatrix Matrix of emission probabilities.
#' @param initialProbs Vector of initial state probabilities.
#' @param X NOT YET IMPLEMENTED. Time-constant covariates for initial state probabilities. 
#' Must be a matrix with dimensions p x k where p is the number of sequences and k is the number of covariates.
#' @param stateNames Optional labels for the observations.
#' @return Object of class \code{HMModel}
#' 
buildHMM<-function(observations,transitionMatrix,emissionMatrix,initialProbs,X,stateNames=NULL,symbolNames=NULL){

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
  
  numberOfSequences<-nrow(observations)
  lengthOfSequences<-ncol(observations)
  symbolNames<-intersect(levels(observations[,1]),as.character(unlist(observations)))
  numberOfSymbols<-length(symbolNames)
  
  if(numberOfStates!=dim(emissionMatrix)[1])
    stop("Number of rows in emissionMatrix is not equal to the number of states.")
  if(numberOfSymbols!=dim(emissionMatrix)[2])
    stop("Number of columns in emissionMatrix is not equal to the number of symbols.")
  if(!isTRUE(all.equal(rep(1,numberOfStates),rowSums(emissionMatrix),check.attributes=FALSE)))
    stop("Emission probabilities in emissionMatrix do not sum to one.")
  dimnames(emissionMatrix)<-list(stateNames=stateNames,symbolNames=symbolNames)
  
  if(!missing(X)){
    warning("Covariates not yet implemented.")
    if(nrow(X)!=p)
      stop("Wrong dimensions in X.")
    numberOfCovariates<-ncol(X)
  }
  
  model<-list(observations=observations,transitionMatrix=transitionMatrix,
              emissionMatrix=emissionMatrix,initialProbs=initialProbs,stateNames=stateNames,
              symbolNames=symbolNames,lengthOfSequences=lengthOfSequences,
              numberOfSequences=numberOfSequences,numberOfSymbols=numberOfSymbols,
              numberOfStates=numberOfStates)
  class(model)<-"HMModel"
  model
}