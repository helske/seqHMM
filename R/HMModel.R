#' Build a Hidden Markov Model
#' 
#' Function HHModel constructs an object of class \code{HMModel}.
#' 

#' @export
#' @useDynLib seqHMM
#' @param observations Observations. Can be either vector (single sequence ), data.frame 
#' (multiple sequences, each row corresponds to one sequence), or list of vectors/data.frames
#' where each component corresponds to one channel.
#' @param transitionMatrix Matrix of transition probabilities.
#' @param emissionMatrix Matrix of emission probabilities or in case multiple channels 
#' list of such matrices.
#' @param initialProbs Vector of initial state probabilities.
#' @param stateNames Optional labels for the observations.
#' @return Object of class \code{HMModel}
#' 
HMModel<-function(observations,transitionMatrix,emissionMatrix,initialProbs,
                  stateNames=NULL){
  
  # determine number of states
  numberOfStates<-length(initialProbs)
  if(is.null(stateNames))
    stateNames<-as.character(1:numberOfStates)
  channelNames<-NULL
  
  if(dim(transitionMatrix)[1]!=dim(transitionMatrix)[2] ||
       numberOfStates!=dim(transitionMatrix)[1])
    stop("Dimensions of transitionMatrix and length of initialProbs do not match.")
  
  if(!isTRUE(all.equal(rowSums(transitionMatrix),rep(1,dim(transitionMatrix)[1]))))
    stop("Transition probabilities in transitionMatrix do not sum to one.")
  if(!isTRUE(all.equal(sum(initialProbs),1)))
    stop("Initial state probabilites do not sum to one.")
  dimnames(transitionMatrix)<-list(from=stateNames,to=stateNames)
  
  # determine the type of HMM from observations:
  if(is.null(dim(observations)) && !is.list(observations)){
    numberOfSequences<-numberOfChannels<-as.integer(1)
    lengthOfSequences<-length(observations)
    symbolNames<-sort(unique(c(as.character(observations))))
    numberOfSymbols<-length(symbolNames)
    
    if(numberOfStates!=dim(emissionMatrix)[1])
      stop("Number of rows in emissionMatrix is not equal to the number of states.")
    if(numberOfSymbols!=dim(emissionMatrix)[2])
      stop("Number of columns in emissionMatrix is not equal to the number of symbols.")
    if(!isTRUE(all.equal(rep(1,numberOfStates),rowSums(emissionMatrix))))
      stop("Emission probabilities in emissionMatrix do not sum to one.")
    dimnames(emissionMatrix)<-list(stateNames=stateNames,symbolNames=symbolNames)
  } else {
    if(is.data.frame(observations)){
      numberOfSequences<-nrow(observations)
      lengthOfSequences<-ncol(observations)
      numberOfChannels<-as.integer(1)
      symbolNames<-sort(unique(as.character(unlist(observations))))
      numberOfSymbols<-length(symbolNames)
      
      
      if(numberOfStates!=dim(emissionMatrix)[1])
        stop("Number of rows in emissionMatrix is not equal to the number of states.")
      if(numberOfSymbols!=dim(emissionMatrix)[2])
        stop("Number of columns in emissionMatrix is not equal to the number of symbols.")
      if(!isTRUE(all.equal(rep(1,numberOfStates),rowSums(emissionMatrix))))
        stop("Emission probabilities in emissionMatrix do not sum to one.")
      dimnames(emissionMatrix)<-list(stateNames=stateNames,symbolNames=symbolNames)
      
      
    } else {
      if(is.list(observations)){
        numberOfChannels<-length(observations)
        if(is.data.frame(observations[[1]])){
          numberOfSequences<-nrow(observations[[1]])
          lengthOfSequences<-ncol(observations[[1]])
          symbolNames<-lapply(observations,function(x) sort(unique(as.character(unlist(x)))))
          numberOfSymbols<-sapply(symbolNames,length)
          
          
        } else {
          numberOfSequences<-as.integer(1)
          lengthOfSequences<-length(observations[[1]])
          symbolNames<-lapply(observations,function(x) sort(unique(c(as.character(x)))))
          numberOfSymbols<-sapply(symbolNames,length)
        }
        if(!isTRUE(all.equal(rep(numberOfStates,numberOfChannels),sapply(emissionMatrix,nrow))))
          stop("Number of rows in emissionMatrix is not equal to the number of states.")
        if(!isTRUE(all.equal(numberOfSymbols,sapply(emissionMatrix,ncol))))
          stop("Number of columns in emissionMatrix is not equal to the number of symbols.")
        if(!isTRUE(all.equal(c(sapply(emissionMatrix,rowSums)),rep(1,numberOfChannels*numberOfStates))))
          stop("Emission probabilities in emissionMatrix do not sum to one.")
        
        channelNames<-names(observations)  
        if(is.null(channelNames))
          channelNames<-as.character(1:numberOfChannels)
        for(i in 1:numberOfChannels)
          dimnames(emissionMatrix[[i]])<-list(stateNames=stateNames,symbolNames=symbolNames[[i]])
        names(emissionMatrix)<-channelNames
      } else stop("Cannot determine the type of HMM from provided observations.")
    }
  }
  
  
  
  model<-list(observations=observations,transitionMatrix=transitionMatrix,
              emissionMatrix=emissionMatrix,initialProbs=initialProbs,stateNames=stateNames,
              symbolNames=symbolNames,channelNames=channelNames,lengthOfSequences=lengthOfSequences,
              numberOfSymbols=numberOfSymbols,numberOfStates=numberOfStates,
              numberOfSequences=numberOfSequences,numberOfChannels=numberOfChannels)
  class(model)<-"HMModel"
  model
}