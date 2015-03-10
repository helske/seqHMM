#' Build a Mixture Hidden Markov Model
#' 
#' Function buildMixHMM constructs a mixture of hidden Markov models.
#' 
#' @import TraMineR
#' @importFrom Rcpp evalCpp
#' @export
#' @useDynLib seqHMM
#' @param observations TraMineR stslist (see \code{\link{seqdef}}) containing
#'   the sequences, or a list of such objects (one for each channel).
#' @param transitionMatrix List of matrices of transition 
#'   probabilities for each model.
#' @param emissionMatrix List which contains matrices of emission probabilities or
#'   a list of such objects (one for each channel) for each model. Note that the
#'   matrices must have dimensions m x s where m is the number of hidden states 
#'   and s is the number of unique symbols (observed states) in the data.
#' @param initialProbs List which contains vectors of initial state 
#'   probabilities for each model.
#' @param X Time-constant covariates for mixtures. Must be a matrix with 
#'   dimensions p x k where p is the number of sequences and k is the number of 
#'   covariates.
#' @param beta An k x l matrix of regression coefficients for time-constant 
#'   covariates for mixture probabilities, where l is the number of models and k
#'   is the number of covariates. A logit-link is used mixture probabilities.
#'   First column is set to zero.
#' @param stateNames List of optional labels for the hidden states
#' @return Object of class \code{mixtureHMModel}
#' @seealso \code{\link{fitMixHMM}} for fitting mixture Hidden Markov models.
#'   

buildMixHMM <- 
  function(observations,transitionMatrix,emissionMatrix,initialProbs, X, beta, stateNames=NULL){
    
    numberOfModels<-length(transitionMatrix)
    if(length(emissionMatrix)!=numberOfModels || length(initialProbs)!=numberOfModels)
      stop("Unequal lengths of transitionMatrix, emissionMatrix and initialProbs.")
    
    model <- vector("list", length = numberOfModels)
    
    # States
    numberOfStates <- unlist(lapply(transitionMatrix,nrow))
    
    if(any(numberOfStates!=unlist(lapply(transitionMatrix,dim))))
      stop("Transition matrices must be square matrices with same dimensions.")
    
    if(is.null(stateNames)){
      stateNames <- vector("list", numberOfModels)
      for(m in 1:numberOfModels){
        stateNames[[m]] <- as.character(1:numberOfStates[m])
      }
    }
    
    if(!all(1==unlist(sapply(transitionMatrix,rowSums))))
      stop("Transition probabilities in transitionMatrix do not sum to one.")
    
    if(!all(1==unlist(sapply(initialProbs,sum))))
      stop("Initial state probabilities do not sum to one.")

    for(i in 1:numberOfModels){

      dimnames(transitionMatrix[[i]]) <- list(from=stateNames[[i]],to=stateNames[[i]])
      # Single channel but emissionMatrix is list of lists  
      if(is.list(emissionMatrix[[i]]) && length(emissionMatrix[[i]])==1)   
        emissionMatrix[[i]] <- emissionMatrix[[i]][[1]]
    }
    
    
    
    
    # Single channel but observations is a list
    if(is.list(observations) && length(observations)==1)
      observations <- observations[[1]]
    
    numberOfChannels <- ifelse(is.list(emissionMatrix[[1]]),length(emissionMatrix[[1]]),1)
    
    if(any(sapply(emissionMatrix,length)!=numberOfChannels))
      stop("Number of channels defined by emission matrices differ from each other.")
    
    if(numberOfChannels>1){
      if(length(observations)!=numberOfChannels){
        stop("Number of channels defined by emissionMatrix differs from one defined by observations.")
      }
      
      
      numberOfSequences<-nrow(observations[[1]])
      lengthOfSequences<-ncol(observations[[1]])
      

      symbolNames<-lapply(observations,alphabet)
      numberOfSymbols<-sapply(symbolNames,length)
      for(i in 1:numberOfModels){
        if(any(lapply(emissionMatrix[[i]],nrow)!=numberOfStates[i]))
          stop("Number of rows in emissionMatrix is not equal to the number of states.")
        
        if(any(numberOfSymbols!=sapply(emissionMatrix[[i]],ncol)))
          stop("Number of columns in emissionMatrix is not equal to the number of symbols.")
        if(!isTRUE(all.equal(c(sapply(emissionMatrix[[i]],rowSums)),
                             rep(1,numberOfChannels*numberOfStates[i]),check.attributes=FALSE)))
          stop("Emission probabilities in emissionMatrix do not sum to one.")

        channelNames<-names(observations)  
        if(is.null(channelNames))
          channelNames<-as.character(1:numberOfChannels)
        for(j in 1:numberOfChannels)
          dimnames(emissionMatrix[[i]][[j]])<-list(stateNames=stateNames[[i]],symbolNames=symbolNames[[j]])
        names(emissionMatrix[[i]])<-channelNames
      }
    } else {
      numberOfChannels <- 1
      channelNames<-NULL
      numberOfSequences<-nrow(observations)
      lengthOfSequences<-ncol(observations)
      symbolNames<-alphabet(observations)
      numberOfSymbols<-length(symbolNames)
      
      for(i in 1:numberOfModels){
        if(numberOfStates[i]!=dim(emissionMatrix[[i]])[1])
          stop("Number of rows in emissionMatrix is not equal to the number of states.")
        if(numberOfSymbols!=dim(emissionMatrix[[i]])[2])
          stop("Number of columns in emissionMatrix is not equal to the number of symbols.")
        if(!isTRUE(all.equal(rep(1,numberOfStates[i]),rowSums(emissionMatrix[[i]]),check.attributes=FALSE)))
          stop("Emission probabilities in emissionMatrix do not sum to one.")
        dimnames(emissionMatrix[[i]])<-list(stateNames=stateNames[[i]],symbolNames=symbolNames)
      }
      
    }
    
    if(!missing(X)){
      if(nrow(X)!=numberOfSequences)
        stop("Wrong dimensions of X.")
      numberOfCovariates<-ncol(X)
      if(missing(beta)){
        beta<-matrix(0,numberOfCovariates,numberOfModels)
      } else {
        if(ncol(beta)!=numberOfModels | nrow(beta)!=numberOfCovariates)
          stop("Wrong dimensions of beta.")
        beta[,1]<-0
      }       
    } else { #Just intercept
      numberOfCovariates <-1
      X <- matrix(1,nrow=numberOfSequences)
      beta <- matrix(0,1,numberOfModels)        
    }
    
    model<-list(observations=observations, transitionMatrix=transitionMatrix,
                emissionMatrix=emissionMatrix, initialProbs=initialProbs,
                beta=beta, X=X,stateNames=stateNames, symbolNames=symbolNames,
                channelNames=channelNames, lengthOfSequences=lengthOfSequences,
                numberOfSequences=numberOfSequences, numberOfModels=numberOfModels,
                numberOfSymbols=numberOfSymbols, numberOfStates=numberOfStates,
                numberOfChannels=numberOfChannels,
                numberOfCovariates=numberOfCovariates)
    class(model)<-"mixHMModel"
    model
  }