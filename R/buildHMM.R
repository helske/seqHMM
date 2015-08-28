#' Build a Hidden Markov Model
#' 
#' Function buildHMM constructs an object of class \code{HMModel}.
#' 
#' @import TraMineR
#' @importFrom Rcpp evalCpp
#' @export
#' @useDynLib seqHMM
#' @param observations TraMineR stslist (see \code{\link{seqdef}}) containing 
#' the sequences, or a list of such objects (one for each channel).
#' @param transitionMatrix A matrix of transition probabilities.
#' @param emissionMatrix A matrix of emission probabilities or a list of such 
#' objects (one for each channel).
#' @param initialProbs A vector of initial state probabilities.
#' @param stateNames A list of optional labels for the hidden states.
#' @param channelNames A vector of optional names for the channels.
#' @return Object of class \code{HMModel}
#' 
#' @examples 
#' require(TraMineR)
#' 
#' data(biofam)
#' biofam <- biofam[1:500,]
#' 
#' # Single-channel data
#' 
#' biofam.seq <- seqdef(
#'   biofam[, 10:25], 
#'   states = c("Parent", "Left", "Married", "Left+Marr",
#'              "Left+Child", "Left+Marr+Child", "Divorced"),
#'   start = 15
#'   )
#' 
#' # Starting values for the emission matrix
#' B <- matrix(NA, nrow = 4, ncol = 7)
#' B[1,] <- seqstatf(biofam.seq[, 1:4])[, 2] + 0.1
#' B[2,] <- seqstatf(biofam.seq[, 5:8])[, 2] + 0.1
#' B[3,] <- seqstatf(biofam.seq[, 9:12])[, 2] + 0.1
#' B[4,] <- seqstatf(biofam.seq[, 13:15])[, 2] + 0.1
#' B <- B / rowSums(B)
#' 
#' # Starting values for the transition matrix
#' A <- matrix(c(0.80, 0.10, 0.05, 0.05,
#'               0.05, 0.80, 0.10, 0.05,
#'               0.05, 0.05, 0.80, 0.10,
#'               0.05, 0.05, 0.10, 0.80), nrow=4, ncol=4, byrow=TRUE)
#' 
#' # Starting values for initial state probabilities
#' initialProbs <- c(0.9, 0.07, 0.02, 0.01)
#' 
#' # Building a hidden Markov model with starting values
#' bHMM <- buildHMM(
#'   observations = biofam.seq, transitionMatrix = A, 
#'   emissionMatrix = B, initialProbs = initialProbs
#' )
#' 
#' #########################################
#' 
#' # Multichannel data
#' 
#' # Building one channel per type of event left, children or married
#' bf <- as.matrix(biofam[, 10:25])
#' children <-  bf==4 | bf==5 | bf==6
#' married <- bf == 2 | bf== 3 | bf==6
#' left <- bf==1 | bf==3 | bf==5 | bf==6
#' 
#' children[children==TRUE] <- "Children"
#' children[children==FALSE] <- "Childless"
#' 
#' married[married==TRUE] <- "Married"
#' married[married==FALSE] <- "Single"
#' 
#' left[left==TRUE] <- "Left home"
#' left[left==FALSE] <- "With parents"
#' 
#' # Building sequence objects
#' child.seq <- seqdef(children)
#' marr.seq <- seqdef(married)
#' left.seq <- seqdef(left)
#' 
#' # Initial values for emission matrices
#' B_child <- matrix(NA, nrow=3, ncol=2)
#' B_child[1,] <- seqstatf(child.seq[,1:5])[,2]+0.1
#' B_child[2,] <- seqstatf(child.seq[,6:10])[,2]+0.1
#' B_child[3,] <- seqstatf(child.seq[,11:15])[,2]+0.1
#' B_child <- B_child/rowSums(B_child)
#' 
#' B_marr <- matrix(NA, nrow=3, ncol=2)
#' B_marr[1,] <- seqstatf(marr.seq[,1:5])[,2]+0.1
#' B_marr[2,] <- seqstatf(marr.seq[,6:10])[,2]+0.1
#' B_marr[3,] <- seqstatf(marr.seq[,11:15])[,2]+0.1
#' B_marr <- B_marr/rowSums(B_marr)
#' 
#' B_left <- matrix(NA, nrow=3, ncol=2)
#' B_left[1,] <- seqstatf(left.seq[,1:5])[,2]+0.1
#' B_left[2,] <- seqstatf(left.seq[,6:10])[,2]+0.1
#' B_left[3,] <- seqstatf(left.seq[,11:15])[,2]+0.1
#' B_left <- B_left/rowSums(B_left)
#' 
#' # Initial values for transition matrix
#' A <- matrix(c(0.9, 0.07, 0.03,
#' 0,    0.9,  0.1,
#' 0,      0,    1), 
#' nrow=3, ncol=3, byrow=TRUE)
#' 
#' # Initial values for initial state probabilities
#' initialProbs <- c(0.9,0.09,0.01)
#' 
#' # Building hidden Markov model with initial parameter values
#' bHMM <- buildHMM(observations=list(child.seq, marr.seq, left.seq), 
#' transitionMatrix=A,
#' emissionMatrix=list(B_child, B_marr, B_left), 
#' initialProbs=initialProbs)
#' 
#' @seealso \code{\link{fitHMM}} for fitting Hidden Markov models.

buildHMM<-function(observations,transitionMatrix,emissionMatrix,initialProbs,
                   stateNames=NULL, channelNames=NULL){
  
  
  if(dim(transitionMatrix)[1]!=dim(transitionMatrix)[2])
    stop("transitionMatrix must be a square matrix.")
  numberOfStates<-nrow(transitionMatrix)
  if(is.null(stateNames))
    stateNames<-as.character(1:numberOfStates)
  
  if(!isTRUE(all.equal(rowSums(transitionMatrix),rep(1,dim(transitionMatrix)[1]),check.attributes=FALSE)))
    stop("Transition probabilities in transitionMatrix do not sum to one.")

  dimnames(transitionMatrix)<-list(from=stateNames,to=stateNames)
  
  if(is.list(emissionMatrix) && length(emissionMatrix)==1){
    emissionMatrix <- emissionMatrix[[1]]   
  }
  if(is.list(observations) && length(observations)==1){
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
    
    if(any(sapply(emissionMatrix,nrow)!=numberOfStates))
      stop("Number of rows in emissionMatrix is not equal to the number of states.")
    if(any(numberOfSymbols!=sapply(emissionMatrix,ncol)))
      stop("Number of columns in emissionMatrix is not equal to the number of symbols.")
    if(!isTRUE(all.equal(c(sapply(emissionMatrix,rowSums)),rep(1,numberOfChannels*numberOfStates),check.attributes=FALSE)))
      stop("Emission probabilities in emissionMatrix do not sum to one.")
    
    if(is.null(channelNames)){
      channelNames<-as.character(1:numberOfChannels)
    }else if(length(channelNames)!=numberOfChannels){
      warning("The length of argument channelNames does not match the number of channels. Names were not used.")
      channelNames<-as.character(1:numberOfChannels)
    }
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
  
  model<-list(observations=observations,transitionMatrix=transitionMatrix,
              emissionMatrix=emissionMatrix,initialProbs=initialProbs,stateNames=stateNames,
              symbolNames=symbolNames,channelNames=channelNames,lengthOfSequences=lengthOfSequences,
              numberOfSequences=numberOfSequences,
              numberOfSymbols=numberOfSymbols,numberOfStates=numberOfStates,
              numberOfChannels=numberOfChannels)
  class(model)<-"HMModel"
  model
}