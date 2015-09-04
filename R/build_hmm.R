#' Build a Hidden Markov Model
#' 
#' Function build_hmm constructs an object of class \code{hmm}.
#' 
#' @import TraMineR
#' @importFrom Rcpp evalCpp
#' @export
#' @useDynLib seqHMM
#' @param observations TraMineR stslist (see \code{\link{seqdef}}) containing 
#' the sequences, or a list of such objects (one for each channel).
#' @param transition_matrix A matrix of transition probabilities.
#' @param emission_matrix A matrix of emission probabilities or a list of such 
#' objects (one for each channel).
#' @param initial_probs A vector of initial state probabilities.
#' @param state_names A list of optional labels for the hidden states.
#' @param channel_names A vector of optional names for the channels.
#' @return Object of class \code{hmm}
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
#' initial_probs <- c(0.9, 0.07, 0.02, 0.01)
#' 
#' # Building a hidden Markov model with starting values
#' bHMM <- build_hmm(
#'   observations = biofam.seq, transition_matrix = A, 
#'   emission_matrix = B, initial_probs = initial_probs
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
#' # Starting values for emission matrices
#' B_child <- matrix(NA, nrow = 3, ncol = 2)
#' B_child[1,] <- seqstatf(child.seq[, 1:5])[, 2] + 0.1
#' B_child[2,] <- seqstatf(child.seq[, 6:10])[, 2] + 0.1
#' B_child[3,] <- seqstatf(child.seq[, 11:15])[, 2] + 0.1
#' B_child <- B_child / rowSums(B_child)
#' 
#' B_marr <- matrix(NA, nrow = 3, ncol = 2)
#' B_marr[1,] <- seqstatf(marr.seq[, 1:5])[, 2] + 0.1
#' B_marr[2,] <- seqstatf(marr.seq[, 6:10])[, 2] + 0.1
#' B_marr[3,] <- seqstatf(marr.seq[, 11:15])[, 2] + 0.1
#' B_marr <- B_marr / rowSums(B_marr)
#' 
#' B_left <- matrix(NA, nrow = 3, ncol = 2)
#' B_left[1,] <- seqstatf(left.seq[, 1:5])[, 2] + 0.1
#' B_left[2,] <- seqstatf(left.seq[, 6:10])[, 2] + 0.1
#' B_left[3,] <- seqstatf(left.seq[, 11:15])[, 2] + 0.1
#' B_left <- B_left / rowSums(B_left)
#' 
#' # Starting values for transition matrix
#' A <- matrix(c(0.9, 0.07, 0.03,
#'                 0,  0.9,  0.1,
#'                 0,    0,    1), nrow = 3, ncol = 3, byrow = TRUE)
#' 
#' # Starting values for initial state probabilities
#' init <- c(0.9, 0.09, 0.01)
#' 
#' # Building hidden Markov model with initial parameter values
#' bHMM <- build_hmm(
#'   observations = list(child.seq, marr.seq, left.seq), 
#'   transition_matrix = A,
#'   emission_matrix = list(B_child, B_marr, B_left), 
#'   initial_probs = init
#'   )
#' 
#' @seealso \code{\link{fit_hmm}} for fitting Hidden Markov models.

build_hmm<-function(observations,transition_matrix,emission_matrix,initial_probs,
                   state_names=NULL, channel_names=NULL){
  
  
  if(dim(transition_matrix)[1]!=dim(transition_matrix)[2])
    stop("transition_matrix must be a square matrix.")
  number_of_states<-nrow(transition_matrix)
  if(is.null(state_names))
    state_names<-as.character(1:number_of_states)
  
  if(!isTRUE(all.equal(rowSums(transition_matrix),rep(1,dim(transition_matrix)[1]),check.attributes=FALSE)))
    stop("Transition probabilities in transition_matrix do not sum to one.")

  dimnames(transition_matrix)<-list(from=state_names,to=state_names)
  
  if(is.list(emission_matrix) && length(emission_matrix)==1){
    emission_matrix <- emission_matrix[[1]]   
  }
  if(is.list(observations) && !inherits(observations, "stslist") && length(observations)==1){
    observations <- observations[[1]]
  }
  
  
  
  if(is.list(emission_matrix)){
    if(length(observations)!=length(emission_matrix)){
      stop("Number of channels defined by emission_matrix differs from one defined by observations.")
    }
    number_of_channels <- length(emission_matrix)
    
    number_of_sequences<-nrow(observations[[1]])
    length_of_sequences<-ncol(observations[[1]])
    
    symbol_names<-lapply(observations,alphabet)
    number_of_symbols<-sapply(symbol_names,length)
    
    if(any(sapply(emission_matrix,nrow)!=number_of_states))
      stop("Number of rows in emission_matrix is not equal to the number of states.")
    if(any(number_of_symbols!=sapply(emission_matrix,ncol)))
      stop("Number of columns in emission_matrix is not equal to the number of symbols.")
    if(!isTRUE(all.equal(c(sapply(emission_matrix,rowSums)),rep(1,number_of_channels*number_of_states),check.attributes=FALSE)))
      stop("Emission probabilities in emission_matrix do not sum to one.")
    
    if(is.null(channel_names)){
      channel_names<-as.character(1:number_of_channels)
    }else if(length(channel_names)!=number_of_channels){
      warning("The length of argument channel_names does not match the number of channels. Names were not used.")
      channel_names<-as.character(1:number_of_channels)
    }
    for(i in 1:number_of_channels)
      dimnames(emission_matrix[[i]])<-list(state_names=state_names,symbol_names=symbol_names[[i]])
    names(emission_matrix)<-channel_names
  } else {
    number_of_channels <- 1
    channel_names<-NULL
    number_of_sequences<-nrow(observations)
    length_of_sequences<-ncol(observations)
    symbol_names<-alphabet(observations)
    number_of_symbols<-length(symbol_names)
    
    if(number_of_states!=dim(emission_matrix)[1])
      stop("Number of rows in emission_matrix is not equal to the number of states.")
    if(number_of_symbols!=dim(emission_matrix)[2])
      stop("Number of columns in emission_matrix is not equal to the number of symbols.")
    if(!isTRUE(all.equal(rep(1,number_of_states),rowSums(emission_matrix),check.attributes=FALSE)))
      stop("Emission probabilities in emission_matrix do not sum to one.")
    dimnames(emission_matrix)<-list(state_names=state_names,symbol_names=symbol_names)
    
  }  
  
  model<-list(observations=observations,transition_matrix=transition_matrix,
              emission_matrix=emission_matrix,initial_probs=initial_probs,
              state_names=state_names,
              symbol_names=symbol_names,channel_names=channel_names,
              length_of_sequences=length_of_sequences,
              number_of_sequences=number_of_sequences,
              number_of_symbols=number_of_symbols,number_of_states=number_of_states,
              number_of_channels=number_of_channels)
  class(model)<-"hmm"
  model
}