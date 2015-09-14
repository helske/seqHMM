#' Build a Hidden Markov Model
#' 
#' Function build_hmm constructs an object of class \code{hmm}.
#' 
#' @importFrom Rcpp evalCpp
#' @importFrom TraMineR alphabet
#' @export
#' @useDynLib seqHMM
#' @param observations TraMineR stslist (see \code{\link[TraMineR]{seqdef}}) containing 
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
#' # Single-channel data
#'
#' data(mvad)
#' 
#' mvad.alphabet <- c("employment", "FE", "HE", "joblessness", "school", 
#'                    "training")
#' mvad.labels <- c("employment", "further education", "higher education", 
#'                  "joblessness", "school", "training")
#' mvad.scodes <- c("EM", "FE", "HE", "JL", "SC", "TR")
#' mvad.seq <- seqdef(mvad, 17:86, alphabet = mvad.alphabet, states = mvad.scodes, 
#'                    labels = mvad.labels, xtstep = 6)
#' 
#' # Starting values for the emission matrix
#' B <- matrix(NA, nrow = 4, ncol = 6)
#' B[1,] <- seqstatf(mvad.seq[, 1:12])[, 2] + 0.1
#' B[2,] <- seqstatf(mvad.seq[, 13:24])[, 2] + 0.1
#' B[3,] <- seqstatf(mvad.seq[, 25:48])[, 2] + 0.1
#' B[4,] <- seqstatf(mvad.seq[, 49:70])[, 2] + 0.1
#' B <- B / rowSums(B)
#' 
#' # Starting values for the transition matrix
#' 
#' A <-  matrix(c(0.80, 0.10, 0.05, 0.05,
#'                0.05, 0.80, 0.10, 0.05,
#'                0.05, 0.05, 0.80, 0.10,
#'                0.05, 0.05, 0.10, 0.80), nrow=4, ncol=4, byrow=TRUE)
#' 
#' # Starting values for initial state probabilities
#' initial_probs <- c(0.3, 0.3, 0.2, 0.2)
#' 
#' # Building a hidden Markov model with starting values
#' bHMM <- build_hmm(
#'   observations = mvad.seq, transition_matrix = A, 
#'   emission_matrix = B, initial_probs = initial_probs
#' )
#' 
#' #########################################
#' 
#' 
#' # Multichannel data
#' 
#' data(biofam3c)
#' 
#' # Building sequence objects
#' child.seq <- seqdef(biofam3c$children)
#' marr.seq <- seqdef(biofam3c$married)
#' left.seq <- seqdef(biofam3c$left)
#' 
#' # Starting values for emission matrices
#' B_child <- matrix(NA, nrow = 3, ncol = 2)
#' B_child[1,] <- seqstatf(child.seq[, 1:5])[, 2] + 0.1
#' B_child[2,] <- seqstatf(child.seq[, 6:10])[, 2] + 0.1
#' B_child[3,] <- seqstatf(child.seq[, 11:15])[, 2] + 0.1
#' B_child <- B_child / rowSums(B_child)
#' 
#' B_marr <- matrix(NA, nrow = 3, ncol = 3)
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
  n_states<-nrow(transition_matrix)
  if(is.null(state_names))
    state_names<-as.character(1:n_states)
  
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
    n_channels <- length(emission_matrix)
    
    n_sequences<-nrow(observations[[1]])
    length_of_sequences<-ncol(observations[[1]])
    
    symbol_names<-lapply(observations,alphabet)
    n_symbols<-sapply(symbol_names,length)
    
    if(any(sapply(emission_matrix,nrow)!=n_states))
      stop("Number of rows in emission_matrix is not equal to the number of states.")
    if(any(n_symbols!=sapply(emission_matrix,ncol)))
      stop("Number of columns in emission_matrix is not equal to the number of symbols.")
    if(!isTRUE(all.equal(c(sapply(emission_matrix,rowSums)),rep(1,n_channels*n_states),check.attributes=FALSE)))
      stop("Emission probabilities in emission_matrix do not sum to one.")
    
    if(is.null(channel_names)){
      channel_names<-as.character(1:n_channels)
    }else if(length(channel_names)!=n_channels){
      warning("The length of argument channel_names does not match the number of channels. Names were not used.")
      channel_names<-as.character(1:n_channels)
    }
    for(i in 1:n_channels)
      dimnames(emission_matrix[[i]])<-list(state_names=state_names,symbol_names=symbol_names[[i]])
    names(emission_matrix)<-channel_names
  } else {
    n_channels <- 1
    channel_names<-NULL
    n_sequences<-nrow(observations)
    length_of_sequences<-ncol(observations)
    symbol_names<-alphabet(observations)
    n_symbols<-length(symbol_names)
    
    if(n_states!=dim(emission_matrix)[1])
      stop("Number of rows in emission_matrix is not equal to the number of states.")
    if(n_symbols!=dim(emission_matrix)[2])
      stop("Number of columns in emission_matrix is not equal to the number of symbols.")
    if(!isTRUE(all.equal(rep(1,n_states),rowSums(emission_matrix),check.attributes=FALSE)))
      stop("Emission probabilities in emission_matrix do not sum to one.")
    dimnames(emission_matrix)<-list(state_names=state_names,symbol_names=symbol_names)
    
  }  
  
  names(initial_probs) <- state_names
  
  model<-list(observations=observations,transition_matrix=transition_matrix,
              emission_matrix=emission_matrix,initial_probs=initial_probs,
              state_names=state_names,
              symbol_names=symbol_names,channel_names=channel_names,
              length_of_sequences=length_of_sequences,
              n_sequences=n_sequences,
              n_symbols=n_symbols,n_states=n_states,
              n_channels=n_channels)
  class(model)<-"hmm"
  model
}