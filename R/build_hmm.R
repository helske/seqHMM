#' Build a Hidden Markov Model
#'
#' Function \code{build_hmm} constructs a hidden Markov model object of class \code{hmm}.
#'
#' The returned model contains some attributes such as \code{nobs} and \code{df},
#' which define the number of observations in the  model and the number of estimable
#' model parameters, used in computing BIC.
#' When computing \code{nobs} for a multichannel model with \eqn{C} channels, 
#' each observed value in a single channel amounts to \eqn{1/C} observation, 
#' i.e. a fully observed time point for a single sequence amounts to one observation. 
#' For the degrees of freedom \code{df}, zero probabilities of the initial model are 
#' defined as structural zeroes.
#' @export
#' @param observations An \code{stslist} object (see \code{\link[TraMineR]{seqdef}}) containing
#' the sequences, or a list of such objects (one for each channel).
#' @param transition_probs A matrix of transition probabilities.
#' @param emission_probs A matrix of emission probabilities or a list of such
#' objects (one for each channel). Emission probabilities should follow the
#' ordering of the alphabet of observations (\code{alphabet(observations)}, returned as \code{symbol_names}).
#' @param initial_probs A vector of initial state probabilities.
#' @param state_names A list of optional labels for the hidden states. If \code{NULL},
#' the state names are taken from the row names of the transition matrix. If this is
#' also \code{NULL}, numbered states are used.
#' @param channel_names A vector of optional names for the channels.
#' @return Object of class \code{hmm} with the following elements:
#' \describe{
#'    \item{\code{observations}}{State sequence object or a list of such objects containing the data.}
#'    \item{\code{transition_probs}}{A matrix of transition probabilities.}
#'    \item{\code{emission_probs}}{A matrix or a list of matrices of emission probabilities.}
#'    \item{\code{initial_probs}}{A vector of initial probabilities.}
#'    \item{\code{state_names}}{Names for hidden states.}
#'    \item{\code{symbol_names}}{Names for observed states.}
#'    \item{\code{channel_names}}{Names for channels of sequence data.}
#'    \item{\code{length_of_sequences}}{(Maximum) length of sequences.}
#'    \item{\code{n_sequences}}{Number of sequences.}
#'    \item{\code{n_symbols}}{Number of observed states (in each channel).}
#'    \item{\code{n_states}}{Number of hidden states.}
#'    \item{\code{n_channels}}{Number of channels.}
#'}
#'
#' @seealso \code{\link{fit_model}} for estimating model parameters; and 
#'   \code{\link{plot.hmm}} for plotting \code{hmm} objects.
#' @examples
#'
#' # Single-channel data
#'
#' data("mvad", package = "TraMineR")
#'
#' mvad_alphabet <- c("employment", "FE", "HE", "joblessness", "school",
#'                    "training")
#' mvad_labels <- c("employment", "further education", "higher education",
#'                  "joblessness", "school", "training")
#' mvad_scodes <- c("EM", "FE", "HE", "JL", "SC", "TR")
#' mvad_seq <- seqdef(mvad, 17:86, alphabet = mvad_alphabet, states = mvad_scodes,
#'                    labels = mvad_labels, xtstep = 6)
#'
#' # Starting values for the emission matrix
#' emiss <- matrix(NA, nrow = 4, ncol = 6)
#' emiss[1,] <- seqstatf(mvad_seq[, 1:12])[, 2] + 1
#' emiss[2,] <- seqstatf(mvad_seq[, 13:24])[, 2] + 1
#' emiss[3,] <- seqstatf(mvad_seq[, 25:48])[, 2] + 1
#' emiss[4,] <- seqstatf(mvad_seq[, 49:70])[, 2] + 1
#' emiss <- emiss / rowSums(emiss)
#'
#' # Starting values for the transition matrix
#'
#' tr <- matrix(
#'   c(0.80, 0.10, 0.05, 0.05,
#'     0.05, 0.80, 0.10, 0.05,
#'     0.05, 0.05, 0.80, 0.10,
#'     0.05, 0.05, 0.10, 0.80),
#'   nrow=4, ncol=4, byrow=TRUE)
#'
#' # Starting values for initial state probabilities
#' init <- c(0.3, 0.3, 0.2, 0.2)
#'
#' # Building a hidden Markov model with starting values
#' init_hmm_mvad <- build_hmm(
#'   observations = mvad_seq, transition_probs = tr,
#'   emission_probs = emiss, initial_probs = init
#' )
#'
#' #########################################
#'
#'
#' # Multichannel data
#'
#' # Three-state three-channel hidden Markov model
#' # See ?hmm_biofam for a five-state version
#'
#' data("biofam3c")
#'
#' # Building sequence objects
#' marr_seq <- seqdef(biofam3c$married, start = 15,
#'   alphabet = c("single", "married", "divorced"))
#' child_seq <- seqdef(biofam3c$children, start = 15,
#'   alphabet = c("childless", "children"))
#' left_seq <- seqdef(biofam3c$left, start = 15,
#'   alphabet = c("with parents", "left home"))
#'
#' # Define colors
#' attr(marr_seq, "cpal") <- c("violetred2", "darkgoldenrod2", "darkmagenta")
#' attr(child_seq, "cpal") <- c("darkseagreen1", "coral3")
#' attr(left_seq, "cpal") <- c("lightblue", "red3")
#'
#' # Starting values for emission matrices
#'
#' emiss_marr <- matrix(NA, nrow = 3, ncol = 3)
#' emiss_marr[1,] <- seqstatf(marr_seq[, 1:5])[, 2] + 1
#' emiss_marr[2,] <- seqstatf(marr_seq[, 6:10])[, 2] + 1
#' emiss_marr[3,] <- seqstatf(marr_seq[, 11:16])[, 2] + 1
#' emiss_marr <- emiss_marr / rowSums(emiss_marr)
#'
#' emiss_child <- matrix(NA, nrow = 3, ncol = 2)
#' emiss_child[1,] <- seqstatf(child_seq[, 1:5])[, 2] + 1
#' emiss_child[2,] <- seqstatf(child_seq[, 6:10])[, 2] + 1
#' emiss_child[3,] <- seqstatf(child_seq[, 11:16])[, 2] + 1
#' emiss_child <- emiss_child / rowSums(emiss_child)
#'
#' emiss_left <- matrix(NA, nrow = 3, ncol = 2)
#' emiss_left[1,] <- seqstatf(left_seq[, 1:5])[, 2] + 1
#' emiss_left[2,] <- seqstatf(left_seq[, 6:10])[, 2] + 1
#' emiss_left[3,] <- seqstatf(left_seq[, 11:16])[, 2] + 1
#' emiss_left <- emiss_left / rowSums(emiss_left)
#'
#' # Starting values for transition matrix
#' trans <- matrix(c(0.9, 0.07, 0.03,
#'                   0,    0.9,  0.1,
#'                   0,      0,    1), 
#'   nrow = 3, ncol = 3, byrow = TRUE)
#'
#' # Starting values for initial state probabilities
#' inits <- c(0.9, 0.09, 0.01)
#'
#' # Building hidden Markov model with initial parameter values
#' init_hmm_bf <- build_hmm(
#'   observations = list(marr_seq, child_seq, left_seq),
#'   transition_probs = trans,
#'   emission_probs = list(emiss_marr, emiss_child, emiss_left),
#'   initial_probs = inits)
#'

build_hmm <- function(observations, transition_probs, emission_probs, initial_probs,
  state_names = NULL, channel_names = NULL){

  if (!is.matrix(transition_probs)) {
    stop(paste("Object provided for transition_probs is not a matrix."))
  }
  if (!is.vector(initial_probs)) {
    stop(paste("Object provided for initial_probs is not a vector."))
  }

  if(dim(transition_probs)[1]!=dim(transition_probs)[2])
    stop("Transition_probs must be a square matrix.")
  n_states <- nrow(transition_probs)

  if (length(initial_probs) != n_states){
    stop(paste("Length of initial_probs is not equal to the number of states."))
  }

  if (is.null(state_names)) {
    if (is.null(state_names <- rownames(transition_probs))) {
      state_names <- paste("State", 1:n_states)
    }
  } else {
    if (length(state_names) != n_states) {
      stop("Length of state_names is not equal to the number of hidden states.")
    }
  }

  if(!isTRUE(all.equal(rowSums(transition_probs),rep(1,dim(transition_probs)[1]),check.attributes=FALSE)))
    stop("Transition probabilities in transition_probs do not sum to one.")

  dimnames(transition_probs)<-list(from=state_names,to=state_names)

  if(is.list(emission_probs) && length(emission_probs)==1){
    emission_probs <- emission_probs[[1]]
  }
  if(is.list(observations) && !inherits(observations, "stslist") && length(observations)==1){
    observations <- observations[[1]]
  }



  if(is.list(emission_probs)){
    if(length(observations)!=length(emission_probs)){
      stop("Number of channels defined by emission_probs differs from one defined by observations.")
    }
    n_channels <- length(emission_probs)
    for (j in 1:n_channels){
      if (!is.matrix(emission_probs[[j]])) {
        stop(paste("Object provided in emission_probs for channel", j, "is not a matrix."))
      }
    }

    if (length(unique(sapply(observations, nrow))) > 1) {
      stop("The number of subjects (rows) is not the same in all channels.")
    }
    if (length(unique(sapply(observations, ncol))) > 1) {
      stop("The length of the sequences (number of columns) is not the same in all channels.")
    }

    n_sequences <- nrow(observations[[1]])
    length_of_sequences <- ncol(observations[[1]])

    symbol_names <- lapply(observations,alphabet)
    n_symbols <- lengths(symbol_names)

    if (any(sapply(emission_probs,nrow) != n_states))
      stop("Number of rows in emission_probs is not equal to the number of states.")
    if (any(n_symbols != sapply(emission_probs,ncol)))
      stop("Number of columns in emission_probs is not equal to the number of symbols.")
    if (!isTRUE(all.equal(c(sapply(emission_probs,rowSums)),rep(1, n_channels * n_states), check.attributes = FALSE)))
      stop("Emission probabilities in emission_probs do not sum to one.")

    if(is.null(channel_names)) {
      if(is.null(channel_names <- names(observations))){
        channel_names <- paste("Channel", 1:n_channels)
      }
    }else if(length(channel_names)!=n_channels){
      warning("The length of argument channel_names does not match the number of channels. Names were not used.")
      channel_names<-paste("Channel", 1:n_channels)
    }
    for(i in 1:n_channels)
      dimnames(emission_probs[[i]])<-list(state_names=state_names,symbol_names=symbol_names[[i]])
    names(emission_probs)<-channel_names
  } else {
    n_channels <- 1
    if (!is.matrix(emission_probs)) {
      stop(paste("Object provided for emission_probs is not a matrix."))
    }
    if (is.null(channel_names)) {
      channel_names <- "Observations"
    }
    n_sequences<-nrow(observations)
    length_of_sequences<-ncol(observations)
    symbol_names<-alphabet(observations)
    n_symbols<-length(symbol_names)

    if(n_states!=dim(emission_probs)[1])
      stop("Number of rows in emission_probs is not equal to the number of states.")
    if(n_symbols!=dim(emission_probs)[2])
      stop("Number of columns in emission_probs is not equal to the number of symbols.")
    if(!isTRUE(all.equal(rep(1,n_states),rowSums(emission_probs),check.attributes=FALSE)))
      stop("Emission probabilities in emission_probs do not sum to one.")
    dimnames(emission_probs)<-list(state_names=state_names,symbol_names=symbol_names)

  }

  names(initial_probs) <- state_names

  if(n_channels > 1){
    nobs <- sum(sapply(observations, function(x) sum(!(x == attr(observations[[1]], "nr") |
        x == attr(observations[[1]], "void") |
        is.na(x)))))/n_channels
  } else {
    nobs <- sum(!(observations == attr(observations, "nr") |
        observations == attr(observations, "void") |
        is.na(observations)))
  }

  model <- structure(list(observations=observations,transition_probs=transition_probs,
    emission_probs=emission_probs,initial_probs=initial_probs,
    state_names=state_names,
    symbol_names=symbol_names,channel_names=channel_names,
    length_of_sequences=length_of_sequences,
    n_sequences=n_sequences,
    n_symbols=n_symbols,n_states=n_states,
    n_channels=n_channels), class = "hmm",
    nobs = nobs,
    df = sum(initial_probs > 0) - 1 + sum(transition_probs > 0) - n_states +
      sum(unlist(emission_probs) > 0) - n_states * n_channels,
    type = "hmm")

  model
}
