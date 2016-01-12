#' Build a Markov Model
#'
#' Function \code{build_mm} is a shortcut for constructing a Markov model as a restricted 
#'   case of an \code{hmm} object.
#' @export
#' @param observations An \code{stslist} object (see \code{\link[TraMineR]{seqdef}}) containing
#' the sequences.
#' @param transition_probs A matrix of transition probabilities. Transition probabilities should follow the
#' ordering of the alphabet of observations (\code{alphabet(observations)}, returned as \code{state_names}).
#' @param initial_probs A vector of initial state probabilities.
#' @return Object of class \code{hmm} with following elements:
#' \describe{
#'    \item{\code{observations}}{State sequence object or a list of such containing the data.}
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
#' @examples
#' # Construct sequence data
#' data("mvad", package = "TraMineR")
#'
#' mvad_alphabet <-
#'   c("employment", "FE", "HE", "joblessness", "school", "training")
#' mvad_labels <- c("employment", "further education", "higher education",
#'   "joblessness", "school", "training")
#' mvad_scodes <- c("EM", "FE", "HE", "JL", "SC", "TR")
#' mvad_seq <- seqdef(mvad, 17:86, alphabet = mvad_alphabet,
#'   states = mvad_scodes, labels = mvad_labels, xtstep = 6)
#'
#' # Define a color palette for the sequence data
#' attr(mvad_seq, "cpal") <- colorpalette[[6]]
#'
#' # Initialize a Markov model
#' init_mm_mvad <- build_mm(observations = mvad_seq,
#'   transition_probs = simulate_transition_probs(6),
#'   initial_probs = rep(1/6,6))
#'
#' # Fit the model
#' mm_mvad <- fit_model(init_mm_mvad)
#'
#' @seealso \code{\link{fit_model}} for estimating model parameters.

build_mm <- function(observations, transition_probs, initial_probs){

  if(!inherits(observations, "stslist")){
    stop("The build_mm function can only be used for single-channel sequence data (as an stslist object). Use the mc_to_sc_data function to convert multiple stslist into single-channel state sequences.")
  }

  n_sequences<-nrow(observations)
  length_of_sequences<-ncol(observations)
  state_names<-alphabet(observations)
  n_states<-length(state_names)

  if (!is.matrix(transition_probs)) {
    stop(paste("Object provided for transition_probs is not a matrix."))
  }
  if (!is.vector(initial_probs)) {
    stop(paste("Object provided for initial_probs is not a vector."))
  }

  if(dim(transition_probs)[1]!=dim(transition_probs)[2])
    stop("Transition_probs must be a square matrix.")

  if (n_states != nrow(transition_probs))
    stop("Number of rows in transition_probs is not equal to the number of distinct symbols in observations.")

  if (length(initial_probs) != n_states){
    stop(paste("Length of initial_probs is not equal to the number of states."))
  }

  if(!isTRUE(all.equal(rowSums(transition_probs),rep(1,dim(transition_probs)[1]),check.attributes=FALSE)))
    stop("Transition probabilities in transition_probs do not sum to one.")

  dimnames(transition_probs) <- list(from=state_names,to=state_names)



  names(initial_probs) <- state_names

    nobs <- sum(!(observations == attr(observations, "nr") |
        observations == attr(observations, "void") |
        is.na(observations)))

    emission_probs <- diag(n_states)
    dimnames(emission_probs) <- dimnames(transition_probs)

    model <- structure(list(observations=observations,transition_probs=transition_probs,
      emission_probs=emission_probs,initial_probs=initial_probs,
      state_names=state_names,
      symbol_names=state_names,channel_names=NULL,
      length_of_sequences=length_of_sequences,
      n_sequences=n_sequences,
      n_symbols=n_states,n_states=n_states,
      n_channels=1), class = "hmm",
      nobs = nobs,
      df = sum(initial_probs > 0) - 1 + sum(transition_probs > 0) - n_states,
      type = "mmm")

  model
}
