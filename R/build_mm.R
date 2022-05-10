#' Build a Markov Model
#'
#' Function \code{build_mm} builds and automatically estimates a Markov model. It is also a shortcut for 
#' constructing a Markov model as a restricted case of an \code{hmm} object.
#' @export
#' @param observations An \code{stslist} object (see \code{\link[TraMineR]{seqdef}}) containing
#' the sequences.
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
#' @details Unlike the other build functions in \code{seqHMM}, the \code{build_mm} function
#' automatically estimates the model parameters. In case of no missing values, 
#' initial and transition probabilities are
#' directly estimated from the observed initial state probabilities and transition counts.
#' In ase of missing values, the EM algorithm is run once.
#'
#' @seealso \code{\link{plot.hmm}} for plotting the model.
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
#' # Estimate the Markov model
#' mm_mvad <- build_mm(observations = mvad_seq)
#'

build_mm <- function(observations){
  
  if(!inherits(observations, "stslist")){
    stop("The build_mm function can only be used for single-channel sequence data (as an stslist object). Use the mc_to_sc_data function to convert multiple stslist into single-channel state sequences.")
  }
  
  n_sequences <- nrow(observations)
  length_of_sequences <- ncol(observations)
  state_names <- alphabet(observations)
  n_states <- length(state_names)
  nobs <- sum(!(observations == attr(observations, "nr") |
                  observations == attr(observations, "void") |
                  is.na(observations)))
  

  if (nobs < prod(dim(observations))) {
    model <- fit_model(build_hmm(observations, 
                                            transition_probs = matrix(1/n_states, n_states, n_states), 
                                            emission_probs = diag(n_states),
                                            initial_probs = rep(1 / n_states, n_states)))$model
    warning("Sequences contain missing/void values, initial and transition probabilities estimated via EM. ")
    initial_probs <- model$initial_probs
    transition_probs <- model$transition_probs
  } else {
    first_timepoint <- suppressMessages(seqdef(observations[observations[, 1] %in% state_names, 1], alphabet = state_names))
    initial_probs <- TraMineR::seqstatf(first_timepoint)[, 2] / 100
    transition_probs <- suppressMessages(TraMineR::seqtrate(observations))
  }
  names(initial_probs) <- state_names
  dimnames(transition_probs) <- list(from = state_names, to = state_names)
  
  emission_probs <- diag(n_states)
  dimnames(emission_probs) <- dimnames(transition_probs)
  
  model <- structure(list(observations=observations,
                          transition_probs=transition_probs,
                          emission_probs=emission_probs,initial_probs=initial_probs,
                          state_names=state_names,
                          symbol_names=state_names, channel_names=NULL,
                          length_of_sequences=length_of_sequences,
                          n_sequences=n_sequences,
                          n_symbols=n_states,n_states=n_states,
                          n_channels=1), class = "hmm",
                     nobs = nobs,
                     df = sum(initial_probs > 0) - 1 + sum(transition_probs > 0) - n_states,
                     type = "mm")
  
  model
}
