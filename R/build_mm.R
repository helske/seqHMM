#' Build a Markov Model
#'
#' Function [build_mm()] builds and automatically estimates a Markov model. It is also a shortcut for
#' constructing a Markov model as a restricted case of an `hmm` object.
#' @export
#' @param observations An `stslist` object (see [TraMineR::seqdef()]) containing
#' the sequences.
#' @return Object of class `hmm` with following elements:
#' * `observations`\cr State sequence object or a list of such containing the data.
#' * `transition_probs`\cr A matrix of transition probabilities.
#' * `emission_probs`\cr A matrix or a list of matrices of emission probabilities.
#' * `initial_probs`\cr A vector of initial probabilities.
#' * `state_names`\cr Names for hidden states.
#' * `symbol_names`\cr Names for observed states.
#' * `channel_names`\cr Names for channels of sequence data
#' * `length_of_sequences`\cr (Maximum) length of sequences.
#' * `sequence_lengths`\cr A vector of sequence lengths.
#' * `n_sequences`\cr Number of sequences.
#' * `n_symbols`\cr Number of observed states (in each channel).
#' * `n_states`\cr Number of hidden states.
#' * `n_channels`\cr Number of channels.
#'
#' @details Unlike the other build functions in `seqHMM`, the [build_mm()] function
#' automatically estimates the model parameters. In case of no missing values,
#' initial and transition probabilities are
#' directly estimated from the observed initial state probabilities and transition counts.
#' In case of missing values, the EM algorithm is run once.
#'
#' Note that it is possible that the data contains a symbol from which there are
#' no transitions anywhere (even to itself), which would lead to a row in
#' transition matrix full of zeros. In this case the [build_mm()]
#' (as well as the EM algorithm) assumes that the
#' the state is absorbing in a way that probability of staying in this state is 1.
#'
#' @seealso [plot.hmm()] for plotting the model.
#'
#' @examples
#' # Construct sequence data
#' data("mvad", package = "TraMineR")
#'
#' mvad_alphabet <-
#'   c("employment", "FE", "HE", "joblessness", "school", "training")
#' mvad_labels <- c(
#'   "employment", "further education", "higher education",
#'   "joblessness", "school", "training"
#' )
#' mvad_scodes <- c("EM", "FE", "HE", "JL", "SC", "TR")
#' mvad_seq <- seqdef(mvad, 15:86,
#'   alphabet = mvad_alphabet,
#'   states = mvad_scodes, labels = mvad_labels, xtstep = 6,
#'   cpal = colorpalette[[6]]
#' )
#'
#' # Estimate the Markov model
#' mm_mvad <- build_mm(observations = mvad_seq)
#'
build_mm <- function(observations) {
  
  observations <- .check_observations(observations)
  n_channels <- attr(observations, "n_channels")
  stopifnot_(
    n_channels == 1,
    paste0(
      "{.fn build_mm} can only be used for single-channel sequence data ",
      "(a {.cls stslist} object). Use {.fn mc_to_sc_data} to convert data ",
      "into single-channel state sequences."
    )
  )
  state_names <- alphabet(observations)
  n_states <- length(state_names)
  
  if (any(observations == attr(observations, "nr"))) {
    model <- build_hmm(
      observations,
      transition_probs = matrix(1 / n_states, n_states, n_states),
      emission_probs = diag(n_states),
      initial_probs = rep(1 / n_states, n_states),
      state_names = state_names)
    model <- fit_model(model)$model
    message_("Sequences contain missing values, initial and transition 
             probabilities estimated via EM. ")
  } else {
    first_timepoint <- suppressMessages(
      seqdef(
        observations[observations[, 1] %in% state_names, 1], 
        alphabet = state_names
      )
    )
    initial_probs <- TraMineR::seqstatf(first_timepoint)[, 2] / 100
    transition_probs <- suppressMessages(TraMineR::seqtrate(observations))
    zeros <- which(rowSums(transition_probs) == 0)
    diag(transition_probs)[zeros] <- 1
    if (length(zeros) > 0) {
      warning_("There are no observed transitions from some of the symbols.")
    }
    model <- build_hmm(
      observations,
      transition_probs = transition_probs,
      emission_probs = diag(n_states),
      initial_probs = initial_probs,
      state_names = state_names)
  }
  model$call <- match.call()
  attr(model, "type") <- "mm"
  model
}
