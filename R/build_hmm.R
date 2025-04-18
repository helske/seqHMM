#' Build a Hidden Markov Model
#'
#' Function `build_hmm` constructs a hidden Markov model object of class `hmm`.
#'
#' The returned model contains some attributes such as `nobs` and `df`,
#' which define the number of observations in the  model and the number of estimable
#' model parameters, used in computing BIC.
#' When computing `nobs` for a multichannel model with \eqn{C} channels,
#' each observed value in a single channel amounts to \eqn{1/C} observation,
#' i.e. a fully observed time point for a single sequence amounts to one observation.
#' For the degrees of freedom `df`, zero probabilities of the initial model are
#' defined as structural zeroes.
#' @export
#' @param observations An `stslist` object (see [TraMineR::seqdef()]) containing
#' the sequences, or a list of such objects (one for each channel).
#' @param n_states A scalar giving the number of hidden states. Not used if starting values for model parameters
#' are given with `initial_probs`, `transition_probs`, or `emission_probs`.
#' @param transition_probs A matrix of transition probabilities.
#' @param emission_probs A matrix of emission probabilities or a list of such
#' objects (one for each channel). Emission probabilities should follow the
#' ordering of the alphabet of observations (`alphabet(observations)`, returned as `symbol_names`).
#' @param initial_probs A vector of initial state probabilities.
#' @param state_names A list of optional labels for the hidden states. If `NULL`,
#' the state names are taken from the row names of the transition matrix. If this is
#' also `NULL`, numbered states are used.
#' @param channel_names A vector of optional names for the channels.
#' @param ... Additional arguments to [simulate_transition_probs()].
#' @return Object of class `hmm` with the following elements:
#' * `observations`\cr State sequence object or a list of such objects containing the data.
#' * `transition_probs`\cr A matrix of transition probabilities.
#' * `emission_probs`\cr A matrix or a list of matrices of emission probabilities.
#' * `initial_probs`\cr A vector of initial probabilities.
#' * `state_names`\cr Names for hidden states.
#' * `symbol_names`\cr Names for observed states.
#' * `channel_names`\cr Names for channels of sequence data.
#' * `length_of_sequences`\cr (Maximum) length of sequences.
#' * `sequence_lengths`\cr A vector of sequence lengths.
#' * `n_sequences`\cr Number of sequences.
#' * `n_symbols`\cr Number of observed states (in each channel).
#' * `n_states`\cr Number of hidden states.
#' * `n_channels`\cr Number of channels.
#'
#' @seealso [fit_model()] for estimating model parameters; and
#'   [plot.hmm()] for plotting `hmm` objects.
#' @examples
#'
#' # Single-channel data
#'
#' data("mvad", package = "TraMineR")
#'
#' mvad_alphabet <- c(
#'   "employment", "FE", "HE", "joblessness", "school",
#'   "training"
#' )
#' mvad_labels <- c(
#'   "employment", "further education", "higher education",
#'   "joblessness", "school", "training"
#' )
#' mvad_scodes <- c("EM", "FE", "HE", "JL", "SC", "TR")
#' mvad_seq <- seqdef(mvad, 15:86,
#'   alphabet = mvad_alphabet, states = mvad_scodes,
#'   labels = mvad_labels, xtstep = 6
#' )
#'
#' # Initializing an HMM with 4 hidden states, random starting values
#' init_hmm_mvad1 <- build_hmm(observations = mvad_seq, n_states = 4)
#'
#' # Starting values for the emission matrix
#' emiss <- matrix(NA, nrow = 4, ncol = 6)
#' emiss[1, ] <- seqstatf(mvad_seq[, 1:12])[, 2] + 1
#' emiss[2, ] <- seqstatf(mvad_seq[, 13:24])[, 2] + 1
#' emiss[3, ] <- seqstatf(mvad_seq[, 25:48])[, 2] + 1
#' emiss[4, ] <- seqstatf(mvad_seq[, 49:70])[, 2] + 1
#' emiss <- emiss / rowSums(emiss)
#'
#' # Starting values for the transition matrix
#'
#' tr <- matrix(
#'   c(
#'     0.80, 0.10, 0.05, 0.05,
#'     0.05, 0.80, 0.10, 0.05,
#'     0.05, 0.05, 0.80, 0.10,
#'     0.05, 0.05, 0.10, 0.80
#'   ),
#'   nrow = 4, ncol = 4, byrow = TRUE
#' )
#'
#' # Starting values for initial state probabilities
#' init <- c(0.3, 0.3, 0.2, 0.2)
#'
#' # HMM with own starting values
#' init_hmm_mvad2 <- build_hmm(
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
#' marr_seq <- seqdef(biofam3c$married,
#'   start = 15,
#'   alphabet = c("single", "married", "divorced"),
#'   cpal = c("violetred2", "darkgoldenrod2", "darkmagenta")
#' )
#' child_seq <- seqdef(biofam3c$children,
#'   start = 15,
#'   alphabet = c("childless", "children"),
#'   cpal = c("darkseagreen1", "coral3")
#' )
#' left_seq <- seqdef(biofam3c$left,
#'   start = 15,
#'   alphabet = c("with parents", "left home"),
#'   cpal = c("lightblue", "red3")
#' )
#'
#' # You could also define the colors using cpal function from TraMineR
#' # cpal(marr_seq) <- c("violetred2", "darkgoldenrod2", "darkmagenta")
#' # cpal(child_seq) <- c("darkseagreen1", "coral3")
#' # cpal(left_seq) <- c("lightblue", "red3")
#'
#' # Left-to-right HMM with 3 hidden states and random starting values
#' set.seed(1010)
#' init_hmm_bf1 <- build_hmm(
#'   observations = list(marr_seq, child_seq, left_seq),
#'   n_states = 3, left_right = TRUE, diag_c = 2
#' )
#'
#'
#' # Starting values for emission matrices
#'
#' emiss_marr <- matrix(NA, nrow = 3, ncol = 3)
#' emiss_marr[1, ] <- seqstatf(marr_seq[, 1:5])[, 2] + 1
#' emiss_marr[2, ] <- seqstatf(marr_seq[, 6:10])[, 2] + 1
#' emiss_marr[3, ] <- seqstatf(marr_seq[, 11:16])[, 2] + 1
#' emiss_marr <- emiss_marr / rowSums(emiss_marr)
#'
#' emiss_child <- matrix(NA, nrow = 3, ncol = 2)
#' emiss_child[1, ] <- seqstatf(child_seq[, 1:5])[, 2] + 1
#' emiss_child[2, ] <- seqstatf(child_seq[, 6:10])[, 2] + 1
#' emiss_child[3, ] <- seqstatf(child_seq[, 11:16])[, 2] + 1
#' emiss_child <- emiss_child / rowSums(emiss_child)
#'
#' emiss_left <- matrix(NA, nrow = 3, ncol = 2)
#' emiss_left[1, ] <- seqstatf(left_seq[, 1:5])[, 2] + 1
#' emiss_left[2, ] <- seqstatf(left_seq[, 6:10])[, 2] + 1
#' emiss_left[3, ] <- seqstatf(left_seq[, 11:16])[, 2] + 1
#' emiss_left <- emiss_left / rowSums(emiss_left)
#'
#' # Starting values for transition matrix
#' trans <- matrix(
#'   c(
#'     0.9, 0.07, 0.03,
#'     0, 0.9, 0.1,
#'     0, 0, 1
#'   ),
#'   nrow = 3, ncol = 3, byrow = TRUE
#' )
#'
#' # Starting values for initial state probabilities
#' inits <- c(0.9, 0.09, 0.01)
#'
#' # HMM with own starting values
#' init_hmm_bf2 <- build_hmm(
#'   observations = list(marr_seq, child_seq, left_seq),
#'   transition_probs = trans,
#'   emission_probs = list(emiss_marr, emiss_child, emiss_left),
#'   initial_probs = inits
#' )
#'
build_hmm <- function(observations, n_states, transition_probs, emission_probs, initial_probs,
                      state_names = NULL, channel_names = NULL, ...) {
  
  observations <- .check_observations(observations, channel_names)
  n_channels <- attr(observations, "n_channels")
  n_symbols <- attr(observations, "n_symbols")
  channel_names <- attr(observations, "channel_names")
  symbol_names <- attr(observations, "symbol_names")
  
  inits_given <- !missing(transition_probs) || !missing(initial_probs) || !missing(emission_probs)
  n_states_given <- !missing(n_states)
  stopifnot_(
    inits_given || n_states_given,
    "Provide either {.arg n_states} or all three of {.arg initial_probs}, ",
    "{.arg transition_probs}, and {.arg emission_probs}."
    )
  if (inits_given) {
    transition_probs <- .check_transition_probs(transition_probs, state_names)
    n_states <- nrow(transition_probs)
    state_names <- rownames(transition_probs)
    initial_probs <- .check_initial_probs(initial_probs, n_states, state_names)
    emission_probs <- .check_emission_probs(
      emission_probs, n_states, n_channels, n_symbols, state_names, 
      symbol_names, channel_names
    )
  } else {
    # Simulate starting values
    transition_probs <- simulate_transition_probs(
      n_states = n_states, 
      n_clusters = 1, ...
    )
    initial_probs <- simulate_initial_probs(n_states = n_states, n_clusters = 1)
    emission_probs <- simulate_emission_probs(
      n_states = n_states, n_symbols = n_symbols, n_clusters = 1
    )
    transition_probs <- .check_transition_probs(transition_probs, state_names)
    state_names <- rownames(transition_probs)
    initial_probs <- .check_initial_probs(initial_probs, n_states, state_names)
    emission_probs <- .check_emission_probs(
      emission_probs, n_states, n_channels, n_symbols, state_names, 
      symbol_names, channel_names
    )
  }
  if (attr(observations, "nobs") == 0) {
    warning_("Sequences contain only missing values.")
  }
  model <- structure(
    list(
      observations = observations, transition_probs = transition_probs,
      emission_probs = emission_probs, initial_probs = initial_probs,
      state_names = as_factor(state_names),
      symbol_names = symbol_names, 
      channel_names = channel_names,
      length_of_sequences = attr(observations, "length_of_sequences"),
      sequence_lengths = attr(observations, "sequence_lengths"),
      n_sequences = attr(observations, "n_sequences"),
      n_symbols = n_symbols, n_states = n_states,
      n_channels = n_channels,
      call = match.call()
    ),
    class = "hmm",
    nobs = attr(observations, "nobs"),
    df = sum(initial_probs > 0) - 1 + sum(transition_probs > 0) - n_states +
      sum(unlist(emission_probs) > 0) - n_states * n_channels,
    type = "hmm"
  )
  
  model
}
