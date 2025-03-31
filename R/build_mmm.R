#' Build a Mixture Markov Model
#'
#' Function [build_mmm()] is a shortcut for constructing a mixture Markov
#' model as a restricted case of an `mhmm` object.
#'
#' @export
#' @param observations An `stslist` object (see [TraMineR::seqdef()]) containing
#'   the sequences.
#' @param n_clusters A scalar giving the number of clusters/submodels
#' (not used if starting values for model parameters are given with
#' `initial_probs` and `transition_probs`).
#' @param transition_probs A list of matrices of transition
#'   probabilities for submodels of each cluster.
#' @param initial_probs A list which contains vectors of initial state
#'   probabilities for submodels of each cluster.
#' @param formula Optional formula of class [formula()] for the
#' mixture probabilities. Left side omitted.
#' @param data A data frame containing the variables used in the formula.
#' Ignored if no formula is provided.
#' @param coefficients An optional \eqn{k x l} matrix of regression coefficients for
#'   time-constant covariates for mixture probabilities, where \eqn{l} is the number
#'   of clusters and \eqn{k} is the number of covariates. A logit-link is used for
#'   mixture probabilities. The first column is set to zero.
#' @param cluster_names A vector of optional names for the clusters.
#' @param ... Additional arguments to `simulate_transition_probs`.
#' @return Object of class `mhmm` with following elements:
#' * `observations`\cr State sequence object or a list of such containing the data.
#' * `transition_probs`\cr A matrix of transition probabilities.
#' * `emission_probs`\cr A matrix or a list of matrices of emission probabilities.
#' * `initial_probs`\cr A vector of initial probabilities.
#' * `coefficients`\cr A matrix of parameter coefficients for covariates (covariates in rows, clusters in columns).
#' * `X`\cr Covariate values for each subject.
#' * `cluster_names`\cr Names for clusters.
#' * `state_names`\cr Names for hidden states.
#' * `symbol_names`\cr Names for observed states.
#' * `channel_names`\cr Names for channels of sequence data
#' * `length_of_sequences`\cr (Maximum) length of sequences.
#' * `sequence_lengths`\cr A vector of sequence lengths.
#' * `n_sequences`\cr Number of sequences.
#' * `n_symbols`\cr Number of observed states (in each channel).
#' * `n_states`\cr Number of hidden states.
#' * `n_channels`\cr Number of channels.
#' * `n_covariates`\cr Number of covariates.
#' * `n_clusters`\cr Number of clusters.
#' 
#' @seealso [fit_model()] for estimating model parameters;
#' [summary.mhmm()] for a summary of a mixture model;
#' [separate_mhmm()] for organizing an `mhmm` object into a list of
#' separate `hmm` objects; and [plot.mhmm()] for plotting
#' mixture models.
#'
#' @examples
#'
#'
#' # Define sequence data
#' data("mvad", package = "TraMineR")
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
#' # Initialize the MMM
#' set.seed(123)
#' mmm_mvad <- build_mmm(
#'   observations = mvad_seq,
#'   n_clusters = 2,
#'   formula = ~male, data = mvad
#' )
#'
#' \dontrun{
#' # Estimate model parameters
#' mmm_mvad <- fit_model(mmm_mvad)$model
#'
#' # Plot model (both clusters in the same plot)
#' require(igraph)
#' plot(mmm_mvad,
#'   interactive = FALSE,
#'   # Modify legend position and properties
#'   with.legend = "right", legend.prop = 0.3, cex.legend = 1.2,
#'   # Define vertex layout
#'   layout = layout_as_star,
#'   # Modify edge properties
#'   edge.label = NA, edge.arrow.size = 0.8, edge.curved = 0.2,
#'   # Modify vertex label positions (initial probabilities)
#'   vertex.label.pos = c("left", "right", "right", "left", "left", "right")
#' )
#'
#' # Summary of the MMM
#' summary(mmm_mvad)
#' }
build_mmm <- function(observations, n_clusters, transition_probs, initial_probs,
                      formula = NULL, data = NULL, coefficients = NULL,
                      cluster_names = NULL, ...) {
  observations <- .check_observations(observations)
  n_channels <- attr(observations, "n_channels")
  stopifnot_(
    n_channels == 1,
    paste0(
      "{.fn build_mmm} can only be used for single-channel sequence data ",
      "(a {.cls stslist} object). Use {.fn mc_to_sc_data} to convert data ",
      "into single-channel state sequences."
    )
  )
  n_symbols <- n_states <- attr(observations, "n_symbols")
  symbol_names <- state_names <- attr(observations, "symbol_names")
  emission_probs <- diag(n_symbols)
  rownames(emission_probs) <- colnames(emission_probs) <- symbol_names
  if (!missing(transition_probs) && !missing(initial_probs)) {
    stopifnot_(
      is.list(transition_probs),
      "{.arg transition_probs} is not a {.cls list}."
    )
    stopifnot_(
      is.list(initial_probs),
      "{.arg initial_probs} is not a {.cls list}."
    )
    n_clusters <- length(transition_probs)
    for (i in seq_len(n_clusters)) {
      transition_probs[[i]] <- .check_transition_probs(
        transition_probs[[i]],
        state_names)
      initial_probs[[i]] <- .check_initial_probs(
        initial_probs[[i]], 
        n_states,
        state_names
      )
    }
    # Simulate starting values
  } else {
    if (missing(n_clusters)) {
      stop_(
        "Provide either {.arg n_clusters} or both {.arg initial_probs} and 
        {.arg transition_probs}."
      )
    }
    n_states <- rep(n_symbols, n_clusters)
    transition_probs <- simulate_transition_probs(
      n_states = n_states,
      n_clusters = n_clusters, ...)
    initial_probs <- simulate_initial_probs(
      n_states = n_states, n_clusters = n_clusters)
  }
  emission_probs <- replicate(n_clusters, emission_probs, simplify = FALSE)
  state_names <- replicate(n_clusters, state_names, simplify = FALSE)
  names(transition_probs) <- names(emission_probs) <- names(initial_probs) <- 
    cluster_names
  model <- build_mhmm(
    observations = observations, transition_probs = transition_probs,
    emission_probs = emission_probs, initial_probs = initial_probs,
    formula = formula, data = data, coefficients = coefficients,
    cluster_names = cluster_names, state_names = state_names)
  model$call <- match.call()
  attr(model, "type") <- "mmm"
  model
}
