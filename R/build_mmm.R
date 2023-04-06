#' Build a Mixture Markov Model
#'
#' Function \code{build_mmm} is a shortcut for constructing a mixture Markov
#' model as a restricted case of an \code{mhmm} object.
#'
#' @export
#' @param observations An \code{stslist} object (see \code{\link[TraMineR]{seqdef}}) containing
#'   the sequences.
#' @param n_clusters A scalar giving the number of clusters/submodels
#' (not used if starting values for model parameters are given with
#' \code{initial_probs} and \code{transition_probs}).
#' @param transition_probs A list of matrices of transition
#'   probabilities for submodels of each cluster.
#' @param initial_probs A list which contains vectors of initial state
#'   probabilities for submodels of each cluster.
#' @param formula Optional formula of class \code{\link{formula}} for the
#' mixture probabilities. Left side omitted.
#' @param data A data frame containing the variables used in the formula.
#' Ignored if no formula is provided.
#' @param coefficients An optional \eqn{k x l} matrix of regression coefficients for
#'   time-constant covariates for mixture probabilities, where \eqn{l} is the number
#'   of clusters and \eqn{k} is the number of covariates. A logit-link is used for
#'   mixture probabilities. The first column is set to zero.
#' @param cluster_names A vector of optional names for the clusters.
#' @param ... Additional arguments to \code{simulate_transition_probs}.
#' @return Object of class \code{mhmm} with following elements:
#' \describe{
#'    \item{\code{observations}}{State sequence object or a list of such containing the data.}
#'    \item{\code{transition_probs}}{A matrix of transition probabilities.}
#'    \item{\code{emission_probs}}{A matrix or a list of matrices of emission probabilities.}
#'    \item{\code{initial_probs}}{A vector of initial probabilities.}
#'    \item{\code{coefficients}}{A matrix of parameter coefficients for covariates (covariates in rows, clusters in columns).}
#'    \item{\code{X}}{Covariate values for each subject.}
#'    \item{\code{cluster_names}}{Names for clusters.}
#'    \item{\code{state_names}}{Names for hidden states.}
#'    \item{\code{symbol_names}}{Names for observed states.}
#'    \item{\code{channel_names}}{Names for channels of sequence data}
#'    \item{\code{length_of_sequences}}{(Maximum) length of sequences.}
#'    \item{\code{n_sequences}}{Number of sequences.}
#'    \item{\code{n_symbols}}{Number of observed states (in each channel).}
#'    \item{\code{n_states}}{Number of hidden states.}
#'    \item{\code{n_channels}}{Number of channels.}
#'    \item{\code{n_covariates}}{Number of covariates.}
#'    \item{\code{n_clusters}}{Number of clusters.}
#' }
#' @seealso \code{\link{fit_model}} for estimating model parameters;
#' \code{\link{summary.mhmm}} for a summary of a mixture model;
#' \code{\link{separate_mhmm}} for organizing an \code{mhmm} object into a list of
#' separate \code{hmm} objects; and \code{\link{plot.mhmm}} for plotting
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
#' mvad_seq <- seqdef(mvad, 17:86,
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
  multichannel <- is_multichannel(observations)
  # Single channel but observations is a list
  if (is.list(observations) && !inherits(observations, "stslist") && length(observations) == 1) {
    observations <- observations[[1]]
    multichannel <- FALSE
  }
  if (multichannel) {
    stop(
      paste0("The 'build_mmm' function can only be used for single-channel ",
             "sequence data (as an stslist object). Use the 'mc_to_sc_data' function ",
             "to convert data into single-channel state sequences."
      )
    )
  }
  n_sequences <- nrow(observations)
  length_of_sequences <- ncol(observations)
  symbol_names <- alphabet(observations)
  n_symbols <- length(symbol_names)


  if (!missing(transition_probs) && !missing(initial_probs)) {
    if (is.list(transition_probs)) {
      n_clusters <- length(transition_probs)
    } else {
      stop("Argument 'transition_probs' is not a list.")
    }
    # Simulate starting values
  } else {
    if (missing(n_clusters)) {
      stop(paste0(
        "Provide either 'n_clusters' or both 'initial_probs' ",
        "and 'transition_probs'."
      ))
    }
    n_states <- rep(n_symbols, n_clusters)
    transition_probs <- simulate_transition_probs(
      n_states = n_states,
      n_clusters = n_clusters, ...)
    initial_probs <- simulate_initial_probs(
      n_states = n_states, n_clusters = n_clusters)
  }
  emission_probs <- replicate(n_clusters, diag(n_symbols), simplify = FALSE)
  state_names <- replicate(n_clusters, symbol_names, simplify = FALSE)

  model <- build_mhmm(
    observations = observations, transition_probs = transition_probs,
    emission_probs = emission_probs, initial_probs = initial_probs,
    formula = formula, data = data, coefficients = coefficients,
    cluster_names = cluster_names, state_names = state_names)
  attr(model, "type") <- "mmm"
  model
}
