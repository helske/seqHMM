#' Build and Estimate a Mixture Non-homogeneous Hidden Markov Model
#'
#' Function `estimate_mnhmm` estimates a hidden Markov model object of class 
#' `mnhmm` where initial, transition, emission, and mixture probabilities 
#' (potentially) depend on covariates.
#' 
#' @inheritParams estimate_nhmm
#' @param n_clusters A positive integer defining the number of clusters 
#' (mixtures).
#' @param cluster_formula of class [formula()] for the
#' mixture probabilities.
#'  @param inits If `inits = "random"` (default), random initial values are 
#' used. Otherwise `inits` should be list of initial values. If coefficients 
#' are given using list components `beta_i_raw`, `beta_s_raw`, `beta_o_raw`, 
#' and `theta_raw`, these are used as is, alternatively initial values can be 
#' given in terms of the initial state, transition, emission, and mixture 
#' probabilities using list components `initial_probs`, `emission_probs`, 
#' `transition_probs`, and `cluster_probs`. These can also be mixed, i.e. you 
#' can give only `initial_probs` and `beta_s_raw`.
#' @param cluster_names A vector of optional labels for the clusters. If this
#' is `NULL` (the default), numbered clusters are used.
#' @return Object of class `mnhmm`.
#' @export
#' @examples
#' data("mvad", package = "TraMineR")
#' 
#' mvad_alphabet <-
#'   c("employment", "FE", "HE", "joblessness", "school", "training")
#' mvad_labels <- c("employment", "further education", "higher education",
#'                  "joblessness", "school", "training")
#' mvad_scodes <- c("EM", "FE", "HE", "JL", "SC", "TR")
#' mvad_seq <- seqdef(mvad, 17:86, alphabet = mvad_alphabet,
#'                    states = mvad_scodes, labels = mvad_labels, xtstep = 6,
#'                    cpal = unname(colorpalette[[6]]))
#' 
#' \dontrun{
#' set.seed(1)
#' fit <- estimate_mnhmm(mvad_seq, n_states = 3, n_clusters = 2)
#' }
estimate_mnhmm <- function(
    observations, n_states, n_clusters, initial_formula = ~1, 
    transition_formula = ~1, emission_formula = ~1, cluster_formula = ~1,
    data = NULL, time = NULL, id = NULL, state_names = NULL, 
    channel_names = NULL, cluster_names = NULL, inits = "random", init_sd = 2, 
    restarts = 1L, threads = 1L, store_data = TRUE, verbose = TRUE, ...) {
  
  model <- build_mnhmm(
    observations, n_states, n_clusters, initial_formula, 
    transition_formula, emission_formula, cluster_formula, data, time, id, 
    state_names, channel_names, cluster_names
    )
  stopifnot_(
    checkmate::test_flag(x = store_data), 
    "Argument {.arg store_data} must be a single {.cls logical} value.")
  if (store_data) {
    model$data <- data
  }
  fit_mnhmm(model, inits, init_sd, restarts, threads, verbose, ...)
}
