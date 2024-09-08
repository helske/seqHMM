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
#' d <- reshape(mvad, direction = "long", varying = list(15:86), 
#'   v.names = "activity")
#' 
#' \dontrun{
#' set.seed(1)
#' fit <- estimate_mnhmm("activity", n_states = 3, n_clusters = 2,
#'   data = d, time = "time", id = "id", 
#'   cluster_formula = ~ male + catholic + gcse5eq + Grammar + 
#'     funemp + fmpr + livboth + Belfast +
#'   N.Eastern + Southern + S.Eastern + Western,
#'   initial_formula = ~ 1, emission_formula =  ~ male + catholic + gcse5eq,
#'   transition_formula = ~ male + gcse5eq, inits = "random"
#'   )
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
