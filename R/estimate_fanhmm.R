#' Build and Estimate a Feedback-Augmented Non-homogeneous Hidden Markov Model
#'
#' Function `estimate_fanhmm` estimates a hidden Markov model object of class 
#' `fanhmm` where transition and/or emission probabilities depend on past 
#' responses.
#' 
#' @inheritParams estimate_nhmm
#' @param autoregression_formula Formula for autoregression \eqn{y_t \to y_{t+1}}.
#' @param feedback_formula Formula for feedback \eqn{y_t \to z_{t+1}}.
#' @param inits If `inits = "random"` (default), random initial values are 
#' used. Otherwise `inits` should be list of initial values. If coefficients 
#' are given using list components `eta_pi`, `eta_A`, `eta_B`, 
#' and `eta_omega`, these are used as is, alternatively initial values 
#' can be given in terms of the initial state, transition, emission, and mixture 
#' probabilities using list components `initial_probs`, `emission_probs`, 
#' `transition_probs`, and `cluster_probs`. These can also be mixed, i.e. you 
#' can give only `initial_probs` and `eta_A`.
#' @param cluster_names A vector of optional labels for the clusters. If this
#' is `NULL` (the default), numbered clusters are used.
#' @return Object of class `mnhmm`.
#' @seealso [estimate_nhmm()] for further details.
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
estimate_fanhmm <- function(
    observations, n_states, initial_formula = ~1, 
    transition_formula = ~1, emission_formula = ~1, autoregression_formula = ~1,
    feedback_formula = ~1,
    data = NULL, time = NULL, id = NULL, state_names = NULL, 
    inits = "random", init_sd = 2, 
    restarts = 0L, lambda = 0, method = "EM-DNM", bound = Inf, 
    control_restart = list(), control_mstep = list(), store_data = TRUE, ...) {
  
  call <- match.call()
  model <- build_fanhmm(
    observations, n_states, initial_formula, 
    transition_formula, emission_formula, autoregression_formula, 
    feedback_formula, data, time, id, state_names
  )
  stopifnot_(
    checkmate::test_flag(x = store_data), 
    "Argument {.arg store_data} must be a single {.cls logical} value."
  )
  if (store_data) {
    model$data <- data
  }
  control <- list(...)
  start_time <- proc.time()
  out <- fit_fanhmm(model, inits, init_sd, restarts, lambda, method, bound, 
                   control, control_restart, control_mstep)
  end_time <- proc.time()
  out$estimation_results$time <- end_time - start_time
  attr(out, "call") <- call
  out
}
