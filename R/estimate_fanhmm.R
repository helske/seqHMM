#' Build and Estimate a Feedback-Augmented Non-homogeneous Hidden Markov Model
#'
#' Function `estimate_fanhmm` estimates a hidden Markov model object of class 
#' `fanhmm` where transition and/or emission probabilities depend on past 
#' responses.
#' 
#' @inheritParams estimate_nhmm
#' @param autoregression_formula Formula for autoregression 
#' \eqn{y_t \to y_{t+1}}.
#' Default intercept means is shorthand for \eqn{y_{t+1} ~ y_{t}}, while 
#' additional terms in formula are interacted with the lagged responses. 
#' If `NULL`, no autoregression is used.
#' @param feedback_formula Formula for feedback \eqn{y_t \to z_{t+1}}.
#' Default intercept means is shorthand for \eqn{z_{t+1} ~ y_{t}}, while 
#' additional terms in formula are interacted with the lagged responses. 
#' If `NULL`, no feedback is used.
#' @param inits If `inits = "random"` (default), random initial values are 
#' used. Otherwise `inits` should be list of initial values. If coefficients 
#' are given using list components `eta_pi`, `eta_A`, `eta_B`, 
#' and `eta_omega`, these are used as is, alternatively initial values 
#' can be given in terms of the initial state, transition, emission, and mixture 
#' probabilities using list components `initial_probs`, `emission_probs`, 
#' `transition_probs`, and `cluster_probs`. These can also be mixed, i.e. you 
#' can give only `initial_probs` and `eta_A`.
#' @return Object of class `fanhmm`.
#' @seealso [estimate_nhmm()] for further details.
#' @export
estimate_fanhmm <- function(
    responses, n_states, initial_formula = ~1, 
    transition_formula = ~1, emission_formula = ~1, autoregression_formula = ~1,
    feedback_formula = ~1,
    data, time, id, state_names = NULL,
    inits = "random", init_sd = 2, 
    restarts = 0L, lambda = 0, method = "EM-DNM", bound = Inf, 
    control_restart = list(), control_mstep = list(), ...) {
  
  call <- match.call()
  model <- build_fanhmm(
    responses, n_states, initial_formula, 
    transition_formula, emission_formula, autoregression_formula, 
    feedback_formula, data, id, time, state_names
  )
  control <- list(...)
  start_time <- proc.time()
  out <- fit_nhmm(model, inits, init_sd, restarts, lambda, method, bound, 
                   control, control_restart, control_mstep)
  end_time <- proc.time()
  out$estimation_results$time <- end_time - start_time
  attr(out, "call") <- call
  out
}
