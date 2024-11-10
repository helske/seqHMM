#' Build and Estimate a Non-homogeneous Hidden Markov Model
#'
#' Function `estimate_nhmm` estimates a hidden Markov model object of class 
#' `nhmm` where initial, transition and emission probabilities 
#' (potentially) depend on covariates.
#' 
#' By default, the model parameters are estimated using LBFGS algorithm of 
#' [nloptr::nloptr()]. The log-likelihood is 
#' scaled by the number of non-missing observations (`nobs(model)`), and the 
#' convergence is claimed when either the absolute or relative change of this 
#' objective function is less than `1e-8`, or the absolute change of the 
#' parameters is less than `1e-8`. The covariate data is standardardized before 
#' optimization.
#' 
#' With multiple runs of optimization (by using the `restarts` argument), it is 
#' possible to parallelize these runs using the `future` package, e.g., by 
#' calling `future::plan(multisession, workers = 2)` before `estimate_nhmm()`.
#' See [future::plan()] for details.
#' 
#' @param observations Either the name of the response variable in `data`, or 
#' an `stslist` object (see [TraMineR::seqdef()]) containing the 
#' sequences. In case of multichannel data, `observations` should be a vector 
#' of response variable names in `data`, or a list of `stslist` objects.
#' @param n_states A positive integer defining the number of hidden states.
#' @param initial_formula of class [formula()] for the
#' initial state probabilities.
#' @param transition_formula of class [formula()] for the
#' state transition probabilities.
#' @param emission_formula of class [formula()] for the
#' state emission probabilities.
#' @param data A data frame containing the variables used in the model 
#' formulas. Can be omitted in case of model with no covariates and observations 
#' given as `stslist` objects.
#' @param time Name of the time index variable in `data`.
#' @param id Name of the id variable in `data` identifying different 
#' sequences.
#' @param state_names A vector of optional labels for the hidden states. If this
#' is `NULL` (the default), numbered states are used.
#' @param channel_names A vector of optional names for the channels. If this
#' is `NULL` (the default), numbered channels are used.
#' @param inits If `inits = "random"` (default), random initial values are 
#' used. Otherwise `inits` should be list of initial values. If coefficients 
#' are given using list components `eta_pi`, `eta_A`, `eta_B`, 
#' these are used as is, alternatively initial values can be given in terms of 
#' the initial state, transition, and emission probabilities using list 
#' components `initial_probs`, `emission_probs`, and `transition_probs`. These 
#' can also be mixed, i.e. you can give only `initial_probs` and `eta_A`.
#' @param init_sd Standard deviation of the normal distribution used to generate
#' random initial values. Default is `2`. If you want to fix the initial values 
#' of the regression coefficients to zero, use `init_sd = 0`.
#' @param restarts Number of times to run optimization using random starting 
#' values (in addition to the final run). Default is 0.
#' @param lambda Penalization factor `lambda` for penalized log-likelihood, where the 
#' penalization is `lambda * sum(parameters^2)/(2 * n_obs)`, where `n_obs` is 
#' the number of non-missing observations.
#' @param method Optimization method used. Default is `"EM"` which uses EM
#' algorithm with L-BFGS in the M-step. Another option is `"DNM"` which uses 
#' direct maximization of the log-likelihood using [nloptr::nloptr()].
#' @param pseudocount A positive scalar to be added for the expected counts of 
#' E-step. Only used in EM algorithm. Default is 1e-4. Larger values can be used 
#' to avoid zero probabilities in initial, transition, and emission 
#' probabilities, i.e. these have similar role as `lambda`.
#' @param store_data If `TRUE` (default), original data frame passed as `data` 
#' is stored to the model object. For large datasets, this can be set to 
#' `FALSE`, in which case you might need to pass the data separately to some 
#' post-prosessing functions.
#' @param ... Additional arguments to [nloptr::nloptr()]. Most importantly,
#' argument `maxeval` defines the maximum number of iterations for optimization.
#' The default is `1000` for restarts and `10000` for the final optimization. 
#' Other useful arguments are `algorithm` (default uses LBFGS),
#' `print_level` (default is `0`, no console output of optimization), and 
#' arguments for adjusting the stopping criteria of the optimization 
#' (see details).
#' @return Object of class `nhmm`.
#' @export
#' @examples
#' data("mvad", package = "TraMineR")
#' 
#' d <- reshape(mvad, direction = "long", varying = list(15:86), 
#'   v.names = "activity")
#' 
#' \dontrun{
#' set.seed(1)
#' fit <- estimate_nhmm("activity", n_states = 3,
#'   data = d, time = "time", id = "id", 
#'   initial_formula = ~ 1, emission_formula =  ~ male + gcse5eq,
#'   transition_formula = ~ male + gcse5eq, inits = "random"
#'   )
#' }
estimate_nhmm <- function(
    observations, n_states, initial_formula = ~1, 
    transition_formula = ~1, emission_formula = ~1, 
    data = NULL, time = NULL, id = NULL, state_names = NULL, channel_names = NULL, 
    inits = "random", init_sd = 2, restarts = 0L, lambda = 0, method = "EM",
    pseudocount = 1e-4, store_data = TRUE, ...) {
  
  call <- match.call()
  
  method <- match.arg(method, c("DNM", "EM"))
  
  model <- build_nhmm(
    observations, n_states, initial_formula, 
    transition_formula, emission_formula, data, time, id, state_names, 
    channel_names
    )
  stopifnot_(
    checkmate::test_flag(x = store_data), 
    "Argument {.arg store_data} must be a single {.cls logical} value."
  )
  stopifnot_(
    checkmate::check_number(lambda, lower = 0), 
    "Argument {.arg lambda} must be a single non-negative {.cls numeric} value."
  )
  if (store_data) {
    model$data <- data
  }
  out <- fit_nhmm(model, inits, init_sd, restarts, lambda, method, pseudocount, 
                  ...)
  attr(out, "call") <- call
  out
}
