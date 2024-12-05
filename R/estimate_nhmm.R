#' Build and Estimate a Non-homogeneous Hidden Markov Model
#'
#' Function `estimate_nhmm` estimates a hidden Markov model object of class 
#' `nhmm` where initial, transition and emission probabilities 
#' (potentially) depend on covariates.
#' 
#' By default, the model parameters are estimated using EM-DNM algorithm 
#' which first runs some iterations (100 by default) of EM algorithm, and then 
#' switches to L-BFGS. Other options include any numerical optimization 
#' algorithm of [nloptr::nloptr()], or plain EM algorithm where the 
#' M-step uses L-BFGS (provided by the NLopt library).
#' 
#' With multiple runs of optimization (by using the `restarts` argument), it is 
#' possible to parallelize these runs using the `future` package, e.g., by 
#' calling `future::plan(multisession, workers = 2)` before `estimate_nhmm()`.
#' See [future::plan()] for details. This is compatible with `progressr` 
#' package, so you can use [progressr::with_progress()] to track the progress 
#' of these multiple runs.
#' 
#' During the estimation, the log-likelihood is scaled by the number of 
#' non-missing observations (`nobs(model)`), and the the covariate data is 
#' standardardized before optimization.
#'  
#' By default, the convergence is claimed when the relative 
#' change of the objective function is less than `1e-8`, or the 
#' relative change of the parameters is less than `1e-6`. These can be changed
#' by passing arguments `ftol_rel` and `xtol_rel` via `...`. These, as well as 
#' arguments `ftol_abs` and `xtol_abs` for absolute changes 
#' (`0` by default), `maxeval` (maximum number of iterations, 1e4 by default), 
#' and `print_level` (default is `0`, no console output of  optimization, 
#' larger values are more verbose), are used by the chosen main optimization 
#' method. The number of initial EM iterations in `EM-DNM` can be set using 
#' argument `maxeval_em_dnm` (default is 100), and algorithm for direct
#' numerical optimization can be defined using argument `algorithm` 
#' (see [nloptr::nloptr()] for possible options).
#' 
#' For controlling these stopping criteria for the multistart phase, argument 
#' `control_restart` takes a list such as `list(ftol_rel = 0.01, print_level = 1)`. 
#' Default are as in the case of main optimization (which is always run once a
#' fter the restarts, using best solution from restarts as initial value),
#' except that `ftol_rel = 1e-6` and that `xtol_rel = 1e-4`.
#' Additionally, same options can be defined separately for the M-step of EM 
#' algorithm via list `control_mstep`. For `control_mstep`, the 
#' default value are `ftol_rel = 1e-6`, `xtol_rel = 1e-4`, and `maxeval = 100`.
#' 
#' @references Steven G. Johnson, The NLopt nonlinear-optimization package, http://github.com/stevengj/nlopt
#' @param observations Either the name of the response variable in `data`, or 
#' an `stslist` object (see [TraMineR::seqdef()]) containing the 
#' sequences. In case of multichannel data, `observations` should be a vector 
#' of response variable names in `data`, or a list of `stslist` objects.
#' @param n_states An integer > 1 defining the number of hidden states.
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
#' penalization is `0.5 * lambda * sum(eta^2)`. Note that with 
#' `method = "L-BFGS"` both objective function (log-likelihood) and 
#' the penalization term is scaled with number of non-missing observations. 
#' Default is `0`, but small values such as `1e-4` can help to ensure numerical 
#' stability of L-BFGS by avoiding extreme probabilities. See also argument 
#' `bound` for hard constraints.
#' @param method Optimization method used. Option `"EM"` uses EM
#' algorithm with L-BFGS in the M-step. Option `"DNM"` uses 
#' direct maximization of the log-likelihood, by default using L-BFGS. Option 
#' `"EM-DNM"` (the default) runs first a maximum of 10 iterations of EM and 
#' then switches to L-BFGS (but other algorithms of NLopt can be used).
#' @param bound Positive value defining the hard bounds for the working 
#' parameters \eqn{\eta}, which are used to avoid extreme probabilities and 
#' corresponding numerical issues especially in the M-step of EM algorithm. 
#' Default is 50, i.e., \eqn{-50<\eta<50}. Note that he bounds are not enforced 
#' for M-step in intercept-only case with `lambda = 0`.
#' @param store_data If `TRUE` (default), original data frame passed as `data` 
#' is stored to the model object. For large datasets, this can be set to 
#' `FALSE`, in which case you might need to pass the data separately to some 
#' post-prosessing functions.
#' @param ... Additional arguments to [nloptr::nloptr()] and EM algorithm. 
#' See details.
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
    data = NULL, time = NULL, id = NULL, state_names = NULL, 
    channel_names = NULL, inits = "random", init_sd = 2, restarts = 0L, 
    lambda = 0, method = "EM-DNM", bound = 50, control_restart = list(), 
    control_mstep = list(), store_data = TRUE, ...) {
  
  call <- match.call()
  model <- build_nhmm(
    observations, n_states, initial_formula, 
    transition_formula, emission_formula, data, time, id, state_names, 
    channel_names
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
  out <- fit_nhmm(model, inits, init_sd, restarts, lambda, method, bound, 
                  control, control_restart, control_mstep)
  end_time <- proc.time()
  out$estimation_results$time <- end_time - start_time
  attr(out, "call") <- call
  out
}
