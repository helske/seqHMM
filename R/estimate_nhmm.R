#' Estimate a Non-homogeneous Hidden Markov Model
#'
#' Function `estimate_nhmm` estimates a non-homogeneous hidden Markov model 
#' (NHMM) where initial, transition, and emission probabilities can depend on 
#' covariates. Transition and emission probabilities can also depend on past 
#' responses, in which case the model is called feedback-augmented NHMM 
#' (FAN-HMM) (Helske, 2025).
#' 
#' In case of FAN-HMM with autoregressive dependency on the observational level, 
#' (i.e. response \eqn{y_t} depend on \eqn{y_{t-1}}), the emission 
#' probabilities at the first time point need special attention. By default,
#' the model is initialized with fixed values for the first time point 
#' (`prior_obs = "fixed"`), meaning that if the input data consists of 
#' time points \eqn{t=1, 2, \ldots}, then the model is defined from \eqn{t=2} 
#' onward and the data on \eqn{t=1} is used only for defining the emission 
#' probabilities at \eqn{t=2}. Note that in this case also the initial state 
#' probabilities correspond to \eqn{t=2}. 
#' 
#' Alternatively, you can define `prior_obs` as a list of vectors, where 
#' the number of vectors is equal to the number of responses, and each vector 
#' gives the prior distribution for the response at \eqn{t=0}. For example,
#' if you have response variables `y` and `x`, where `y` has 3 categories and 
#' `x` 2 categories, you can define 
#' `prior_obs = list(y = c(0.5, 0.3, 0.2), x = c(0.7, 0.3))`.
#' These distributions are then used to marginalize out \eqn{y_0} and 
#' \eqn{x_0} in the relevant emission probabilities.
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
#' change of the objective function is less than `1e-12`, the absolute change 
#' is less than `1e-8` or the relative or absolute change of the working 
#' parameters eta is less than `1e-6`. These can be changed
#' by passing arguments `ftol_rel`, `ftol_abs`, `xtol_rel`, and `xtol_abs` 
#' via `...`. These, as well as, `maxeval` (maximum number of iterations, 
#' 1e4 by default), and `print_level` (default is `0`, no console output, 
#' larger values are more verbose), are used by the chosen main optimization 
#' method. The number of initial EM iterations in `EM-DNM` can be set using 
#' argument `maxeval_em_dnm` (default is 100), and algorithm for direct
#' numerical optimization can be defined using argument `algorithm` 
#' (see [nloptr::nloptr()] for possible options).
#' 
#' For controlling these stopping criteria for the multistart phase, argument 
#' `control_restart` takes a list such as `list(ftol_rel = 0.01, print_level = 1)`. 
#' Default are as in the case of main optimization (which is always run once
#' after the restarts, using best solution from restarts as initial value)
#' Additionally, same options can be defined separately for the M-step of EM 
#' algorithm via list `control_mstep`. For `control_mstep`, the 
#' default values are `ftol_rel = 1e-10`, and `maxeval = 1000`, and otherwise 
#' identical to previous defaults above.
#' 
#' @references 
#' Helske, J (2025). Feedback-augmented Non-homogeneous Hidden Markov Models for 
#' Longitudinal Causal Inference. arXiv preprint. <doi:10.48550/arXiv.2503.16014>.
#' 
#' Johnson, SG. The NLopt nonlinear-optimization package, http://github.com/stevengj/nlopt.
#' 
#' @param n_states An integer > 1 defining the number of hidden states.
#' @param initial_formula of class [formula()] for the
#' initial state probabilities. Left-hand side of the formula should be empty.
#' @param transition_formula of class [formula()] for the state transition 
#' probabilities. Left-hand side of the formula should be empty.
#' @param emission_formula of class [formula()] for the
#' state emission probabilities, or a list of such formulas in case of multiple 
#' response variables. The left-hand side of formulas define the responses. 
#' For multiple responses having same formula, you can use a form
#' `c(y1, y2) ~ x`, where `y1` and `y2` are the response variables.
#' @param data A data frame containing the variables used in the model 
#' formulas.
#' @param time Name of the time index variable in `data`.
#' @param id Name of the id variable in `data` identifying different 
#' sequences.
#' @param prior_obs Either `"fixed"` or a list of vectors given the prior 
#' distributions for the responses at time "zero". See details.
#' @param state_names A vector of optional labels for the hidden states. If this
#' is `NULL` (the default), numbered states are used.
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
#' @param bound Positive value defining the hard lower and upper bounds for the 
#' working parameters \eqn{\eta}, which are used to avoid extreme probabilities and 
#' corresponding numerical issues especially in the M-step of EM algorithm. 
#' Default is `Inf´, i.e., no bounds. Note that he bounds are not enforced 
#' for M-step in intercept-only case with `lambda = 0`.
#' @param control_restart Controls for restart steps, see details.
#' @param control_mstep Controls for M-step of EM algorithm, see details.
#' @param check_rank If `TRUE`, the rank of the design matrices are 
#' checked for identifiability issues. Default is `NULL`, in which case checks 
#' are performed only if the number of sequences is 1000 or less, as the QR 
#' decomposition quickly becomes computationally demanding. If check is not 
#' performed, a warning is given, which can be circumvented by explicitly 
#' using `check_rank = FALSE`.
#' @param ... Additional arguments to [nloptr::nloptr()] and EM algorithm. 
#' See details.
#' @return Object of class `nhmm` or `fanhmm`.
#' @export
#' @examples
#' data("mvad", package = "TraMineR")
#' 
#' d <- reshape(mvad, direction = "long", varying = list(15:86), 
#'   v.names = "activity")
#' 
#' \dontrun{
#' set.seed(1)
#' fit <- estimate_nhmm(n_states = 3,
#'   data = d, time = "time", id = "id", 
#'   emission_formula = activity ~ gcse5eq, initial_formula = ~ 1, 
#'   transition_formula = ~ male + gcse5eq,
#'   method = "DNM", maxeval = 2 # very small number of iterations for CRAN
#'   )
#' }
estimate_nhmm <- function(
    n_states, emission_formula, initial_formula = ~1, 
    transition_formula = ~1, 
    data, time, id, lambda = 0, prior_obs = "fixed", state_names = NULL, 
    inits = "random", init_sd = 2, restarts = 0L, 
    method = "EM-DNM", bound = Inf, control_restart = list(), 
    control_mstep = list(), check_rank = NULL, ...) {
  
  call <- match.call()
  model <- build_nhmm(
    n_states, emission_formula, initial_formula, transition_formula, 
    data, id, time, state_names, scale = TRUE, prior_obs, 
    check = check_rank
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
