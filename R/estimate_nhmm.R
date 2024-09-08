#' Build and Estimate a Non-homogeneous Hidden Markov Model
#'
#' Function `estimate_nhmm` estimates a hidden Markov model object of class 
#' `nhmm` where initial, transition and emission probabilities 
#' (potentially) depend on covariates.
#'
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
#' are given using list components `beta_i_raw`, `beta_s_raw`, `beta_o_raw`, 
#' these are used as is, alternatively initial values can be given in terms of 
#' the initial state, transition, and emission probabilities using list 
#' components `initial_probs`, `emission_probs`, and `transition_probs`. These 
#' can also be mixed, i.e. you can give only `initial_probs` and `beta_s_raw`.
#' @param init_sd Standard deviation of the normal distribution used to generate
#' random initial values. Default is `2`. If you want to fix the initial values 
#' of the regression coefficients to zero, use `init_sd = 0`.
#' @param restarts Number of times to run optimization using random starting 
#' values. Default is 1.
#' @param threads Number of parallel threads for optimization with restarts. 
#' Default is 1.
#' @param store_data If `TRUE` (default), original data frame passed as `data` 
#' is stored to the model object. For large datasets, this can be set to 
#' `FALSE`, in which case you might need to pass the data separately to some 
#' post-prosessing functions.
#' @param verbose If `TRUE` (default), print progress messages during the final 
#' optimization (after restarts). Note that due to the design of `rstan`, the 
#' reason for the termination of optimization (e.g., reached the maximum number 
#' of iterations) is only available by the final message show with 
#' `verbose = TRUE`.
#' @param ... Additional arguments to [rstan::optimizing()]. Most importantly,
#' argument `iter` defines the maximum number of iterations for optimization.
#' The default is `2000`.
#' @return Object of class `nhmm`.
#' @export
#' @examples
#' data("mvad", package = "TraMineR")
#' 
#' mvad_alphabet <-
#'   c("employment", "FE", "HE", "joblessness", "school", "training")
#' mvad_labels <- c("employment", "further education", "higher education",
#'                  "joblessness", "school", "training")
#' mvad_scodes <- c("EM", "FE", "HE", "JL", "SC", "TR")
#' mvad_seq <- seqdef(mvad, 15:86, alphabet = mvad_alphabet,
#'                    states = mvad_scodes, labels = mvad_labels, xtstep = 6,
#'                    cpal = unname(colorpalette[[6]]))
#' 
#' set.seed(1)
#' \dontrun{
#' fit <- estimate_nhmm(mvad_seq, n_states = 3)
#' }
estimate_nhmm <- function(
    observations, n_states, initial_formula = ~1, 
    transition_formula = ~1, emission_formula = ~1, 
    data = NULL, time = NULL, id = NULL, state_names = NULL, channel_names = NULL, 
    inits = "random", init_sd = 2, restarts = 1L, threads = 1L, 
    store_data = TRUE, verbose = TRUE, ...) {
  
  model <- build_nhmm(
    observations, n_states, initial_formula, 
    transition_formula, emission_formula, data, time, id, state_names, 
    channel_names
    )
  stopifnot_(
    checkmate::test_flag(x = store_data), 
    "Argument {.arg store_data} must be a single {.cls logical} value.")
  if (store_data) {
    model$data <- data
  }
  fit_nhmm(model, inits, init_sd, restarts, threads, verbose, ...)
}
