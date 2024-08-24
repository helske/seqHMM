#' Build and Estimate a Non-homogeneous Hidden Markov Model
#'
#' Function `estimate_nhmm` estimates a hidden Markov model object of class 
#' `nhmm` where initial, transition and emission probabilities 
#' (potentially) depend on covariates.
#'
#' 
#' @param observations An `stslist` object 
#' (see [TraMineR::seqdef()]) containing the sequences.
#' @param n_states A positive integer defining the number of hidden states.
#' @param initial_formula of class [formula()] for the
#' initial state probabilities.
#' @param transition_formula of class [formula()] for the
#' state transition probabilities.
#' @param emission_formula of class [formula()] for the
#' state emission probabilities.
#' @param data A data frame containing the variables used in the transition and 
#' emission formulas. Data should be sorted so the first T rows corresponds to 
#' observations of the first sequence and so forth.
#' @param data0 A data frame containing the variables used in the initial 
#' state formula. Data should be sorted so that the first row corresponds to 
#' covariates of the first sequence and so forth.
#' @param state_names A vector of optional labels for the hidden states. If this
#' is `NULL` (the default), numbered states are used.
#' @param channel_names A vector of optional names for the channels. If this
#' is `NULL` (the default), numbered channels are used.
#' @param inits Optional initial values for the initial state, transition, 
#' emission, and mixture probabilities. Either a list with `initial_probs`, 
#' `emission_probs`, `transition_probs`, `cluster_probs`, or `"random"`.
#' @param init_sd Standard deviation of the normal distribution used to generate
#' random initial values. Default is `2`. If you want to fix the initial values 
#' of the regression coefficients to zero, use `init_sd = 0`.
#' @param restarts Number of times to run optimization using random starting 
#' values. Default is 1.
#' @param threads Number of parallel threads for optimization with restarts. 
#' Default is 1.
#' @param ... Additional arguments to [rstan::optimizing()].
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
#' mvad_seq <- seqdef(mvad, 17:86, alphabet = mvad_alphabet,
#'                    states = mvad_scodes, labels = mvad_labels, xtstep = 6,
#'                    cpal = unname(colorpalette[[6]]))
#' 
#' set.seed(1)
#' fit <- estimate_nhmm(mvad_seq, n_states = 3)
estimate_nhmm <- function(
    observations, n_states, initial_formula = ~1, 
    transition_formula = ~1, emission_formula = ~1, 
    data = NULL, data0 = NULL, state_names = NULL, channel_names = NULL, 
    inits = "random", init_sd = 2, restarts = 1L, threads = 1L, ...) {
  
  model <- build_nhmm(
    observations, n_states, initial_formula, 
    transition_formula, emission_formula, data, data0, state_names, 
    channel_names
    )
  fit_nhmm(model, inits, init_sd, restarts, threads, ...)
}
