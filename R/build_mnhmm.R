#' Build a Hidden Markov Model with Covariates
#'
#' Function \code{build_nhmm} constructs a non-homogenous hidden Markov model object of class \code{mnhmm}.
#'
#' The returned model contains some attributes such as \code{nobs} and \code{df},
#' which define the number of observations in the  model and the number of estimable
#' model parameters, used in computing BIC.
#' When computing \code{nobs} for a multichannel model with \eqn{C} channels,
#' each observed value in a single channel amounts to \eqn{1/C} observation,
#' i.e. a fully observed time point for a single sequence amounts to one observation.
#' For the degrees of freedom \code{df}, zero probabilities of the initial model are
#' defined as structural zeroes.
#' @export
#' @param observations An \code{stslist} object 
#' (see \code{\link[TraMineR]{seqdef}}) containing the sequences.
#' @param id_variable Name of the grouping variable in the data.
#' @param time_variable Name of the time index variable in the data.
#' @param initial_formula of class \code{\link{formula}} for the
#' initial state probabilities.
#' @param transition_formula of class \code{\link{formula}} for the
#' state transition probabilities.
#' @param emission_formula of class \code{\link{formula}} for the
#' state emission probabilities.
#' @param data A data frame containing the variables used in the transition and 
#' emission formulas.
#' @param init_data A data frame containing the variables used in the initial 
#' state formula.
#'
#' @return Object of class \code{mnhmm}.
#' @examples

build_mnhmm <- function(
    observations, n_states, n_clusters, initial_formula, 
    transition_formula, emission_formula, data = NULL, 
    data0 = NULL,
    state_names = NULL, channel_names = NULL) {
  
  observations <- check_observations(observations, channel_names)
  n_channels <- attr(observations, "n_channels")
  n_sequences <- attr(observations, "n_sequences")
  length_of_sequences <- attr(observations, "length_of_sequences")
  n_symbols <- attr(observations, "n_symbols")
  
  if (is.null(state_names)) {
    state_names <- paste("State", seq_len(n_states))
  } else {
    if (length(state_names) != n_states) {
      stop("Length of 'state_names' is not equal to the number of hidden states.")
    }
  }
  vars <- c("alpha_i", "beta_i", "beta_s", "beta_o", "rho", "A", "B")
  # this could be optimized so to avoid duplicate variables 
  # (i.e. use one common X in Stan with suitable subsetting of columns)
  if (is.null(init_data)) {
    X_i <- matrix(0, n_sequences, 0L)
  } else {
    X_i <- model.matrix.lm(initial_formula, data = init_data, na.action = na.pass)
    X_i[is.na(X_i)] <- 0
    if (ncol(X_i) < 1L) {
      stop("Initial state probability formula should contain at least one term.")
    }
    cnames <- colnames(X_i)
    if ("(Intercept)" %in% cnames) {
      X_i <- X_i[, cnames != "(Intercept)", drop = FALSE]
    }
  }
  K_i <- ncol(X_i)
  X_i <- t(X_i)
  if (K_i == 0L) vars <- vars[-2L]
  if (is.null(data)) {
    X_s <- X_o <- array(1, c(length_of_sequences, 1, n_sequences))
    K_s <- K_o <- 1L
  } else {
    X_s <- model.matrix.lm(transition_formula, data = data, na.action = na.pass)
    X_s[is.na(X_s)] <- 0
    if (ncol(X_s) < 1L) {
      stop("State transition probability formula should contain at least one term.")
    }
    K_s <- ncol(X_s)
    X_s <- t(X_s)
    dim(X_s) <- c(length_of_sequences, K_s, n_sequences)
    
    X_o <- model.matrix.lm(transition_formula, data = data, na.action = na.pass)
    X_o[is.na(X_o)] <- 0
    cnames <- colnames(X_o)
    if (ncol(X_o) < 1L) {
      stop("Emission probability formula should at least one term.")
    }
    K_o <- ncol(X_o)
    X_o <- t(X_o)
    dim(X_o) <- c(length_of_sequences, K_o, n_sequences)
  }
  structure(
    list(
      observations = observations, 
      X_initial = X_i, X_transition = X_s, X_emission = X_o,
      state_names = state_names,
      symbol_names = attr(observations, "symbol_names"),
      channel_names = channel_names,
      length_of_sequences = length_of_sequences,
      n_sequences = n_sequences,
      n_symbols = n_symbols, 
      n_states = n_states,
      n_channels = n_channels
    ),
    class = "mnhmm",
    nobs = attr(observations, "nobs"),
    df = n_states^2 - 1L + n_symbols * (n_states - 1),
    type = "mnhmm"
  )
}
