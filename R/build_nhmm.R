#' Build a Non-homogeneous Hidden Markov Model
#'
#' @noRd
build_nhmm <- function(
    observations, n_states, initial_formula, 
    transition_formula, emission_formula, 
    data, data0, state_names, channel_names) {
  
  stopifnot_(
    inherits(initial_formula, "formula"), 
    "Argument {.arg initial_formula} must be a {.cls formula} object.")
  stopifnot_(
    inherits(transition_formula, "formula"), 
    "Argument {.arg transition_formula} must be a {.cls formula} object.")
  stopifnot_(
    inherits(emission_formula, "formula"), 
    "Argument {.arg emission_formula} must be a {.cls formula} object.")
  
  observations <- check_observations(observations, channel_names)
  channel_names <- attr(observations, "channel_names")
  n_channels <- attr(observations, "n_channels")
  n_sequences <- attr(observations, "n_sequences")
  length_of_sequences <- attr(observations, "length_of_sequences")
  n_symbols <- attr(observations, "n_symbols")
  n_states <- as.integer(n_states)
  if (is.null(state_names)) {
    state_names <- paste("State", seq_len(n_states))
  } else {
    if (length(state_names) != n_states) {
      stop("Length of 'state_names' is not equal to the number of hidden states.")
    }
  }
  if (intercept_only(initial_formula)) {
    init_type <- "c"
    vars <- "pi"
    n_pars <- n_states - 1L
    X_i <- matrix(1, n_sequences, 1)
    coef_names_initial <- "(Intercept)"
  } else {
    stopifnot_(
      is.data.frame(data0), 
      "Argument {.arg data0} must be a {.cls data.frame} object."
    )
    # ensure there is an intercept for which we define ordering constraint
    X_i <- model.matrix.lm(
      initial_formula, data = data0, na.action = na.pass
    )
    coef_names_initial <- colnames(X_i)
    X_i[is.na(X_i)] <- 0
    init_type <- "v"
    vars <- c("alpha_i", "beta_i")
    n_pars <- (n_states - 1L) * ncol(X_i)
  }
  
  if (intercept_only(transition_formula)) {
    A_type <- "c"
    vars <- "A"
    n_pars <- n_pars + n_states * (n_states - 1L)
    X_s <- array(1, c(length_of_sequences, n_sequences, 1L))
    coef_names_transition <- "(Intercept)"
  } else {
    stopifnot_(
      is.data.frame(data),
      "Argument {.arg data} must be a {.cls data.frame} object."
    )
    X_s <- model.matrix.lm(
      transition_formula, data = data, na.action = na.pass
    )
    coef_names_transition <- colnames(X_s)
    X_s[is.na(X_s)] <- 0
    dim(X_s) <- c(length_of_sequences, n_sequences, ncol(X_s))
    A_type <- "v"
    vars <- c(vars, "beta_s")
    n_pars <- n_pars + n_states * (n_states - 1L) * dim(X_s)[3]
  }
  
  if (intercept_only(emission_formula)) {
    B_type <- "c"
    vars <- "B"
    n_pars <- n_pars + n_channels * n_states * (n_symbols - 1L)
    X_o <- array(1, c(length_of_sequences, n_sequences, 1L))
    coef_names_emission <- "(Intercept)"
  } else {
    stopifnot_(
      is.data.frame(data), 
      "Argument {.arg data} must be a {.cls data.frame} object."
    )
    X_o <- model.matrix.lm(
      emission_formula, data = data, na.action = na.pass
    )
    coef_names_emission <- colnames(X_i)
    X_o[is.na(X_o)] <- 0
    dim(X_o) <- c(length_of_sequences, n_sequences, ncol(X_o))
    B_type <- "v"
    vars <- c(vars, "beta_o")
    n_pars <- n_pars + n_channels * n_states * (n_symbols - 1L) * dim(X_o)[3]
  }
  multichannel <- ifelse(n_channels > 1, "multichannel", "")
  structure(
    list(
      observations = observations, 
      X_initial = X_i, X_transition = X_s, X_emission = X_o,
      initial_formula = initial_formula, 
      transition_formula = transition_formula,
      emission_formula = emission_formula,
      state_names = state_names,
      symbol_names = attr(observations, "symbol_names"),
      channel_names = channel_names,
      length_of_sequences = length_of_sequences,
      n_sequences = n_sequences,
      n_symbols = n_symbols,
      n_states = n_states,
      n_channels = n_channels,
      coef_names_initial = coef_names_initial,
      coef_names_transition = coef_names_transition,
      coef_names_emission = coef_names_emission,
    ),
    class = "nhmm",
    nobs = attr(observations, "nobs"),
    df = n_pars,
    type = paste0(multichannel, "nhmm_", init_type, A_type, B_type)
  )
}
