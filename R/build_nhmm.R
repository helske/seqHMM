#' Build a Non-homogeneous Hidden Markov Model
#'
#' @noRd
build_nhmm <- function(
    observations, n_states, initial_formula = NULL, 
    transition_formula = NULL, emission_formula = NULL, 
    data = NULL, data0 = NULL, state_names = NULL, channel_names = NULL) {
  
  observations <- check_observations(observations, channel_names)
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
  if (is.null(initial_formula)) {
    init_type <- "c"
    vars <- "pi"
    n_pars <- n_states - 1L
    X_i <- NULL
  } else {
    if (inherits(initial_formula, "formula") && is.data.frame(data0)) {
      # ensure there is an intercept for which we define ordering constraint
      X_i <- model.matrix.lm(
        update(initial_formula, ~ . + 1), 
        data = data0, na.action = na.pass)
      X_i[is.na(X_i)] <- 0
      cnames <- colnames(X_i)
      # remove intercept as that is handled separately
      if ("(Intercept)" %in% cnames) {
        X_i <- X_i[, cnames != "(Intercept)", drop = FALSE]
      }
      init_type <- "v"
      vars <- c("alpha_i", "beta_i")
      n_pars <- (n_states - 1L) * (ncol(X_i) + 1L)
    } else {
      if (!inherits(initial_formula, "formula")) {
        stop(
          paste0("Argument 'initial_formula' should be either 'formula' ",
                 "object or NULL.")
        )
      }
      if (!is.data.frame(data0)) {
        stop("Argument 'data0' should be a data.frame or NULL.")
      }
    }
  }
  if (is.null(transition_formula)) {
    A_type <- "c"
    vars <- "A"
    n_pars <- n_pars + n_states * (n_states - 1L)
    X_s <- NULL
  } else {
    if (inherits(transition_formula, "formula") && is.data.frame(data)) {
      X_s <- model.matrix.lm(
        transition_formula, data = data, na.action = na.pass
      )
      X_s[is.na(X_s)] <- 0
      dim(X_s) <- c(length_of_sequences, n_sequences, ncol(X_s))
      A_type <- "v"
      vars <- c(vars, "beta_s")
      n_pars <- n_pars + n_states * (n_states - 1L) * dim(X_s)[3]
    } else {
      if (!inherits(transition_formula, "formula")) {
        stop(
          paste0("Argument 'transition_formula' should be either 'formula' ",
                 "object or NULL.")
        )
      }
      if (!is.data.frame(data)) {
        stop("Argument 'data' should be a data.frame or NULL.")
      }
    }
  }
  if (is.null(emission_formula)) {
    B_type <- "c"
    vars <- "B"
    n_pars <- n_pars + n_channels * n_states * (n_symbols - 1L)
    X_o <- NULL
  } else {
    if (inherits(emission_formula, "formula") && is.data.frame(data)) {
      X_o <- model.matrix.lm(
        emission_formula, data = data, na.action = na.pass
      )
      X_o[is.na(X_o)] <- 0
      dim(X_o) <- c(length_of_sequences, n_sequences, ncol(X_p))
      B_type <- "v"
      vars <- c(vars, "beta_o")
      n_pars <- n_pars + n_channels * n_states * (n_symbols - 1L) * dim(X_o)[3]
    } else {
      if (!inherits(emission_formula, "formula")) {
        stop(
          paste0("Argument 'emission_formula' should be either 'formula' ",
                 "object or NULL.")
        )
      }
      if (!is.data.frame(data)) {
        stop("Argument 'data' should be a data.frame or NULL.")
      }
    }
  }
  multichannel <- ifelse(n_channels > 1, "multichannel", "")
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
    class = "nhmm",
    nobs = attr(observations, "nobs"),
    df = n_pars,
    type = paste0(multichannel, "nhmm_", init_type, A_type, B_type)
  )
}
