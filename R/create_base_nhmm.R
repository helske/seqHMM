# Creates the Model Components for Non-homogeneous Hidden Markov Model
#' @noRd
create_base_nhmm <- function(data, id_var, time_var, n_states, state_names, 
                             emission_formula, initial_formula, transition_formula, 
                             cluster_formula = NA, cluster_names = "", 
                             scale = TRUE, prior_obs = "fixed", check = NULL, 
                             drop_levels = TRUE) {
  
  # avoid CRAN check warnings due to NSE
  .Ti <- y <- NULL
  stopifnot_(!missing(data), "Argument {.arg data} is missing.")
  stopifnot_(
    !missing(n_states) && checkmate::test_int(x = n_states, lower = 2L), 
    "Argument {.arg n_states} must be a single integer larger than 1."
  )
  stopifnot_(
    !missing(emission_formula),
    "Argument {.arg emission_formula} is missing."
  )
  if (inherits(emission_formula, "formula")) {
    responses <- get_responses(emission_formula)
    C <- length(responses)
    if (C > 1L) {
      rhs <- deparse1(emission_formula[[3L]])
      emission_formula <- lapply(
        responses, \(y) stats::as.formula(
          paste(y, " ~ ", rhs), 
          env = environment(emission_formula)
        )
      )
    } else {
      emission_formula <- list(emission_formula)
    }
  } else {
    responses <- vapply(emission_formula, get_responses, allow_mv = FALSE, "")
    C <- length(responses)
  }
  stopifnot_(
    C == n_unique(responses), 
    "Response names in {.arg responses} should be unique."
  )
  names(emission_formula) <- responses
  D <- length(cluster_names)
  mixture <- D > 1L
  stopifnot_(
    inherits(initial_formula, "formula"), 
    "Argument {.arg initial_formula} must be a {.cls formula} object.")
  stopifnot_(
    inherits(transition_formula, "formula"), 
    "Argument {.arg transition_formula} must be a {.cls formula} object.")
  stopifnot_(
    !mixture || inherits(cluster_formula, "formula"), 
    "Argument {.arg cluster_formula} must be a {.cls formula} object.")
  
  n_states <- as.integer(n_states)
  if (is.null(state_names)) {
    state_names <- as_factor(paste("State", seq_len(n_states)))
    if (mixture) {
      state_names <- replicate(D, state_names, simplify = FALSE)
      names(state_names) <- cluster_names
    }
  } else {
    if (mixture) {
      names_is_vec <- !is.list(state_names) && length(state_names) == n_states
      stopifnot_(
        length(state_names) == D || names_is_vec,
        paste0(
          "For MNHMMs, {.arg state_names} should be a list of length ", 
          "{D}, the number of clusters, or a vector of length 
          {n_states}, number of hidden states."
        )
      )
      if (names_is_vec) {
        state_names <- replicate(D, as_factor(state_names), simplify = FALSE)
      } else {
        state_names <- lapply(seq_len(n_states), \(i) {
          stopifnot_(
            length(state_names[[i]]) == n_states,
            paste0(
              "Length of {.arg state_names[[{i}]]} is not equal to ",
              "{n_states}, the number of hidden states."
            )
          )
          as_factor(state_names[[i]])
        })
      }
    } else {
      stopifnot_(
        length(state_names) == n_states,
        paste0(
          "Length of {.arg state_names} is not equal to {n_states}, the number", 
          " of hidden states."
        )
      )
      state_names <- as_factor(state_names)
    }
    names(state_names) <- cluster_names
  }
  data <- .check_data(data, id_var, time_var, responses)
  data <- fill_time(data, id_var, time_var)
  sequence_lengths <- data[, .Ti[1], by = id_var, 
                           env = list(id_var = id_var),
                           showProgress = FALSE]$V1
  data[, .Ti := NULL]
  symbol_names <- stats::setNames(
    lapply(responses, \(y) as_factor(levels(data[[y]]))), responses
  )
  M <- lengths(symbol_names)
  
  initial_formula <- check_formula(
    initial_formula, responses, "pi", data, id_var
  )
  transition_formula <- check_formula(
    transition_formula, responses, "A", data, id_var
  )
  feedback <- character(0)
  if (!is.null(resp <- attr(transition_formula, "responses"))) {
    feedback <- resp
  }
  autoregression <- stats::setNames(logical(C), responses)
  for (y in responses) {
    emission_formula[[y]] <- check_formula(
      emission_formula[[y]], responses, "B", data, id_var
    )
    if (!is.null(attr(emission_formula[[y]], "responses"))) {
      autoregression[y] <- TRUE
    }
  }
  autoregression <- responses[autoregression]
  fanhmm <- length(feedback) > 0L || length(autoregression) > 0L
  if (fanhmm) {
    W_X_B <- list()
    if (length(autoregression) > 0L && identical(prior_obs, "fixed")) {
      .idx <- setdiff(
        seq_row(data), 
        cumsum(c(1, head(sequence_lengths, -1)))
      )
      data <- data[.idx]
      sequence_lengths <- sequence_lengths - 1L
    } 
  }
  length_of_sequences <- n_unique(data[[time_var]])
  n_sequences <- n_unique(data[[id_var]])
  n_obs <- sum(!is.na(data[, y, env = list(y = I(responses))])) / C
  if (n_obs == 0) {
    warning_("Responses contain only missing values.")
  }
  if (drop_levels) {
    setdroplevels(data)
  }
  if (is.null(check)) {
    check <- n_sequences <= 1000
    if (!check) {
      warning_(
        c("Number of sequences is more than 1000, 
          disabling the rank check of design matrices.",
          i = "Explicitly use {.arg check_rank = FALSE} to 
          avoid this warning.")
      )
    }
  }
  X_pi <- model_matrix_initial_formula(
    initial_formula, data, n_sequences, n_states, id_var, scale = scale, 
    check = check
  )
  np_pi <- (n_states - 1L) * nrow(X_pi)
  X_A <- model_matrix_transition_formula(
    transition_formula, data, n_sequences, length_of_sequences, n_states, 
    id_var, time_var, sequence_lengths, scale = scale, check = check
  )
  np_A <- n_states * (n_states - 1L) * nrow(X_A[[1]])
  attr(transition_formula, "responses") <- NULL
  X_B <- stats::setNames(vector("list", C), responses)
  for (y in responses) {
    X_B[[y]] <- model_matrix_emission_formula(
      emission_formula[[y]], data, n_sequences, length_of_sequences, n_states, 
      M[y], id_var, time_var, sequence_lengths, scale = scale, check = check
    )
    attr(emission_formula[[y]], "responses") <- NULL
  }
  np_B <- sum(n_states * (M - 1L) * vapply(X_B, \(x) nrow(x[[1]]), 1L))
  if (mixture) {
    X_omega <- model_matrix_cluster_formula(
      cluster_formula, data, n_sequences, D, id_var, scale = scale, 
      check = check
    )
    np_omega <- (D - 1L) * nrow(X_omega)
  }
  if (length(autoregression) > 0L && !identical(prior_obs, "fixed")) {
    stopifnot_(
      is.list(prior_obs) && length(prior_obs) == C,
      c(x = "Argument {.arg prior_obs} must be {.val fixed}, or a list of 
          length {C}, the number of responses.",
        i = "Each element of the list must be a vector defining the prior 
          distribution of the response at time 'zero'."
      )
    )
    for (i in seq_len(C)) {
      stopifnot_(
        is_pmf(prior_obs[[i]], M[i]),
        x = "Argument {.arg prior_obs[[{i}]]} must be a valid probability 
        vector of length {M[i]}."
      )
    }
    prior_obs <- c(joint_probability(prior_obs))
    W_X_B <- create_W_X_B(
      data, id_var, time_var, symbol_names, n_sequences, emission_formula, 
      n_states, X_B
    )
  } else {
    prior_obs <- 0L
  }
  
  structure(
    list(
      responses = responses, 
      time_variable = time_var,
      id_variable = id_var,
      X_pi = X_pi, 
      X_A = X_A, 
      X_B = X_B,
      X_omega = if (mixture) X_omega else NULL,
      initial_formula = initial_formula, 
      transition_formula = transition_formula,
      emission_formula = emission_formula,
      cluster_formula = if (mixture) cluster_formula else NULL,
      state_names = state_names,
      symbol_names = symbol_names,
      cluster_names = cluster_names,
      length_of_sequences = length_of_sequences,
      sequence_lengths = sequence_lengths,
      n_sequences = n_sequences,
      n_states = n_states,
      n_symbols = M,
      n_channels = C,
      n_clusters = D,
      K_pi = nrow(X_pi),
      K_A = nrow(X_A[[1]]),
      K_B = vapply(X_B, \(x) nrow(x[[1]]), 1L),
      K_omega = if (mixture) nrow(X_omega) else NULL,
      data = data,
      feedback = feedback,
      autoregression = autoregression,
      W_X_B = if (fanhmm) W_X_B,
      prior_obs = if (fanhmm) prior_obs
    ),
    class = c(if (fanhmm) "fanhmm", ifelse(D > 1, "mnhmm", "nhmm")),
    nobs = n_obs,
    df = D * (np_pi + np_A + np_B) + if (mixture) np_omega else 0,
    np_pi = D * np_pi,
    np_A = D * np_A,
    np_B = D * np_B,
    np_omega = if (mixture) np_omega else 0
  )
}
