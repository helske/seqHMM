# Creates the Model Components for Non-homogeneous Hidden Markov Model
#' @noRd
create_base_nhmm <- function(observations, data, time, id, n_states, 
                             state_names, channel_names,
                             initial_formula, transition_formula, 
                             emission_formula, cluster_formula = NA, 
                             cluster_names = "") {
  
  stopifnot_(
    !missing(n_states) && checkmate::test_int(x = n_states, lower = 2L), 
    "Argument {.arg n_states} must be a single integer larger than one."
  )
  stopifnot_(
    inherits(initial_formula, "formula"), 
    "Argument {.arg initial_formula} must be a {.cls formula} object.")
  stopifnot_(
    inherits(transition_formula, "formula"), 
    "Argument {.arg transition_formula} must be a {.cls formula} object.")
  stopifnot_(
    inherits(emission_formula, "formula"), 
    "Argument {.arg emission_formula} must be a {.cls formula} object.")
  n_clusters <- length(cluster_names)
  mixture <- n_clusters > 1L
  stopifnot_(
    !mixture || inherits(cluster_formula, "formula"), 
    "Argument {.arg cluster_formula} must be a {.cls formula} object.")
  
  n_states <- as.integer(n_states)
  if (is.null(state_names)) {
    state_names <- paste("State", seq_len(n_states))
    if (mixture) {
      state_names <- replicate(n_clusters, state_names, simplify = FALSE)
      names(state_names) <- cluster_names
    }
  } else {
    if (mixture) {
      names_is_vec <- !is.list(state_names) && length(state_names) == n_states
      stopifnot_(
        length(state_names) == n_clusters || names_is_vec,
        paste0(
          "For MNHMMs, {.arg state_names} should be a list of length ", 
          "{n_clusters}, the number of clusters, or a vector of length 
          {n_states}, number of hidden states."
        )
      )
      if (names_is_vec) {
        state_names <- replicate(n_clusters, state_names, simplify = FALSE)
      } else {
        lapply(seq_len(n_states), function(i) {
          stopifnot_(
            length(state_names[[i]]) == n_states,
            paste0(
              "Length of {.arg state_names[[{i}]]} is not equal to ",
              "{n_states}, the number of hidden states."
            )
          )
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
    }
  }
  
  icpt_only_i <- intercept_only(initial_formula)
  icpt_only_s <- intercept_only(transition_formula)
  icpt_only_o <- intercept_only(emission_formula)
  icpt_only_d <- ifelse(mixture, intercept_only(cluster_formula), TRUE)
  y_in_data <- checkmate::test_character(observations)
  if (!icpt_only_i || !icpt_only_s || !icpt_only_o || !icpt_only_d || y_in_data) {
    data <- .check_data(data, time, id)
    ids <- unique(data[[id]])
    times <- unique(data[[time]])
    n_sequences <- length(ids)
    length_of_sequences <- length(times)
    if (y_in_data) {
      channel_names <- observations
      observations <- lapply(
        observations,
        function(y) {
          stopifnot_(
            !is.null(data[[y]]), 
            "Can't find response variable {.var {y}} in {.arg data}."
          )
          x <- suppressMessages(
            seqdef(matrix(
              data[[y]], 
              n_sequences, 
              length_of_sequences, byrow = TRUE),
              id = ids
            )
          )
          colnames(x) <- sort(times)
          x
        }
      )
      observations <- .check_observations(observations, channel_names)
    }
  }
  if (!y_in_data) {
    observations <- .check_observations(observations, channel_names, 
                                        nhmm = TRUE)
    n_sequences <- attr(observations, "n_sequences")
    length_of_sequences <- attr(observations, "length_of_sequences")  
  }
  n_channels <- attr(observations, "n_channels")
  n_symbols <- attr(observations, "n_symbols")
  sequence_lengths <- attr(observations, "sequence_lengths")
  pi <- model_matrix_initial_formula(
    initial_formula, data, n_sequences, length_of_sequences, n_states, 
    time, id
  )
  A <- model_matrix_transition_formula(
    transition_formula, data, n_sequences, length_of_sequences, n_states, 
    time, id, sequence_lengths
  )
  B <- model_matrix_emission_formula(
    emission_formula, data, n_sequences, length_of_sequences, n_states, 
    n_symbols, n_channels, time, id, sequence_lengths
  )
  if (mixture) {
    omega <- model_matrix_cluster_formula(
      cluster_formula, data, n_sequences, n_clusters, time, id
    )
  } else {
    omega <- list(n_pars = 0, iv = FALSE, X_mean = NULL, X_sd = NULL)
  }
  list(
    model = list(
      observations = observations, 
      time_variable = if (is.null(time)) "time" else time,
      id_variable = if (is.null(id)) "id" else id,
      X_pi = pi$X, 
      X_A = A$X, 
      X_B = B$X,
      X_omega = if(mixture) omega$X else NULL,
      initial_formula = pi$formula, 
      transition_formula = A$formula,
      emission_formula = B$formula,
      cluster_formula = if(mixture) omega$formula else NULL,
      state_names = state_names,
      symbol_names = attr(observations, "symbol_names"),
      channel_names = attr(observations, "channel_names"),
      cluster_names = cluster_names,
      length_of_sequences = length_of_sequences,
      sequence_lengths = sequence_lengths,
      n_sequences = n_sequences,
      n_states = n_states,
      n_symbols = attr(observations, "n_symbols"),
      n_channels = attr(observations, "n_channels"),
      n_clusters = n_clusters
    ),
    extras = list(
      np_pi = pi$n_pars,
      np_A = A$n_pars,
      np_B = B$n_pars,
      np_omega = omega$n_pars,
      multichannel = ifelse(n_channels > 1, "multichannel_", ""),
      intercept_only = icpt_only_i && icpt_only_s && icpt_only_o && icpt_only_d
    )
  )
}
