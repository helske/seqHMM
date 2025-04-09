# Creates the Model Components for Non-homogeneous Hidden Markov Model
#' @noRd
create_base_nhmm <- function(data, id_var, time_var, n_states, state_names, 
                             emission_formula, initial_formula, transition_formula, 
                             cluster_formula = NA, cluster_names = "", 
                             scale = TRUE, fanhmm = FALSE) {
  # avoid CRAN check warnings due to NSE
  tmax <- y <- NULL
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
        responses, \(y) as.formula(
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
    C == length(unique(responses)), 
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
    state_names <- paste("State", seq_len(n_states))
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
        state_names <- replicate(D, state_names, simplify = FALSE)
      } else {
        lapply(seq_len(n_states), \(i) {
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
    names(state_names) <- cluster_names
  }
  data <- .check_data(data, id_var, time_var, responses)
  ids <- unique(data[[id_var]])
  times <- unique(data[[time_var]])
  n_sequences <- length(ids)
  symbol_names <- setNames(lapply(responses, \(y) levels(data[[y]])), responses)
  M <- lengths(symbol_names)
  n_obs <- sum(!is.na(data[, y, env = list(y = I(responses))])) / C
  if (n_obs == 0) {
    warning_("Responses contain only missing values.")
  }
  data[, tmax := max(time_var), by = id_var, 
       env = list(id_var = id_var, time_var = time_var)]
  data <- fill_time(data, id_var, time_var)
  sequence_lengths <- data[, 
                           sum(time_var <= tmax, na.rm = TRUE), 
                           by = id_var, 
                           env = list(
                             id_var = id_var, 
                             time_var = time_var
                           )]$V1
  data[, tmax := NULL]
  length_of_sequences <- length(unique(data[[time_var]]))
  if (fanhmm) {
    for (y in responses) {
      lag_obs <- paste0("lag_", y)
      y1 <- data[[y]][1] #the value is not used anywhere
      data[, lag_obs := shift(y, type = "lag", fill = y1), by = id, 
           env = list(id = id_var, y = y, lag_obs = lag_obs, y1 = y1)]
    }
  }
  X_pi <- model_matrix_initial_formula(
    initial_formula, data, n_sequences, n_states, id_var, scale = scale
  )
  np_pi <- (n_states - 1L) * nrow(X_pi)
  X_A <- model_matrix_transition_formula(
    transition_formula, data, n_sequences, length_of_sequences, n_states, 
    id_var, time_var, sequence_lengths, scale = scale
  )
  np_A <- n_states * (n_states - 1L) * nrow(X_A)
  
  X_B <- stats::setNames(
    lapply(responses, \(y) {
      model_matrix_emission_formula(
        emission_formula[[y]], data, n_sequences, length_of_sequences, n_states, 
        M[y], id_var, time_var, sequence_lengths, scale = scale
      )
    }),
    responses
  )
  np_B <- sum(n_states * (M - 1L) * vapply(X_B, \(x) nrow(x), 1L))
  if (mixture) {
    X_omega <- model_matrix_cluster_formula(
      cluster_formula, data, n_sequences, D, id_var, scale = scale
    )
    np_omega <- (D - 1L) * nrow(X_omega)
  }
  
  list(
    model = list(
      responses = responses, 
      time_variable = time_var,
      id_variable = id_var,
      X_pi = X_pi, 
      X_A = X_A, 
      X_B = X_B,
      X_omega = if(mixture) X_omega else NULL,
      initial_formula = initial_formula, 
      transition_formula = transition_formula,
      emission_formula = emission_formula,
      cluster_formula = if(mixture) cluster_formula else NULL,
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
      data = data
    ),
    extras = list(
      n_obs = n_obs,
      np_pi = np_pi,
      np_A = np_A,
      np_B = np_B,
      np_omega = if(mixture) np_omega else NULL
    )
  )
}
