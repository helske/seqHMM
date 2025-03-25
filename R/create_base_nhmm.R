# Creates the Model Components for Non-homogeneous Hidden Markov Model
#' @noRd
create_base_nhmm <- function(responses, data, id_var, time_var, n_states, state_names, 
                             initial_formula, transition_formula, 
                             emission_formula, cluster_formula = NA, 
                             cluster_names = "", scale = TRUE, 
                             fanhmm = FALSE) {
  # avoid CRAN check warnings due to NSE
  tmax <- y <- NULL
  stopifnot_(!missing(data), "Argument {.arg data} is missing.")
  stopifnot_(
    !missing(n_states) && checkmate::test_int(x = n_states, lower = 2L), 
    "Argument {.arg n_states} must be a single integer larger than 1."
  )
  stopifnot_(
    !missing(responses) && checkmate::test_character(x = responses), 
    "Argument {.arg responses} must be a character vector defining the response 
    variable(s) in the {.arg data}."
  )
  stopifnot_(
    length(responses) == length(unique(responses)), 
    "Response names in {.arg responses} should be unique."
  )
  n_clusters <- length(cluster_names)
  mixture <- n_clusters > 1L
  
  stopifnot_(
    inherits(initial_formula, "formula"), 
    "Argument {.arg initial_formula} must be a {.cls formula} object.")
  stopifnot_(
    inherits(transition_formula, "formula"), 
    "Argument {.arg transition_formula} must be a {.cls formula} object.")
  stopifnot_(
    inherits(emission_formula, "formula"), 
    "Argument {.arg emission_formula} must be a {.cls formula} object.")
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
  
  n_channels <- length(responses)
  data <- .check_data(data, id_var, time_var, responses)
  ids <- unique(data[[id_var]])
  times <- unique(data[[time_var]])
  n_sequences <- length(ids)
  symbol_names <- lapply(responses, function(y) levels(data[[y]]))
  n_symbols <- lengths(symbol_names)
  if (n_channels == 1) symbol_names <- symbol_names[[1]]
  n_obs <- sum(!is.na(data[, y, env = list(y = I(responses))])) / n_channels
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
    lag_obs <- paste0("lag_", responses)
    data[[lag_obs]] <- group_lag(data, id_var, responses)
    ar <- lag_obs %in% attr(stats::terms(emission_formula), "term.labels")
  } else {
    ar <- FALSE
  }
  if (ar) n_obs <- n_obs - n_sequences
  
  pi <- model_matrix_initial_formula(
    initial_formula, data, n_sequences, n_states, id_var, scale = scale
  )
  A <- model_matrix_transition_formula(
    transition_formula, data, n_sequences, length_of_sequences, n_states, 
    id_var, time_var, sequence_lengths, scale = scale
  )
  B <- model_matrix_emission_formula(
    emission_formula, data, n_sequences, length_of_sequences, n_states, 
    n_symbols, id_var, time_var, sequence_lengths, scale = scale, 
    autoregression = ar
  )
  if (mixture) {
    omega <- model_matrix_cluster_formula(
      cluster_formula, data, n_sequences, n_clusters, id_var, scale = scale
    )
  } else {
    omega <- list(n_pars = 0, iv = FALSE, X_mean = NULL)
  }
  
  list(
    model = list(
      responses = responses, 
      time_variable = time_var,
      id_variable = id_var,
      X_pi = pi$X, 
      X_A = A$X, 
      X_B = B$X,
      X_omega = if(mixture) omega$X else NULL,
      initial_formula = pi$formula, 
      transition_formula = A$formula,
      emission_formula = B$formula,
      cluster_formula = if(mixture) omega$formula else NULL,
      state_names = state_names,
      symbol_names = symbol_names,
      cluster_names = cluster_names,
      length_of_sequences = length_of_sequences,
      sequence_lengths = sequence_lengths,
      n_sequences = n_sequences,
      n_states = n_states,
      n_symbols = n_symbols,
      n_channels = n_channels,
      n_clusters = n_clusters,
      data = data
    ),
    extras = list(
      n_obs = n_obs,
      np_pi = pi$n_pars,
      np_A = A$n_pars,
      np_B = B$n_pars,
      np_omega = omega$n_pars,
      multichannel = ifelse(n_channels > 1, "multichannel_", ""),
      intercept_only = icpt_only_i && icpt_only_s && icpt_only_o && icpt_only_d
    )
  )
}
