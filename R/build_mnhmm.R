#' Build a Mixture Non-homogeneous Hidden Markov Model
#' @noRd
build_mnhmm <- function(
    observations, n_states, n_clusters, initial_formula, 
    transition_formula, emission_formula, cluster_formula,
    data, data0, state_names, channel_names, cluster_names) {
  
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
    inherits(cluster_formula, "formula"), 
    "Argument {.arg cluster_formula} must be a {.cls formula} object.")
  
  observations <- .check_observations(observations, channel_names)
  channel_names <- attr(observations, "channel_names")
  n_channels <- attr(observations, "n_channels")
  n_sequences <- attr(observations, "n_sequences")
  length_of_sequences <- attr(observations, "length_of_sequences")
  n_symbols <- attr(observations, "n_symbols")
  check_positive_integer(n_states, "n_states")
  n_states <- as.integer(n_states)
  check_positive_integer(n_clusters, "n_clusters")
  n_clusters <- as.integer(n_clusters)
  if (is.null(state_names)) {
    state_names <- paste("State", seq_len(n_states))
  } else {
    stopifnot_(
      length(state_names) == n_states,
      "Length of {.arg state_names} is not equal to the number of hidden 
      states."
    )
  }
  if (is.null(cluster_names)) {
    cluster_names <- paste("Cluster", seq_len(n_clusters))
  } else {
    stopifnot_(
      length(cluster_names) == n_clusters,
      "Length of {.arg cluster_names} is not equal to the number of clusters."
    )
  }
  if (intercept_only(initial_formula)) {
    init_type <- "c"
    n_pars <- n_states - 1L
    X_i <- matrix(1, n_sequences, 1)
    coef_names_initial <- "(Intercept)"
  } else {
    stopifnot_(
      is.data.frame(data0), 
      "If {.arg initial_formula} is provided, {.arg data0} must be a 
      {.cls data.frame} object."
    )
    X_i <- model.matrix.lm(
      initial_formula, data = data0, na.action = na.pass
    )
    coef_names_initial <- colnames(X_i)
    X_i[is.na(X_i)] <- 0
    init_type <- "v"
    n_pars <- (n_states - 1L) * ncol(X_i)
  }
  
  if (intercept_only(transition_formula)) {
    A_type <- "c"
    n_pars <- n_pars + n_states * (n_states - 1L)
    X_s <- array(1, c(length_of_sequences, n_sequences, 1L))
    coef_names_transition <- "(Intercept)"
  } else {
    stopifnot_(
      is.data.frame(data),
      "If {.arg transition_formula} is provided, {.arg data} must be a 
      {.cls data.frame} object."
    )
    X_s <- model.matrix.lm(
      transition_formula, data = data, na.action = na.pass
    )
    coef_names_transition <- colnames(X_s)
    X_s[is.na(X_s)] <- 0
    dim(X_s) <- c(length_of_sequences, n_sequences, ncol(X_s))
    A_type <- "v"
    n_pars <- n_pars + n_states * (n_states - 1L) * dim(X_s)[3]
  }
  
  if (intercept_only(emission_formula)) {
    B_type <- "c"
    n_pars <- n_pars + n_channels * n_states * (n_symbols - 1L)
    X_o <- array(1, c(length_of_sequences, n_sequences, 1L))
    coef_names_emission <- "(Intercept)"
  } else {
    stopifnot_(
      is.data.frame(data), 
      "If {.arg emission_formula} is provided, {.arg data} must be a 
      {.cls data.frame} object."
    )
    X_o <- model.matrix.lm(
      emission_formula, data = data, na.action = na.pass
    )
    coef_names_emission <- colnames(X_i)
    X_o[is.na(X_o)] <- 0
    dim(X_o) <- c(length_of_sequences, n_sequences, ncol(X_o))
    B_type <- "v"
    n_pars <- n_pars + n_channels * n_states * (n_symbols - 1L) * dim(X_o)[3]
  }
  
  if (intercept_only(cluster_formula)) {
    theta_type <- "c"
    n_pars <- n_states - 1L
    X_d <- matrix(1, n_sequences, 1)
    coef_names_cluster <- "(Intercept)"
  } else {
    stopifnot_(
      is.data.frame(data0), 
      "If {.arg cluster_formula} is provided, {.arg data0} must be a 
      {.cls data.frame} object."
    )
    X_d <- model.matrix.lm(
      cluster_formula, data = data0, na.action = na.pass
    )
    coef_names_cluster <- colnames(X_i)
    X_d[is.na(X_d)] <- 0
    theta_type <- "v"
  }
  n_pars <- n_clusters * n_pars
  multichannel <- ifelse(n_channels > 1, "multichannel_", "")
  structure(
    list(
      observations = observations, 
      X_initial = X_i, X_transition = X_s, X_emission = X_o, X_cluster = X_d,
      initial_formula = initial_formula, 
      transition_formula = transition_formula,
      emission_formula = emission_formula,
      cluster_formula = cluster_formula,
      state_names = state_names,
      symbol_names = attr(observations, "symbol_names"),
      channel_names = channel_names,
      length_of_sequences = length_of_sequences,
      n_sequences = n_sequences,
      n_symbols = n_symbols,
      n_states = n_states,
      n_channels = n_channels,
      n_clusters = n_clusters,
      coef_names_initial = coef_names_initial,
      coef_names_transition = coef_names_transition,
      coef_names_emission = coef_names_emission,
      coef_names_cluster = coef_names_cluster
    ),
    class = "mnhmm",
    nobs = attr(observations, "nobs"),
    df = n_pars,
    type = paste0(multichannel, "mnhmm_", init_type, A_type, B_type, theta_type)
  )
}
