#' Average Marginal Predictions for Non-homogenous Hidden Markov Models
#' 
#' @param model A Hidden Markov Model of class `nhmm` or `mnhmm`.
#' @param variable Name of the variable of interest.
#' @param values Vector containing one or two values for `variable`.
#' @param marginalize_B_over Character string defining the dimensions over which
#' emission probabilities are marginalized. Default is `"sequences"`.
#' @param newdata Optional data frame which is used for marginalization.
#' @param nsim Non-negative integer defining the number of samples from the 
#' normal approximation of the model parameters used in 
#' computing the approximate quantiles of the estimates. If `0`, only point 
#' estimates are returned.
#' @param probs Vector defining the quantiles of interest. Default is 
#' `c(0.025, 0.5, 0.975)`.
#' @export
average_marginal_prediction <- function(
    model, variable, values, marginalize_B_over = "sequences", newdata = NULL, 
    nsim = 0, probs = c(0.025, 0.5, 0.975)) {
  # avoid warnings of NSEs
  cluster <- state <- estimate <- state_from <- state_to <- time_var <- 
    channel <- observation <- NULL
  stopifnot_(
    inherits(model, "nhmm") || inherits(model, "mnhmm"),
    "Argument {.arg model} must be a {.cls nhmm} or {.cls mnhmm} object."
  )
  marginalize_B_over <- match.arg(
    marginalize_B_over, c("sequences", "states", "clusters"), several.ok = TRUE)
  stopifnot_(
    marginalize_B_over != "clusters" || model$n_clusters > 1,
    "Cannot marginalize over clusters as {.arg model} is not a {.cls mnhmm} object."
  )
  stopifnot_(
    checkmate::test_count(nsim),
    "Argument {.arg nsim} should be a single non-negative integer."
  )
  stopifnot_(
    checkmate::test_string(x = variable), 
    "Argument {.arg variable} must be a single character string."
  )
  stopifnot_(
    length(values) != 2, 
    "Argument {.arg values} should contain two values for 
    variable {.var variable}.")
  if (is.null(newdata)) {
    time <- model$time_variable
    id <- model$id_variable
    stopifnot_(
      is.data.frame(newdata), 
      "Argument {.arg newdata} must be a {.cls data.frame} object."
    )
    stopifnot_(
      !is.null(newdata[[id]]), 
      "Can't find grouping variable {.var {id}} in {.arg newdata}."
    )
    stopifnot_(
      !is.null(newdata[[time]]), 
      "Can't find time index variable {.var {time}} in {.arg newdata}."
    )
    stopifnot_(
      !is.null(newdata[[variable]]), 
      "Can't find time variable {.var {variable}} in {.arg newdata}."
    )
  } else {
    stopifnot_(
      !is.null(model$data),
      "Model does not contain original data and argument {.arg newdata} is 
      {.var NULL}."
    )
    newdata <- model$data
  }
  # use same RNG seed so that the same samples of coefficients are drawn
  newdata[[variable]] <- values[1]
  seed <- sample(.Machine$integer.max, 1)
  set.seed(seed)
  pred <- predict(model, newdata, nsim, return_samples = TRUE)
  if (length(values) == 2) {
    newdata[[variable]] <- values[2]
    set.seed(seed)
    pred2 <- predict(model, newdata, nsim, return_samples = TRUE)
    pred <- mapply("-", pred, pred2, SIMPLIFY = FALSE)
  }
  T <- model$length_of_sequences
  N <- model$n_sequences
  S <- model$n_states
  M <- model$n_symbols
  C <- model$n_channels
  D <- model$n_clusters
  if (C == 1) {
    ids <- rownames(model$observations)
    times <- colnames(model$observations)
    symbol_names <- list(model$symbol_names)
  } else {
    ids <- rownames(model$observations[[1]])
    times <- colnames(model$observations[[1]])
    symbol_names <- model$symbol_names
  }
  marginalize <- c(
    switch(
      marginalize_B_over,
      "clusters" = c("cluster", "state", "id"),
      "states" = c("state", "id"),
      "sequences" = "id"), 
    "time", "channel", "observation")
  
  pi <- data.frame(
    cluster = rep(model$cluster_names, each = S * N),
    id = rep(ids, each = S),
    state = model$state_names,
    estimate = unlist(pred$pi)
  ) |> 
    dplyr::group_by(cluster, state) |>
    dplyr::summarise(estimate = mean(estimate))
  
  A <- data.frame(
    cluster = rep(model$cluster_names, each = S^2 * T * N),
    id = rep(ids, each = S^2 * T),
    time = rep(times, each = S^2),
    state_from = model$state_names,
    state_to = rep(model$state_names, each = S),
    estimate = unlist(pred$A)
  )  |> 
    dplyr::group_by(cluster, time, state_from, state_to) |>
    dplyr::summarise(estimate = mean(estimate)) |> 
    dplyr::rename(!!time_var := time)
  
  B <- data.frame(
    cluster =  rep(model$cluster_names, each = S * sum(M) * T * N),
    id = unlist(lapply(seq_len(C), function(i) rep(ids, each = S * M[i] * T))),
    time = unlist(lapply(seq_len(C), function(i) rep(times, each = S * M[i]))),
    state = model$state_names,
    channel = unlist(lapply(seq_len(C), function(i) {
      rep(model$channel_names[i], each = S * M[i]* T * N)
    })),
    observation = unlist(lapply(seq_len(C), function(i) {
      rep(symbol_names[[i]], each = S)
    })),
    estimate = unlist(pred$B)
  ) |> 
    dplyr::group_by(across(all_of(marginalize))) |>
    dplyr::summarise(estimate = mean(estimate)) |> 
    dplyr::rename(!!time_var := time)
  
  if (D > 1) {
    omega <- data.frame(
      cluster = model$cluster_names,
      id = rep(ids, each = D),
      estimate = c(pred$omega)
    )
    out <- list(omega = omega, pi = pi, A = A, B = B)  
  } else {
    out <- list(pi = pi, A = A, B = B)
  }
  class(out) <- "amp"
  attr(out, "seed") <- seed
  attr(out, "marginalize_B_over") <- marginalize_B_over
  out
}
