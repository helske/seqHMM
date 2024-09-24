#' Average Marginal Effects for Non-homogenous Hidden Markov Models
#' 
#' @param model A Hidden Markov Model of class `nhmm` or `mnhmm`.
#' @param variable Name of the variable of interest.
#' @param values Vector containing one or two values for `variable`.
#' @param marginalize_B_over Character string defining the dimensions over which
#' emission probabilities are marginalized. Default is `"sequences"`. Other 
#' options are `"states"` and `"clusters"`.
#' @param newdata Optional data frame which is used for marginalization.
#' @param ... Further arguments passed to specific methods.
#' @rdname ame
#' @export
ame <- function(model, variable, values, ...) {
  UseMethod("ame", model)
}
#' @rdname ame
#' @export
ame.nhmm <- function(
    model, variable, values, marginalize_B_over = "sequences", newdata = NULL, 
    ...) {
  # avoid warnings of NSEs
  cluster <- state <- estimate <- state_from <- state_to <- time_var <- 
    channel <- observation <- NULL
  stopifnot_(
    attr(model, "intercept_only") == FALSE,
    "Model does not contain any covariates."
  )
  marginalize_B_over <- match.arg(
    marginalize_B_over, c("sequences", "states", "clusters"))
  stopifnot_(
    marginalize_B_over != "clusters",
    "Cannot marginalize over clusters as {.arg model} is not a {.cls mnhmm} object."
  )
  stopifnot_(
    checkmate::test_string(x = variable), 
    "Argument {.arg variable} must be a single character string."
  )
  stopifnot_(
    length(values) == 2, 
    "Argument {.arg values} should contain two values for 
    variable {.var variable}.")
  time <- model$time_variable
  id <- model$id_variable
  if (!is.null(newdata)) {
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
  pred <- get_probs(model, newdata)
  newdata[[variable]] <- values[2]
  pred2 <- get_probs(model, newdata)
  pars <- c("initial_probs", "transition_probs", "emission_probs")
  for (i in pars) {
    pred[[i]]$estimate <- pred[[i]]$estimate - pred2[[i]]$estimate
  }
  channel <- if (model$n_channels > 1) "channel" else NULL
  group_by_B <- switch(
    marginalize_B_over,
    "states" = c("time", channel, "observation"),
    "sequences" = c("time", "state", channel, "observation")
  )
  
  qs <- function(x, probs) {
    x <- quantile(x, probs)
    names(x) <- paste0("q", 100 * probs)
    as.data.frame(t(x))
  }
  out <- list()
  out$initial_probs <- pred$initial_probs |> 
    dplyr::group_by(state) |>
    dplyr::summarise(estimate = mean(estimate)) |> 
    dplyr::ungroup()
  
  out$transition_probs <- pred$transition_probs |> 
    dplyr::group_by(time, state_from, state_to) |>
    dplyr::summarise(estimate = mean(estimate)) |> 
    dplyr::ungroup() |> 
    dplyr::rename(!!time := time)
  
  out$emission_probs <- pred$emission_probs |> 
    dplyr::group_by(dplyr::pick(group_by_B)) |>
    dplyr::summarise(estimate = mean(estimate)) |> 
    dplyr::ungroup() |> 
    dplyr::rename(!!time := time)
  
  class(out) <- "amp"
  attr(out, "model") <- "nhmm"
  attr(out, "marginalize_B_over") <- marginalize_B_over
  out
}

#' @rdname ame
#' @export
ame.mnhmm <- function(
    model, variable, values, marginalize_B_over = "sequences", newdata = NULL, 
    ...) {
  # avoid warnings of NSEs
  cluster <- state <- estimate <- state_from <- state_to <- time_var <- 
    channel <- observation <- NULL
  stopifnot_(
    attr(model, "intercept_only") == FALSE,
    "Model does not contain any covariates."
  )
  marginalize_B_over <- match.arg(
    marginalize_B_over, c("sequences", "states", "clusters"), several.ok = TRUE)
  stopifnot_(
    marginalize_B_over != "clusters" || model$n_clusters > 1,
    "Cannot marginalize over clusters as {.arg model} is not a {.cls mnhmm} object."
  )
  stopifnot_(
    checkmate::test_string(x = variable), 
    "Argument {.arg variable} must be a single character string."
  )
  stopifnot_(
    length(values) == 2, 
    "Argument {.arg values} should contain two values for 
    variable {.var variable}.")
  time <- model$time_variable
  id <- model$id_variable
  if (!is.null(newdata)) {
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
  pred <- get_probs(model, newdata)
  
  newdata[[variable]] <- values[2]
  pred2 <- get_probs(model, newdata)
  pars <- c("initial_probs", "transition_probs", "emission_probs",
            "cluster_probs")
  for (i in pars) {
    pred[[i]]$estimate <- pred[[i]]$estimate - pred2[[i]]$estimate
  }
  channel <- if (model$n_channels > 1) "channel" else NULL
  group_by_B <- switch(
    marginalize_B_over,
    "clusters" = c("time", channel, "observation"),
    "states" = c("cluster", "time", channel, "observation"),
    "sequences" = c("cluster", "time", "state", channel, "observation")
  )
  
  out <- list()
  out$initial_probs <- pred$initial_probs |> 
    dplyr::group_by(cluster, state) |>
    dplyr::summarise(estimate = mean(estimate)) |> 
    dplyr::ungroup()
  
  out$transition_probs <- pred$transition_probs |> 
    dplyr::group_by(cluster, time, state_from, state_to) |>
    dplyr::summarise(estimate = mean(estimate)) |> 
    dplyr::ungroup() |> 
    dplyr::rename(!!time := time)
  
  out$emission_probs <- pred$emission_probs |> 
    dplyr::group_by(dplyr::pick(group_by_B)) |>
    dplyr::summarise(estimate = mean(estimate)) |> 
    dplyr::ungroup() |> 
    dplyr::rename(!!time := time)
  
  out$cluster_probs <- pred$cluster_probs |> 
    dplyr::group_by(cluster) |>
    dplyr::summarise(estimate = mean(estimate)) |> 
    dplyr::ungroup()
  
  class(out) <- "amp"
  attr(out, "model") <- "mnhmm"
  attr(out, "marginalize_B_over") <- marginalize_B_over
  out
}
