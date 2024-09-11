#' Average Marginal Effects for Non-homogenous Hidden Markov Models
#' 
#' @param model A Hidden Markov Model of class `nhmm` or `mnhmm`.
#' @param variable Name of the variable of interest.
#' @param values Vector containing one or two values for `variable`.
#' @param marginalize_B_over Character string defining the dimensions over which
#' emission probabilities are marginalized. Default is `"sequences"`. Other 
#' options are `"states"` and `"clusters"`.
#' @param newdata Optional data frame which is used for marginalization.
#' @param nsim Non-negative integer defining the number of samples from the 
#' normal approximation of the model parameters used in 
#' computing the approximate quantiles of the estimates. If `0`, only point 
#' estimates are returned.
#' @param probs Vector defining the quantiles of interest. Default is 
#' `c(0.025, 0.5, 0.975)`.
#' @rdname ame
#' @export
ame <- function(model, variable, values, ...) {
  UseMethod("ame", model)
}
#' @rdname ame
#' @export
ame.nhmm <- function(
    model, variable, values, marginalize_B_over = "sequences", newdata = NULL, 
    nsim = 0, probs = c(0.025, 0.5, 0.975)) {
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
    checkmate::test_count(nsim),
    "Argument {.arg nsim} should be a single non-negative integer."
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
  seed <- sample(.Machine$integer.max, 1)
  set.seed(seed)
  pred <- predict(model, newdata, nsim, return_samples = TRUE, 
                  dontchange_colnames = TRUE)
  newdata[[variable]] <- values[2]
  set.seed(seed)
  pred2 <- predict(model, newdata, nsim, return_samples = TRUE, 
                   dontchange_colnames = TRUE)
  pars <- c("initial_probs", "transition_probs", "emission_probs")
  for (i in pars) {
    pred[[i]]$estimate <- pred[[i]]$estimate - pred2[[i]]$estimate
  }
  if (nsim > 0) {
    for (i in pars) {
      pred$samples[[i]]$estimate <- pred$samples[[i]]$estimate - 
        pred2$samples[[i]]$estimate
    }
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
  out$initial_probs <- cbind(
    pred$initial_probs |> 
      dplyr::group_by(state) |>
      dplyr::summarise(estimate = mean(estimate)) |> 
      dplyr::ungroup(),
    pred$samples$initial_probs |> 
      dplyr::group_by(state, replication) |>
      dplyr::summarise(estimate = mean(estimate)) |>
      dplyr::group_by(state) |>
      dplyr::summarise(qs(estimate, probs)) |>
      dplyr::ungroup() |> 
      dplyr::select(-state)
  )
  
  out$transition_probs <- cbind(
    pred$transition_probs |> 
      dplyr::group_by(time, state_from, state_to) |>
      dplyr::summarise(estimate = mean(estimate)) |> 
      dplyr::ungroup(),
    pred$samples$transition_probs |> 
      dplyr::group_by(time, state_from, state_to, replication) |>
      dplyr::summarise(estimate = mean(estimate)) |>
      dplyr::group_by(time, state_from, state_to) |>
      dplyr::summarise(qs(estimate, probs)) |>
      dplyr::ungroup() |> 
      dplyr::select(-c(time, state_from, state_to))
  ) |>  dplyr::rename(!!time := time)
  
  out$emission_probs <- cbind(
    pred$emission_probs |> 
      dplyr::group_by(dplyr::pick(group_by_B)) |>
      dplyr::summarise(estimate = mean(estimate)) |> 
      dplyr::ungroup(),
    pred$samples$emission_probs |> 
      dplyr::group_by(dplyr::pick(c(group_by_B, "replication"))) |>
      dplyr::summarise(estimate = mean(estimate)) |>
      dplyr::group_by(dplyr::pick(group_by_B)) |>
      dplyr::summarise(qs(estimate, probs)) |>
      dplyr::ungroup() |> 
      dplyr::select(dplyr::starts_with("q"))
  ) |> dplyr::rename(!!time := time)
  
  class(out) <- "amp"
  attr(out, "model") <- "nhmm"
  attr(out, "seed") <- seed
  attr(out, "marginalize_B_over") <- marginalize_B_over
  out
}

#' @rdname ame
#' @export
ame.mnhmm <- function(
    model, variable, values, marginalize_B_over = "sequences", newdata = NULL, 
    nsim = 0, probs = c(0.025, 0.5, 0.975)) {
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
    checkmate::test_count(nsim),
    "Argument {.arg nsim} should be a single non-negative integer."
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
  seed <- sample(.Machine$integer.max, 1)
  set.seed(seed)
  pred <- predict(model, newdata, nsim, return_samples = TRUE, 
                  dontchange_colnames = TRUE)
 
  newdata[[variable]] <- values[2]
  set.seed(seed)
  pred2 <- predict(model, newdata, nsim, return_samples = TRUE, 
                   dontchange_colnames = TRUE)
  pars <- c("initial_probs", "transition_probs", "emission_probs",
            "cluster_probs")
  for (i in pars) {
    pred[[i]]$estimate <- pred[[i]]$estimate - pred2[[i]]$estimate
  }
  if (nsim > 0) {
    for (i in pars) {
      pred$samples[[i]]$estimate <- pred$samples[[i]]$estimate - 
        pred2$samples[[i]]$estimate
    }
  }
  channel <- if (model$n_channels > 1) "channel" else NULL
  group_by_B <- switch(
    marginalize_B_over,
    "clusters" = c("time", channel, "observation"),
    "states" = c("cluster", "time", channel, "observation"),
    "sequences" = c("cluster", "time", "state", channel, "observation")
  )
  
  qs <- function(x, probs) {
    x <- quantile(x, probs)
    names(x) <- paste0("q", 100 * probs)
    as.data.frame(t(x))
  }
  out <- list()
  out$initial_probs <- cbind(
    pred$initial_probs |> 
      dplyr::group_by(cluster, state) |>
      dplyr::summarise(estimate = mean(estimate)) |> 
      dplyr::ungroup(),
    pred$samples$initial_probs |> 
      dplyr::group_by(cluster, state, replication) |>
      dplyr::summarise(estimate = mean(estimate)) |>
      dplyr::group_by(cluster, state) |>
      dplyr::summarise(qs(estimate, probs)) |>
      dplyr::ungroup() |> 
      dplyr::select(-c(cluster, state))
  )
  
  out$transition_probs <- cbind(
    pred$transition_probs |> 
      dplyr::group_by(cluster, time, state_from, state_to) |>
      dplyr::summarise(estimate = mean(estimate)) |> 
      dplyr::ungroup(),
    pred$samples$transition_probs |> 
      dplyr::group_by(cluster, time, state_from, state_to, replication) |>
      dplyr::summarise(estimate = mean(estimate)) |>
      dplyr::group_by(cluster, time, state_from, state_to) |>
      dplyr::summarise(qs(estimate, probs)) |>
      dplyr::ungroup() |> 
      dplyr::select(-c(cluster, time, state_from, state_to))
  ) |> dplyr::rename(!!time := time)
  
  out$emission_probs <- cbind(
    pred$emission_probs |> 
      dplyr::group_by(dplyr::pick(group_by_B)) |>
      dplyr::summarise(estimate = mean(estimate)) |> 
      dplyr::ungroup(),
    pred$samples$emission_probs |> 
      dplyr::group_by(dplyr::pick(c(group_by_B, "replication"))) |>
      dplyr::summarise(estimate = mean(estimate)) |>
      dplyr::group_by(dplyr::pick(group_by_B)) |>
      dplyr::summarise(qs(estimate, probs)) |>
      dplyr::ungroup() |> 
      dplyr::select(dplyr::starts_with("q"))
  ) |> dplyr::rename(!!time := time)
  
  out$cluster_probs <- cbind(
    pred$cluster_probs |> 
      dplyr::group_by(cluster) |>
      dplyr::summarise(estimate = mean(estimate)) |> 
      dplyr::ungroup(),
    pred$samples$cluster_probs |> 
      dplyr::group_by(cluster, replication) |>
      dplyr::summarise(estimate = mean(estimate)) |>
      dplyr::group_by(cluster) |>
      dplyr::summarise(qs(estimate, probs)) |>
      dplyr::ungroup() |> 
      dplyr::select(-cluster)
  )
  
  class(out) <- "amp"
  attr(out, "model") <- "mnhmm"
  attr(out, "seed") <- seed
  attr(out, "marginalize_B_over") <- marginalize_B_over
  out
}
