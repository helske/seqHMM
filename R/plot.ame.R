#' Visualize Average Marginal Effects
#' 
#' @param x Output from [ame_param() or ame_prob()].
#' @param type Type of plot to create. One of `"initial"`, `"transition"`,
#'  `"emission"`, or `"cluster"`.
#' @param probs A numeric vector of length 2 with the lower and upper limits for 
#' confidence intervals. Default is `c(0.025, 0.975)`. If the limits are not 
#' found in the input object `x`, an error is thrown.
#' @param ... Ignored.
#' @export
plot.ame <- function(x, type, probs = c(0.025, 0.975), ...) {
  
  type <- match.arg(type, c("initial", "transition", "emission", "cluster"))
  
  cluster <- time <- state <- state_from <- state_to <- observation <- 
    estimate <- NULL
  stopifnot_(
    checkmate::test_numeric(
      x = probs, lower = 0, upper = 1, any.missing = FALSE, min.len = 2L, 
      max.len = 2
    ),
    "Argument {.arg probs} must be a {.cls numeric} vector with two values
      between 0 and 1."
  )
  lwr <- paste0("q", 100 * probs[1])
  upr <- paste0("q", 100 * probs[2])
  stopifnot_(
    all(c(lwr, upr) %in% names(x$initial)),
    paste0(
      "The probabilities in {.arg probs} are not available in the {.arg x}. ",
      "Run {.fn ame} with {.arg probs} as ",
      "`c({paste(probs, collapse = ', ')}) ."
    )
  )
  if (type == "initial") {
    p <- ggplot(x$initial_probs, aes(state, estimate)) + 
      geom_pointrange(aes(ymin = .data[[lwr]], ymax = .data[[upr]]))
    if (!is.null(x$initial_probs$cluster)) {
      p <- p + facet_wrap(~ cluster)
    }
  }
  if (type == "transition") {
    p <- ggplot(x$transition_probs, aes(time, estimate)) +
      geom_pointrange(
        aes(ymin = .data[[lwr]], ymax = .data[[upr]], colour = state_to)
      )
    if (is.null(x$transition_probs$cluster)) {
      p <- p + facet_wrap(~ state_from)
    } else {
      p <- p + facet_wrap(~ state_from + cluster)
    }
  }
  if (type == "emission") {
    p <- ggplot(x$emission_probs, aes(time, estimate)) +
      geom_pointrange(
        aes(ymin = .data[[lwr]], ymax = .data[[upr]], colour = observation)
      )
    if (is.null(x$emission_probs$cluster)) {
      p <- p + facet_wrap(~ state)
    } else {
      p <- p + facet_wrap(~ state + cluster)
    }
  }
  if (type == "cluster") {
    p <- ggplot(x$cluster_probs, aes(cluster, estimate)) + 
      geom_pointrange(aes(ymin = .data[[lwr]], ymax = .data[[upr]]))
  }
  
}