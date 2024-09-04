#' Visualize Average Marginal Effects
#' 
#' @importFrom ggplot2 ggplot aes geom_pointrange geom_ribbon geom_line facet_wrap
#' @param x Output from [amp()].
#' @param alpha Transparency level for [ggplot2::geom_ribbon()].
plot.amp <- function(x, type, probs = c(0.025, 0.975), alpha = 0.25) {
  type <- match.arg(type, c("initial", "transition", "emission", "cluster"))
  
  cluster <- time <- state <- state_from <- state_to <- observation <- NULL
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
  if (type == "initial") {
    p <- ggplot(x$initial, aes(estimate, state)) + 
      geom_pointrange(aes(ymin = .data[[lwr]], ymax = .data[[upr]]))
    if (!is.null(cluster)) {
      p <- p + facet_wrap(~ cluster)
    }
  }
  if (type == "transition") {
    p <- ggplot(x$transition, aes(estimate, time)) +
      geom_ribbon(
        aes(ymin = .data[[lwr]], ymax = .data[[upr]], fill = state_to),
        alpha = alpha
      ) +
      geom_line(aes(colour = state_to)) +
      facet_wrap(~ state_from)
    if (!is.null(cluster)) {
      p <- p + facet_wrap(~ cluster)
    }
  }
  if (type == "emission") {
    p <- ggplot(x$emission, aes(estimate, time)) +
      geom_ribbon(
        aes(ymin = .data[[lwr]], ymax = .data[[upr]], fill = observation),
        alpha = alpha
      ) +
      geom_line(aes(colour = observation)) +
      facet_wrap(~ state)
    if (!is.null(cluster)) {
      p <- p + facet_wrap(~ cluster)
    }
  }
  if (type == "cluster") {
    p <- ggplot(x$theta, aes(estimate, cluster)) + 
      geom_pointrange(aes(ymin = .data[[lwr]], ymax = .data[[upr]]))
  }
  
}