check_deprecated_args <- function(x) {
  if (!is.null(x$withlegend))
    stop("Argument `withlegend` is deprecated. Use `with.legend` instead.")
}