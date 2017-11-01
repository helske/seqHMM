check_deprecated_args <- function(x) {
  if (!is.null(x$withlegend))
    warning("Argument `withlegend` is deprecated. Use `with.legend` instead. ")
}