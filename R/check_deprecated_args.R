check_deprecated_args <- function(x) {
  if (!is.null(x$withlegend)) {
    warning("Argument `withlegend` is deprecated. Use `with.legend` instead.")
  } else if (!is.null(x$withlegend) && with.legend == "auto") {
    warning("with.legend and withlegend cannot be specified together. withlegend is deprecated, use with.legend instead.")
  }
}