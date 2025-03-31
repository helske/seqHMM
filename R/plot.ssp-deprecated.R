#' Stack Multichannel Sequence Plots and/or Most Probable Paths Plots from Hidden Markov
#' Models
#'
#' Function `plot.ssp` plots stacked sequence plots from `ssp` objects defined with
#' [ssp()].
#'
#' @export
#' @param x An `ssp` object.
#' @param ... Ignored.
#'
#' @seealso [ssp()] for more examples and information on defining the plot before using
#'   `plot.ssp`; [ssplot()] for straight plotting of `ssp` objects;
#'   and [gridplot()] for plotting multiple `ssp` objects.
#'
#' @references Helske S. and Helske J. (2019). Mixture Hidden Markov Models for Sequence Data: The seqHMM Package in R,
#' Journal of Statistical Software, 88(3), 1-32. doi:10.18637/jss.v088.i03
plot.ssp <- function(x, ...) {
  .Deprecated("stacked_sequence_plot")
  plot.new()
  grid.newpage()
  savepar <- par(no.readonly = TRUE)
  do.call(SSPlotter, args = x)
  par(savepar)
}
