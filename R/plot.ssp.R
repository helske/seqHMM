#' Stack Multichannel Sequence Plots and/or Most Probable Paths Plots from Hidden Markov
#' Models
#'
#' Function \code{plot.ssp} plots stacked sequence plots from \code{ssp} objects defined with
#' \code{\link{ssp}}.
#'
#' @export
#' @param x An \code{ssp} object.
#' @param ... Ignored.
#'
#' @seealso \code{\link{ssp}} for more examples and information on defining the plot before using
#'   \code{plot.ssp}; \code{\link{ssplot}} for straight plotting of \code{ssp} objects;
#'   and \code{\link{gridplot}} for plotting multiple \code{ssp} objects.
#'
#' @examples
#'
#' data("biofam3c")
#'
#' ## Building sequence objects
#' child_seq <- seqdef(biofam3c$children, start = 15)
#' marr_seq <- seqdef(biofam3c$married, start = 15)
#' left_seq <- seqdef(biofam3c$left, start = 15)
#'
#' ## Choosing colors
#' attr(child_seq, "cpal") <- c("#66C2A5", "#FC8D62")
#' attr(marr_seq, "cpal") <- c("#AB82FF", "#E6AB02", "#E7298A")
#' attr(left_seq, "cpal") <- c("#A6CEE3", "#E31A1C")
#'
#'
#' # Plotting state distribution plots of observations
#' ssp1 <- ssp(list(child_seq, marr_seq, left_seq))
#' plot(ssp1)


plot.ssp <- function(x, ...) {
  plot.new()
  grid.newpage()
  savepar <- par(no.readonly = TRUE)
  do.call(SSPlotter,args = x)
  par(savepar)
}


