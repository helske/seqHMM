#' Plot colorpalettes
#' 
#' Function \code{plot_colors} plots colors and their given names (as defined by the 
#' user) for easy visualization of a colorpalette.
#' 
#' @export
#' 
#' @param x A vector of colors.
#' 
#' @seealso See e.g. the \code{\link{colorpalette}} data and \code{\link{RColorBrewer}} 
#' package for ready-made color palettes.
#'
#' @examples
#' plot_colors(colorpalette[[10]])
#' 
#' plot_colors(1:7)
#' 
#' plot_colors(c("yellow", "orange", "red", "purple", "blue", "green"))
#' 
#' plot_colors(rainbow(15))


plot_colors <- function(x) {
  
  if(!all(isColor(x))){
    stop("Please provide a vector of colors.")
  }
  par(mai=c(0.1, max(strwidth(x, "inch") + 0.4, na.rm = TRUE), 0.1, 0.4))
  barplot(rep(1, length(x)), col = rev(x), space = 0.2, axes = FALSE, 
          names.arg = rev(x), cex.names = 0.8, horiz = TRUE, las = 1)  
}