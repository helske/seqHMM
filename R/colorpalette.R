#' Color palettes
#'
#' A list containing ready defined color palettes with distinct colors using
#' iWantHue. By default, `seqHMM` uses these palettes when assigning colors.
#'
#' @format A list with 200 color palettes.
#'
#' @source iWantHue web page <https://medialab.github.io/iwanthue/>
#'
#' @seealso [plot_colors()] for visualization of color palettes.
#' Implementations of iWantHue for R:
#'   \itemize{
#'     \item <https://github.com/hoesler/rwantshue>
#'     \item <https://github.com/johnbaums/hues>
#'   }
#' @docType data
#' @keywords datasets
#' @name colorpalette
#' @examples
#' data("colorpalette")
#' # Color palette with 9 colors
#' colorpalette[[9]]
#' # Color palette with 24 colors
#' colorpalette[[24]]
#'
NULL
