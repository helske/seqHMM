#' Interactive Plotting for Mixed Hidden Markov Model (mhmm)
#'
#' Function `plot.mhmm` plots a directed graph of the parameters of each model
#' with pie charts of emission probabilities as vertices/nodes.
#'
#' @export
#'
#' @param x A hidden Markov model object of class `mhmm` created with
#'   [build_mhmm()] (or [build_mmm()] or
#'   [build_lcm()]). Multichannel
#'   `mhmm` objects are automatically transformed into single-channel objects.
#'   See function [mc_to_sc()] for more information on the
#'   transformation.
#' @param interactive Whether to plot each cluster in succession or in a grid.
#'   Defaults to `TRUE`, i.e. clusters are plotted one after another.
#'
#' @param ask If `TRUE` and `which.plots` is NULL,
#'   `plot.mhmm` operates in interactive mode, via [utils::menu()].
#'   Defaults to `FALSE`. Ignored if `interactive = FALSE`.
#' @param which.plots The number(s) of the requested cluster(s) as an integer
#'   vector. The default `NULL` produces all plots.
#'
#' @param nrow,ncol Optional arguments to arrange plots in a grid. Ignored if
#'   `interactive = TRUE`.
#' @param byrow Controls the order of plotting in a grid. Defaults to `FALSE`,
#'   i.e. plots are arranged column-wise. Ignored if `interactive = TRUE`.
#' @param row.prop Sets the proportions of the row heights of the grid. The default
#'   value is `"auto"` for even row heights. Takes a vector of values from
#'   0 to 1, with values summing to 1. Ignored if `interactive = TRUE`.
#' @param col.prop Sets the proportion of the column heights of the grid. The default
#'   value is `"auto"` for even column widths. Takes a vector of values
#'   from 0 to 1, with values summing to 1. Ignored if `interactive = TRUE`.
#'
#' @param layout specifies the layout of vertices (nodes). Accepts a
#'   numerical matrix, a [igraph::layout_()] function (without quotation marks),
#'   or either of the predefined options `"horizontal"` (the
#'   default) and `"vertical"`. Options `"horizontal"` and
#'   `"vertical"` position vertices at the same horizontal or vertical
#'   line. A two-column numerical matrix can be used to give x and y coordinates of
#'   the vertices. The [igraph::layout_()] functions available in the
#'   `igraph` package offer other automatic layouts for graphs.
#' @param pie Are vertices plotted as pie charts of emission probabilities?
#'   Defaults to TRUE.
#' @param vertex.size Size of vertices, given as a scalar or numerical
#'   vector. The default value is 40.
#' @param vertex.label Labels for vertices. Possible options include
#'   `"initial.probs"`, `"names"`, `NA`, and a character or
#'   numerical vector. The default `"initial.probs"` prints the initial
#'   probabilities of the model and `"names"` prints the names of the
#'   hidden states as labels. `NA` prints no labels.
#' @param vertex.label.dist Distance of the label of the vertex from its
#'   center. The default value `"auto"` places the label outside the
#'   vertex.
#' @param vertex.label.pos Positions of vertex labels, relative to
#'   the center of the vertex. A scalar or numerical vector giving
#'   position(s) as radians or one of `"bottom"` (`pi/2` as radians),
#'   `"top"` (`-pi/2`), `"left"` (`pi`), or
#'   `"right"` (`0`).
#' @param vertex.label.family,edge.label.family Font family to be used for
#'   vertex/edge labels. See argument `family` in [par()] for
#'   more information.
#' @param loops Defines whether transitions back to same states are plotted.
#' @param edge.curved Defines whether to plot curved edges (arcs, arrows)
#'   between vertices. A logical or numerical vector or scalar. Numerical
#'   values specify curvatures of edges. The default value `TRUE`
#'   gives curvature of 0.5 to all edges. See [igraph::igraph.plotting()] for
#'   more information.
#' @param edge.label Labels for edges. Possible options include
#'   `"auto"`, `NA`, and a character or numerical vector. The
#'   default `"auto"` prints transition probabilities as edge labels.
#'   `NA` prints no labels.
#' @param edge.width Width(s) for edges. The default `"auto"` determines
#'   widths according to transition probabilities between hidden states.
#'   Other possibilities are a scalar or a numerical vector of widths.
#' @param cex.edge.width An expansion factor for edge widths. Defaults to 1.
#' @param edge.arrow.size Size of the arrow in edges (constant). Defaults to 1.5.
#' @param label.signif Rounds labels of model parameters to specified number
#'   of significant digits, 2 by default. Ignored for user-given labels.
#' @param label.scientific Defines if scientific notation should be used to
#'   describe small numbers. Defaults to `FALSE`, e.g. 0.0001 instead of
#'   1e-04. Ignored for user-given labels.
#' @param label.max.length Maximum number of digits in labels of model
#'   parameters. Ignored for user-given labels.
#' @param trim Scalar between 0 and 1 giving the highest probability of
#'   transitions that are plotted as edges, defaults to 1e-15.
#' @param combine.slices Scalar between 0 and 1 giving the highest probability
#'   of emission probabilities that are combined into one state. The dafault
#'   value is 0.05.
#' @param combined.slice.color Color of the combined slice that includes
#'   the smallest emission probabilities (only if argument
#'   `"combine.slices"` is greater than 0). The default color is white.
#' @param combined.slice.label The label for combined states (when argument
#'   `"combine.slices"` is greater than 0) to appear in the legend.
#' @param with.legend Defines if and where the legend of state colors is
#'   plotted. Possible values include `"bottom"` (the default),
#'   `"top"`, `"left"`, and `"right"`. `FALSE` omits the
#'   legend.
#' @param ltext Optional description of (combined) observed states to appear
#'   in the legend. A vector of character strings. See [TraMineR::seqplot()] for
#'   more information.
#' @param legend.prop Proportion used for plotting the legend. A scalar between
#'   0 and 1, defaults to 0.5.
#' @param cex.legend Expansion factor for setting the size of the font for
#'   labels in the legend. The default value is 1. Values lesser than 1 will
#'   reduce the size of the font, values greater than 1 will increase the size.
#' @param ncol.legend The number of columns for the legend. The default value
#'   `"auto"` sets the number of columns automatically.
#' @param cpal Optional color palette for (combinations of) observed states.
#'   The default value `"auto"` uses automatic color palette. Otherwise a
#'   vector of length `x$n_symbols` is given, i.e. the argument requires a color
#'   specified for all (combinations of) observed states even if they are not
#'   plotted (if the probability is less than `combine.slices`).
#' @param main Optional main titles for plots. The default `"auto"` uses
#' `cluster_names` as titles, `NULL` prints no titles.
#' @param withlegend Deprecated. Use `with.legend` instead.
#' @param ... Other parameters passed on to [igraph::plot.igraph()] such as
#'   `vertex.color`, `vertex.label.cex`, or `edge.lty`.
#'
#' @seealso [build_mhmm()] and [fit_model()] for building and
#'   fitting mixture hidden Markov models; [igraph::plot.igraph()] for plotting
#'   directed graphs; and [mhmm_biofam()] and [mhmm_mvad()] for
#'   the models used in examples.
#'
#' @references Helske S. and Helske J. (2019). Mixture Hidden Markov Models for Sequence Data: The seqHMM Package in R,
#' Journal of Statistical Software, 88(3), 1-32. doi:10.18637/jss.v088.i03
#'
#' @examples
#' # Loading mixture hidden Markov model (mhmm object)
#' # of the biofam data
#' data("mhmm_biofam")
#'
#' # Plotting only the first cluster
#' plot(mhmm_biofam, which.plots = 1)
#'
#' if (interactive()) {
#'   # Plotting each cluster (change with Enter)
#'   plot(mhmm_biofam)
#'
#'   # Choosing the cluster (one at a time)
#'   plot(mhmm_biofam, ask = TRUE)
#'
#'   # Loading MHMM of the mvad data
#'   data("mhmm_mvad")
#'
#'   # Plotting models in the same graph (in a grid)
#'   # Note: the plotting window must be high enough!
#'   set.seed(123)
#'   plot(mhmm_mvad,
#'     interactive = FALSE,
#'     # automatic layout, legend on the right-hand side
#'     layout = layout_nicely, with.legend = "right",
#'     # Smaller and less curved edges
#'     edge.curved = 0.2, cex.edge.width = 0.5, edge.arrow.size = 0.7,
#'     vertex.label.pos = -4 * pi / 5, vertex.label.dist = 5
#'   )
#' }
#'
plot.mhmm <- function(x, interactive = TRUE,
                      ask = FALSE, which.plots = NULL,
                      nrow = NA, ncol = NA, byrow = FALSE,
                      row.prop = "auto", col.prop = "auto",
                      layout = "horizontal", pie = TRUE,
                      vertex.size = 40, vertex.label = "initial.probs",
                      vertex.label.dist = "auto", vertex.label.pos = "bottom",
                      vertex.label.family = "sans",
                      loops = FALSE, edge.curved = TRUE, edge.label = "auto",
                      edge.width = "auto", cex.edge.width = 1,
                      edge.arrow.size = 1.5, edge.label.family = "sans",
                      label.signif = 2, label.scientific = FALSE, label.max.length = 6,
                      trim = 1e-15,
                      combine.slices = 0.05, combined.slice.color = "white",
                      combined.slice.label = "others",
                      with.legend = "bottom", ltext = NULL, legend.prop = 0.5,
                      cex.legend = 1, ncol.legend = "auto", cpal = "auto",
                      main = "auto", withlegend, ...) {

  if (interactive) {
    do.call(mHMMplotint, c(list(
      x = x, ask = ask, which.plots = which.plots, layout = layout, pie = pie,
      vertex.size = vertex.size, vertex.label = vertex.label,
      vertex.label.dist = vertex.label.dist, vertex.label.pos = vertex.label.pos,
      vertex.label.family = vertex.label.family,
      loops = loops, edge.curved = edge.curved, edge.label = edge.label,
      edge.width = edge.width, cex.edge.width = cex.edge.width,
      edge.arrow.size = edge.arrow.size, edge.label.family = edge.label.family,
      label.signif = label.signif, label.scientific = label.scientific,
      label.max.length = label.max.length,
      trim = trim,
      combine.slices = combine.slices, combined.slice.color = combined.slice.color,
      combined.slice.label = combined.slice.label,
      with.legend = with.legend, ltext = ltext, legend.prop = legend.prop,
      cex.legend = cex.legend, ncol.legend = ncol.legend, cpal = cpal,
      main = main
    ), list(...)))
  } else {
    args <- as.list(match.call())[-1]
    args$ask <- args$interactive <- NULL
    args$x <- x
    do.call(mHMMplotgrid, args)
  }
}
