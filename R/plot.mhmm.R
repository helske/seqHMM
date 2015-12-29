#' Interactive Plotting for Mixed Hidden Markov Model (mhmm)
#'
#' Function \code{plot.mhmm} plots a directed graph of the parameters of each model
#' with pie charts of emission probabilities as vertices/nodes.
#'
#' @export
#'
#' @param x A hidden Markov model object of class \code{mhmm} created with
#'   \code{\link{build_mhmm}} (or \code{\link{build_mmm}} or
#'   \code{\link{build_lcm}}). Multichannel
#'   \code{mhmm} objects are automatically transformed into single-channel objects.
#'   See function \code{\link{mc_to_sc}} for more information on the
#'   transformation.
#' @param interactive Whether to plot each cluster in succession or in a grid.
#'   Defaults to \code{TRUE}, i.e. clusters are plotted one after another.
#'
#' @param ask If \code{TRUE} and \code{which.plots} is NULL,
#'   \code{plot.mhmm} operates in interactive mode, via \code{\link{menu}}.
#'   Defaults to \code{FALSE}. Ignored if \code{interactive = FALSE}.
#' @param which.plots The number(s) of the requested cluster(s) as an integer
#'   vector. The default \code{NULL} produces all plots.
#'
#' @param nrow,ncol Optional arguments to arrange plots in a grid. Ignored if
#'   \code{interactive = TRUE}.
#' @param byrow Controls the order of plotting in a grid. Defaults to \code{FALSE},
#'   i.e. plots are arranged column-wise. Ignored if \code{interactive = TRUE}.
#' @param row.prop Sets the proportions of the row heights of the grid. The default
#'   value is \code{"auto"} for even row heights. Takes a vector of values from
#'   0 to 1, with values summing to 1. Ignored if \code{interactive = TRUE}.
#' @param col.prop Sets the proportion of the column heights of the grid. The default
#'   value is \code{"auto"} for even column widths. Takes a vector of values
#'   from 0 to 1, with values summing to 1. Ignored if \code{interactive = TRUE}.
#'
#' @param layout specifies the layout of vertices (nodes). Accepts a
#'   numerical matrix, a \code{\link[igraph]{layout_}} function (without quotation marks),
#'   or either of the predefined options \code{"horizontal"} (the
#'   default) and \code{"vertical"}. Options \code{"horizontal"} and
#'   \code{"vertical"} position vertices at the same horizontal or vertical
#'   line. A two-column numerical matrix can be used to give x and y coordinates of
#'   the vertices. The \code{\link[igraph]{layout_}} functions available in the
#'   \code{igraph} package offer other automatic layouts for graphs.
#' @param pie Are vertices plotted as pie charts of emission probabilities?
#'   Defaults to TRUE.
#' @param vertex.size Size of vertices, given as a scalar or numerical
#'   vector. The default value is 40.
#' @param vertex.label Labels for vertices. Possible options include
#'   \code{"initial.probs"}, \code{"names"}, \code{NA}, and a character or
#'   numerical vector. The default \code{"initial.probs"} prints the initial
#'   probabilities of the model and \code{"names"} prints the names of the
#'   hidden states as labels. \code{NA} prints no labels.
#' @param vertex.label.dist Distance of the label of the vertex from its
#'   center. The default value \code{"auto"} places the label outside the
#'   vertex.
#' @param vertex.label.pos Positions of vertex labels, relative to
#'   the center of the vertex. A scalar or numerical vector giving
#'   position(s) as radians or one of \code{"bottom"} (\code{pi/2} as radians),
#'   \code{"top"} (\code{-pi/2}), \code{"left"} (\code{pi}), or
#'   \code{"right"} (\code{0}).
#' @param vertex.label.family,edge.label.family Font family to be used for
#'   vertex/edge labels. See argument \code{family} in \code{\link{par}} for
#'   more information.
#' @param loops Defines whether transitions back to same states are plotted.
#' @param edge.curved Defines whether to plot curved edges (arcs, arrows)
#'   between vertices. A logical or numerical vector or scalar. Numerical
#'   values specify curvatures of edges. The default value \code{TRUE}
#'   gives curvature of 0.5 to all edges. See \code{\link{igraph.plotting}} for
#'   more information.
#' @param edge.label Labels for edges. Possible options include
#'   \code{"auto"}, \code{NA}, and a character or numerical vector. The
#'   default \code{"auto"} prints transition probabilities as edge labels.
#'   \code{NA} prints no labels.
#' @param edge.width Width(s) for edges. The default \code{"auto"} determines
#'   widths according to transition probabilities between hidden states.
#'   Other possibilities are a scalar or a numerical vector of widths.
#' @param cex.edge.width An expansion factor for edge widths. Defaults to 1.
#' @param edge.arrow.size Size of the arrow in edges (constant). Defaults to 1.5.
#' @param label.signif Rounds labels of model parameters to specified number
#'   of significant digits, 2 by default. Ignored for user-given labels.
#' @param label.scientific Defines if scientific notation should be used to
#'   describe small numbers. Defaults to \code{FALSE}, e.g. 0.0001 instead of
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
#'   \code{"combine.slices"} is greater than 0). The default color is white.
#' @param combined.slice.label The label for combined states (when argument
#'   \code{"combine.slices"} is greater than 0) to appear in the legend.
#' @param withlegend Defines if and where the legend of state colors is
#'   plotted. Possible values include \code{"bottom"} (the default),
#'   \code{"top"}, \code{"left"}, and \code{"right"}. \code{FALSE} omits the
#'   legend.
#' @param ltext Optional description of (combined) observed states to appear
#'   in the legend. A vector of character strings. See \code{\link{seqplot}} for
#'   more information.
#' @param legend.prop Proportion used for plotting the legend. A scalar between
#'   0 and 1, defaults to 0.5.
#' @param cex.legend Expansion factor for setting the size of the font for
#'   labels in the legend. The default value is 1. Values lesser than 1 will
#'   reduce the size of the font, values greater than 1 will increase the size.
#' @param ncol.legend The number of columns for the legend. The default value
#'   \code{"auto"} sets the number of columns automatically.
#' @param cpal Optional color palette for (combinations of) observed states.
#'   The default value \code{"auto"} uses automatic color palette. Otherwise a
#'   vector of length \code{x$n_symbols} is given, i.e. the argument requires a color
#'   specified for all (combinations of) observed states even if they are not
#'   plotted (if the probability is less than \code{combine.slices}).
#' @param main Optional main titles for plots. The default \code{"auto"} uses
#' \code{cluster_names} as titles, \code{NULL} prints no titles.
#' @param ... Other parameters passed on to \code{\link{plot.igraph}} such as
#'   \code{vertex.color}, \code{vertex.label.cex}, or \code{edge.lty}.
#'
#' @seealso \code{\link{build_mhmm}} and \code{\link{fit_model}} for building and
#'   fitting mixture hidden Markov models; \code{\link{plot.igraph}} for plotting
#'   directed graphs; and \code{\link{mhmm_biofam}} and \code{\link{mhmm_mvad}} for
#'   the models used in examples.
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
#'   # Note: plotting window must be high enough!
#'   set.seed(123)
#'   plot(mhmm_mvad, interactive = FALSE,
#'     # automatic layout, legend on the right-hand side
#'     layout = layout_nicely, withlegend = "right",
#'     # Smaller and less curved edges
#'     edge.curved = 0.2, cex.edge.width = 0.5, edge.arrow.size = 0.7,
#'     vertex.label.pos = "bottom")
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
                      withlegend = "bottom", ltext = NULL, legend.prop = 0.5,
                      cex.legend = 1, ncol.legend = "auto", cpal = "auto",
                      main = "auto", ...){

  if (interactive) {
    do.call(mHMMplotint, c(list(x = x, ask = ask, which.plots = which.plots, layout = layout, pie = pie,
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
                                withlegend = withlegend, ltext = ltext, legend.prop = legend.prop,
                                cex.legend = cex.legend, ncol.legend = ncol.legend, cpal = cpal,
                                main = main), list(...)))
  } else {
    args <- as.list(match.call())[-1]
    args$ask <- args$interactive <- NULL
    do.call(mHMMplotgrid, args)
  }

}
