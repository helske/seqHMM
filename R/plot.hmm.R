#' Plot hidden Markov models
#'
#' Function \code{plot.hmm} plots a directed graph with pie charts of
#' emission probabilities as vertices/nodes.
#'
#' @export
#' @param x A hidden Markov model object of class \code{hmm} created with
#'   \code{\link{build_hmm}} (or \code{\link{build_mm}}). Multichannel
#'   \code{hmm} objects are automatically transformed into single-channel objects.
#'   See function \code{\link{mc_to_sc}} for more information on the
#'   transformation.
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
#' @param main Main title for the plot. Omitted by default.
#' @param ... Other parameters passed on to \code{\link{plot.igraph}} such as
#'   \code{vertex.color}, \code{vertex.label.cex}, or \code{edge.lty}.
#'
#' @seealso \code{\link{build_hmm}} and \code{\link{fit_model}} for building and
#'   fitting Hidden Markov models, \code{\link{mc_to_sc}} for transforming
#'   multistate \code{hmm} objects into single-channel objects,
#'   \code{\link{hmm_biofam}} and \code{\link{hmm_mvad}} for information on the models
#'   used in the examples, and
#'   \code{\link{plot.igraph}} for the general plotting function of directed graphs.
#'
#' @examples
#' # Multichannel data, left-to-right model
#'
#' # Loading a HMM of the biofam data
#' data("hmm_biofam")
#'
#' # Plotting hmm object
#' plot(hmm_biofam)
#'
#' # Plotting HMM with
#' plot(hmm_biofam,
#'   # varying curvature of edges
#'   edge.curved = c(0, -0.7, 0.6, 0.7, 0, -0.7, 0),
#'   # legend with two columns and less space
#'   ncol.legend = 2, legend.prop = 0.4,
#'   # new label for combined slice
#'   combined.slice.label = "States with probability < 0.05")
#'
#' # Plotting HMM with given coordinates
#' plot(hmm_biofam,
#'   # layout given in 2x5 matrix
#'   # x coordinates in the first column
#'   # y coordinates in the second column
#'   layout = matrix(c(1, 3, 3, 5,  3,
#'                     0, 0, 1, 0, -1), ncol = 2),
#'   # larger vertices
#'   vertex.size = 50,
#'   # straight edges
#'   edge.curved = FALSE,
#'   # thinner edges and arrows
#'   cex.edge.width = 0.5, edge.arrow.size = 1,
#'   # varying positions for vertex labels (initial probabilities)
#'   vertex.label.pos = c(pi, pi/2, -pi/2, 0, pi/2),
#'   # different legend properties
#'   withlegend = "top", legend.prop = 0.3, cex.legend = 1.1,
#'   # Fix axes to the right scale
#'   xlim = c(0.5, 5.5), ylim = c(-1.5, 1.5), rescale = FALSE,
#'   # all states (not combining states with small probabilities)
#'   combine.slices = 0,
#'   # legend with two columns
#'   ncol.legend = 2)
#'
#' # Plotting HMM with own color palette
#' plot(hmm_biofam, cpal = 1:10,
#'   # States with emission probability less than 0.2 removed
#'   combine.slices = 0.2,
#'   # legend with two columns
#'   ncol.legend = 2)
#'
#' # Plotting HMM without pie graph and with a layout function
#' require("igraph")
#' # Setting the seed for a random layout
#' set.seed(1234)
#' plot(hmm_biofam,
#'   # Without pie graph
#'   pie = FALSE,
#'   # Using an automatic layout function from igraph
#'   layout = layout_nicely,
#'   vertex.size = 30,
#'   # Straight edges and probabilities of moving to the same state
#'   edge.curved = FALSE, loops = TRUE,
#'   # Labels with three significant digits
#'   label.signif = 3,
#'   # Fixed edge width
#'   edge.width = 1,
#'   # Remove edges with probability less than 0.01
#'   trim = 0.01,
#'   # Hidden state names as vertex labels
#'   vertex.label = "names",
#'   # Labels insidde vertices
#'   vertex.label.dist = 0,
#'   # Fix x-axis (more space on the right-hand side)
#'   xlim = c(-1, 1.3))
#'
#'
#' # Single-channel data, unrestricted model
#'
#' # Loading a hidden Markov model of the mvad data (hmm object)
#' data("hmm_mvad")
#'
#' # Plotting HMM
#' plot(hmm_mvad)
#'
#' require("igraph")
#' plot(hmm_mvad,
#'   # Layout in circle (layout function from igraph)
#'   layout = layout_in_circle,
#'   # Less curved edges with smaller arrows, no labels
#'   edge.curved = 0.2, edge.arrow.size = 0.9, edge.label = NA,
#'   # Positioning vertex labels (initial probabilities)
#'   vertex.label.pos = c("right", "right", "left", "left", "right"),
#'   # Less space for the legend
#'   legend.prop = 0.3)



plot.hmm  <- function(x, layout = "horizontal", pie = TRUE,
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
                      main = NULL, ...){


  # Saving and changing marginals
  oldPar  <- par(no.readonly = TRUE)

  if (is.null(main)) {
    par(mar = c(0.5, 0.5, 0.5, 0.5))
  } else {
    par(mai = c(0, 0, 1, 0))
  }
  on.exit(par(oldPar), add = TRUE)
  on.exit(par(mfrow = c(1, 1)), add = TRUE)

  dots  <- list(...)

  labelprint  <- function(z, labs) {
    if (labs == TRUE && (z > 0.001 || z == 0)) {
      labs  <- FALSE
    }
    if (z < 10 ^ -(label.max.length)) {
      z  <- prettyNum(signif(round(z, digits = label.max.length), digits = label.signif), scientific = labs)
    } else {
      z  <- prettyNum(signif(z, digits = label.signif), scientific = labs)
    }
  }

  if (!is.matrix(layout) && !is.function(layout)) {
    if (!(layout %in% c("horizontal", "vertical"))) {
      stop("Argument layout only accepts numerical matrices, igraph layout functions, or strings \"horizontal\" and \"vertical\".")
    }
  }

  if (!is.numeric(vertex.label.pos)) {
    choices <- c("bottom", "top", "left", "right")
    ind <- pmatch(vertex.label.pos, choices, duplicates.ok = TRUE)
    if (any(is.na(ind))) {
      stop("Argument vertex.label.pos only accepts values \"bottom\", \"top\", \"left\", \"right\" or a numerical vector.")
    }
    vertex.label.pos <- choices[ind]
  }

  choices <- c(TRUE, FALSE, "bottom", "top", "left", "right")
  ind <- pmatch(withlegend, choices)
  if (is.na(ind)) {
    stop("Argument withlegend must be one of TRUE, FALSE, \"bottom\", \"right\", \"top\", or \"left\".")
  }
  withlegend <- choices[ind]
  if (withlegend %in% c(TRUE, "auto")){
    withlegend <- "bottom"
  }


  # Convert multichannel models to single-channel
  if (x$n_channels > 1) {
    x  <- mc_to_sc(x)
  }

  # No slices -> no legends needed
  if (pie == FALSE && withlegend != FALSE) {
    withlegend  <- FALSE
  }


  # Positions of vertex labels
  if (!is.numeric(vertex.label.pos)) {
    vpos <- numeric(length(vertex.label.pos))
    for(i in 1:length(vertex.label.pos)){
      if (vertex.label.pos[i] == "bottom") {
        vpos[i]  <- pi/2
      } else if(vertex.label.pos[i] == "top") {
        vpos[i]  <- -pi/2
      } else if(vertex.label.pos[i] == "left") {
        vpos[i]  <- pi
      } else {
        vpos[i]  <- 0
      }
    }
    vertex.label.pos <- vpos
  }

  # Vertex labels
  if (length(vertex.label) == 1 && !is.na(vertex.label) && vertex.label != FALSE) {
    if (vertex.label == "initial.probs") {
      vertex.label  <- sapply(x$initial_probs, labelprint, labs = label.scientific)
    } else if(vertex.label == "names") {
      vertex.label  <- x$state_names
    }
  } else if (length(vertex.label) != length(x$state_names)) {
    warning("The length of the vector provided for the argument \"vertex.label\" does not match the number of hidden states.")
    vertex.label  <- rep(vertex.label, length.out = length(x$state_names))
  }

  # Vertex label distances
  if (is.character(vertex.label.dist)) {
    ind <- pmatch(vertex.label.dist, "auto")
    if (is.na(ind)) {
      stop("Argument vertex.label.dist only accepts the value \"auto\" or a numerical vector.")
    }
    vertex.label.dist  <- vertex.size * 0.4 / 10
  } else if (length(vertex.label.dist) > 1 && length(vertex.label.dist) != x$n_states){
    warning("The length of the vector provided for the argument \"vertex.label.dist\" does not match the number of edges.")
    vertex.label.dist  <- rep(vertex.label.dist, length.out = length(x$n_states))
  }


  # Trimming (remove small transition probablities from plot)
  transM  <- x$transition_probs
  transM[transM < trim]  <- 0

  # Adjacency matrix (which edges to plot)
  edges  <- transM
  edges[edges > 0]  <- 1
  # Remove transitions back to the same state
  if(!loops){
    diag(edges)  <- 0
  }

  # Vector of non-zero transition probabilities
  transitions  <- transM
  if (loops == FALSE && length(transitions) > 1) {
    diag(transitions)  <- 0
  }
  transitions  <- t(transitions)[t(transitions) > 0]

  # Edge labels
  if (!is.na(edge.label) && edge.label != FALSE) {
    if (length(edge.label) == 1 && (edge.label == "auto" || edge.label == TRUE)) {
      edge.label  <- sapply(transitions, labelprint, labs = label.scientific)
    } else if (length(edge.label) > 1 && length(edge.label) != length(transitions)) {
      warning("The length of the vector provided for the argument \"edge.label\" does not match the number of edges.")
      edge.label  <- rep(edge.label, length.out = length(transitions))
    }
  }


  # Edge widths
  if (is.character(edge.width)) {
    ind <- pmatch(edge.width, "auto")
    if (is.na(ind)) {
      stop("Argument edge.width only accepts the value \"auto\" or a numerical vector.")
    }
    edge.width  <- transitions * (7 / max(transitions)) * cex.edge.width
  } else if (length(edge.width) > 1 && edge.width != length(transitions)) {
    warning("The length of the vector provided for the argument \"edge.width\" does not match the number of edges.")
    edge.width  <- rep(edge.width, length.out = length(transitions))
  }

  # Defining the graph structure
  g1  <- graph.adjacency(edges, mode = "directed")

  # Layout of the graph
  if (is.function(layout)) {
    glayout  <- layout(g1)
  } else if (is.matrix(layout)) {
    glayout  <- layout
  } else {
    if (layout == "horizontal") {
      glayout  <- layout_on_grid(g1, width = x$n_states)
    } else if(layout == "vertical") {
      glayout  <- layout_on_grid(g1, width = 1)
    }
  }


  # Colors for the (combinations of) observed states
  if (length(cpal) == 1 && cpal == "auto") {
    pie.colors  <- attr(x$observations, "cpal")
  } else if(length(cpal) != ncol(x$emiss)) {
    warning("The length of the vector provided for argument cpal does not match the number of observed states. Automatic color palette was used.")
    pie.colors  <- attr(x$observations, "cpal")
  } else if(!all(isColor(cpal))) {
    stop(paste("Please provide a vector of colors for argument cpal or use value \"auto\" for automatic color palette."))
  } else {
    pie.colors  <- cpal
  }
  if (withlegend != FALSE) {
    pie.colors.l  <- pie.colors
  }

  # Legend position and number of columns
  if (withlegend != FALSE && pie == TRUE) {
    if (!is.null(ltext)) {
      if (length(ltext) != x$n_symbols) {
        stop("The length of the argument ltext does not match the number of (combined) observed states.")
      }
      # ltext = NULL
    } else {
      ltext  <- x$symbol_names
    }
    if (withlegend == "bottom") {
      graphics::layout(matrix(1:2, nrow = 2), heights = c(1 - legend.prop, legend.prop))
    } else if (withlegend == "right") {
      graphics::layout(matrix(1:2, nrow = 1), widths = c(1 - legend.prop, legend.prop))
    } else if (withlegend == "left") {
      graphics::layout(matrix(2:1, nrow = 1), widths = c(legend.prop, 1 - legend.prop))
    } else {
      graphics::layout(matrix(2:1, nrow = 2), widths = c(legend.prop, 1 - legend.prop))
    }
    par(cex = 1)
  }

  # Defining rescale, xlim, ylim if not given
  if (!is.matrix(layout) && !is.function(layout)) {
    if (layout == "horizontal") {
      if (hasArg(rescale)) {
        rescale  <- dots$rescale
      } else {
        rescale <- FALSE
      }
      if (hasArg(xlim)) {
        xlim  <- dots$xlim
      } else {
        if (rescale == TRUE) {
          xlim  <- c(-1, 1)
        } else {
          xlim  <- c(-0.1, ncol(transM) - 1 + 0.1)
        }
      }
      if (hasArg(ylim)) {
        ylim  <- dots$ylim
      } else {
        if (rescale == TRUE) {
          ylim  <- c(-1, 1)
        } else {
          ylim  <- c(-0.5, 0.5)
        }
      }
      dots[["xlim"]]  <- NULL
      dots[["ylim"]]  <- NULL
      dots[["rescale"]]  <- NULL
    } else if (layout == "vertical") {
      if (hasArg(rescale)) {
        rescale  <- dots$rescale
      } else {
        rescale <- FALSE
      }
      if (hasArg(xlim)) {
        xlim  <- dots$xlim
      } else {
        if (rescale == TRUE) {
          xlim  <- c(-1, 1)
        } else {
          xlim  <- c(-0.5, 0.5)
        }
      }
      if (hasArg(ylim)) {
        ylim  <- dots$ylim
      } else {
        if (rescale == TRUE) {
          ylim  <- c(-1, 1)
        } else {
          ylim  <- c(-0.1, ncol(transM) - 1 + 0.1)
        }
      }
      dots[["xlim"]]  <- NULL
      dots[["ylim"]]  <- NULL
      dots[["rescale"]]  <- NULL
    }
  }


  # Plotting graph
  if (pie == TRUE) {
    pie.values  <- lapply(seq_len(nrow(transM)), function(i) x$emission_probs[i,])
    # If slices are combined
    if (combine.slices > 0 &&
        !all(unlist(pie.values)[unlist(pie.values) > 0] > combine.slices)) {
      if (withlegend != FALSE) {
        pie.colors.l  <- NULL
        lt  <- NULL
      }
      for (i in 1:x$n_states) {
        # How much probability for combined slice
        cs.prob  <- sum(pie.values[[i]][pie.values[[i]] < combine.slices])
        # Remove small probabilities form pies
        pie.values[[i]][pie.values[[i]] < combine.slices]  <- 0
        # Colors and labels for large slices
        pie.values[[i]]  <- c(pie.values[[i]], cs.prob)
        # Texts and colors for legends
        if (withlegend != FALSE) {
          pie.colors.l  <- c(pie.colors.l, pie.colors[pie.values[[i]][1:(length(pie.values[[i]]) - 1)] >= combine.slices])
          lt  <- c(lt, ltext[pie.values[[i]][1:(length(pie.values[[i]]) - 1)] >= combine.slices])
        }
      }
      if (withlegend != FALSE) {
        ltext  <- c(unique(lt), combined.slice.label)
        pie.colors.l  <- c(unique(pie.colors.l), combined.slice.color)
      }
      if (ncol.legend == "auto") {
        if (withlegend == "bottom" || withlegend == "top") {
          ncol.legend  <- ceiling(length(pie.colors) / 4)
        } else {
          ncol.legend  <- 1
        }
      }
      pie.colors <- c(pie.colors, combined.slice.color)
      # Slices not combined
    } else {
      if (ncol.legend == "auto") {
        if (withlegend == "bottom" || withlegend == "top") {
          ncol.legend  <- ceiling(ncol(x$emission_probs) / 4)
        } else {
          ncol.legend  <- 1
        }
      }
    }

    if (!is.matrix(layout) && !is.function(layout) &&
        (layout == "horizontal" || layout == "vertical")) {
      do.call(plot.igraph2, c(list(g1, layout = glayout,
                                   vertex.shape = "pie", vertex.pie = pie.values,
                                   vertex.pie.color = list(pie.colors),
                                   vertex.size = vertex.size,
                                   vertex.label = vertex.label, vertex.label.dist = vertex.label.dist,
                                   vertex.label.degree = vertex.label.pos,
                                   vertex.label.family = vertex.label.family,
                                   edge.curved = edge.curved, edge.width = edge.width,
                                   edge.label = edge.label,
                                   edge.label.family = edge.label.family,
                                   edge.arrow.size = edge.arrow.size,
                                   xlim = xlim, ylim = ylim, rescale = rescale,
                                   main = main), dots))
    } else {
      do.call(plot.igraph2, c(list(g1, layout = glayout,
                                   vertex.shape = "pie", vertex.pie = pie.values,
                                   vertex.pie.color = list(pie.colors),
                                   vertex.size = vertex.size,
                                   vertex.label = vertex.label, vertex.label.dist = vertex.label.dist,
                                   vertex.label.degree = vertex.label.pos,
                                   vertex.label.family = vertex.label.family,
                                   edge.curved = edge.curved, edge.width = edge.width,
                                   edge.label = edge.label,
                                   edge.label.family = edge.label.family,
                                   edge.arrow.size = edge.arrow.size,
                                   main = main), dots))
    }
  # pie = FALSE
  } else {
    if(!is.matrix(layout) && !is.function(layout) && (layout == "horizontal" || layout == "vertical")){
      do.call(plot.igraph2, c(list(g1, layout = glayout,
                                   vertex.size = vertex.size,
                                   vertex.label = vertex.label, vertex.label.dist = vertex.label.dist,
                                   vertex.label.degree = vertex.label.pos,
                                   vertex.label.family = vertex.label.family,
                                   edge.curved = edge.curved, edge.width = edge.width,
                                   edge.label = edge.label,
                                   edge.label.family = edge.label.family,
                                   xlim = xlim, ylim = ylim, rescale = rescale,
                                   main = main), dots))
    } else {
      do.call(plot.igraph2, c(list(g1, layout = glayout,
                                   vertex.size = vertex.size,
                                   vertex.label = vertex.label, vertex.label.dist = vertex.label.dist,
                                   vertex.label.degree = vertex.label.pos,
                                   vertex.label.family = vertex.label.family,
                                   edge.curved = edge.curved, edge.width = edge.width,
                                   edge.label = edge.label,
                                   edge.label.family = edge.label.family,
                                   main = main), dots))
    }
  }


  # Plotting legend
  if (withlegend != FALSE && pie == TRUE) {
    seqlegend(x$observations, cpal = pie.colors.l, ltext = ltext,
              position = "center", fontsize = cex.legend, ncol = ncol.legend,
              with.missing = FALSE)
  }

  par(mfrow = c(1, 1))
}
