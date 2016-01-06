# Used in plot.mhmm through mHMMplotgrid

HMMplot <- function(x, layout = "horizontal", pie = TRUE,
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
                    legend.pos = "center", main = "auto", ...){


  dots <- list(...)

  labelprint <- function(z, labs) {
    if (labs == TRUE && (z > 0.001 || z == 0)) {
      labs <- FALSE
    }
    if (z < 10^-(label.max.length)) {
      z <- prettyNum(signif(round(z, digits = label.max.length), digits = label.signif), scientific = labs)
    } else {
      z <- prettyNum(signif(z, digits = label.signif), scientific = labs)
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
    x <- mc_to_sc(x)
  }

  # No slices -> no legends needed
  if (pie == FALSE && withlegend != FALSE) {
    withlegend <- FALSE
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
      vertex.label <- sapply(x$initial_probs, labelprint, labs = label.scientific)
    } else if(vertex.label == "names") {
      vertex.label <- x$state_names
    }
  } else if (length(vertex.label) != length(x$state_names)) {
    warning("The length of the vector provided for the argument \"vertex.label\" does not match the number of hidden states.")
    vertex.label <- rep(vertex.label, length.out = length(x$state_names))
  }

  # Vertex label distances
  if (is.character(vertex.label.dist)) {
    ind <- pmatch(vertex.label.dist, "auto")
    if (is.na(ind)) {
      stop("Argument vertex.label.dist only accepts the value \"auto\" or a numerical vector.")
    }
    vertex.label.dist <- vertex.size * 0.4 / 10
  } else if (length(vertex.label.dist) > 1 && length(vertex.label.dist) != x$n_states){
    warning("The length of the vector provided for the argument \"vertex.label.dist\" does not match the number of edges.")
    vertex.label.dist <- rep(vertex.label.dist, length.out = length(x$n_states))
  }


  # Trimming (remove small transition probablities from plot)
  transM <- x$transition_probs
  transM[transM < trim] <- 0

  # Adjacency matrix (which edges to plot)
  edges <- transM
  edges[edges > 0] <- 1
  # Remove transitions back to the same state
  if(!loops){
    diag(edges) <- 0
  }

  # Vector of non-zero transition probabilities
  transitions <- transM
  if (loops == FALSE && length(transitions) > 1) {
    diag(transitions) <- 0
  }
  transitions <- t(transitions)[t(transitions) > 0]

  # Edge labels
  if (!is.na(edge.label) && edge.label != FALSE) {
    if (length(edge.label) == 1 && (edge.label == "auto" || edge.label == TRUE)) {
      edge.label <- sapply(transitions, labelprint, labs = label.scientific)
    } else if (length(edge.label) > 1 && length(edge.label) != length(transitions)) {
      warning("The length of the vector provided for the argument \"edge.label\" does not match the number of edges.")
      edge.label <- rep(edge.label, length.out = length(transitions))
    }
  }


  # Edge widths
  if (is.character(edge.width)) {
    ind <- pmatch(edge.width, "auto")
    if (is.na(ind)) {
      stop("Argument edge.width only accepts the value \"auto\" or a numerical vector.")
    }
    edge.width <- transitions * (7 / max(transitions)) * cex.edge.width
  } else if (length(edge.width) > 1 && edge.width != length(transitions)) {
    warning("The length of the vector provided for the argument \"edge.width\" does not match the number of edges.")
    edge.width <- rep(edge.width, length.out = length(transitions))
  }

  # Defining the graph structure
  g1 <- graph.adjacency(edges, mode = "directed")

  # Layout of the graph
  if (is.function(layout)) {
    glayout <- layout(g1)
  } else if (is.matrix(layout)) {
    glayout <- layout
  } else {
    if (layout == "horizontal") {
      glayout <- layout_on_grid(g1, width = x$n_states)
    } else if(layout == "vertical") {
      glayout <- layout_on_grid(g1, width = 1)
    }
  }


  # Colors for the (combinations of) observed states
  if (length(cpal) == 1 && cpal == "auto") {
    pie.colors <- attr(x$observations, "cpal")
  } else if(length(cpal) != ncol(x$emiss)) {
    warning("The length of the vector provided for argument cpal does not match the number of observed states. Automatic color palette was used.")
    pie.colors <- attr(x$observations, "cpal")
  } else if(!all(isColor(cpal))) {
    stop(paste("Please provide a vector of colors for argument cpal or use value \"auto\" for automatic color palette."))
  } else {
    pie.colors <- cpal
  }
  if (withlegend != FALSE) {
    pie.colors.l <- pie.colors
  }

  # Legend position and number of columns
  if (withlegend != FALSE && pie == TRUE) {
    if (!is.null(ltext)) {
      if (length(ltext) != x$n_symbols) {
        stop("The length of the argument ltext does not match the number of (combined) observed states.")
      }
      # ltext = NULL
    } else {
      ltext <- x$symbol_names
    }
  }

  # Defining rescale, xlim, ylim if not given
  if (!is.matrix(layout) && !is.function(layout)) {
    if (layout == "horizontal") {
      if (hasArg(rescale)) {
        rescale <- dots$rescale
      } else {
        rescale <- FALSE
      }
      if (hasArg(xlim)) {
        xlim <- dots$xlim
      } else {
        if (rescale == TRUE) {
          xlim <- c(-1, 1)
        } else {
          xlim <- c(-0.1, ncol(transM) - 1 + 0.1)
        }
      }
      if (hasArg(ylim)) {
        ylim <- dots$ylim
      } else {
        if (rescale == TRUE) {
          ylim <- c(-1, 1)
        } else {
          ylim <- c(-0.5, 0.5)
        }
      }
      dots[["xlim"]] <- NULL
      dots[["ylim"]] <- NULL
      dots[["rescale"]] <- NULL
    } else if (layout == "vertical") {
      if (hasArg(rescale)) {
        rescale <- dots$rescale
      } else {
        rescale <- FALSE
      }
      if (hasArg(xlim)) {
        xlim <- dots$xlim
      } else {
        if (rescale == TRUE) {
          xlim <- c(-1, 1)
        } else {
          xlim <- c(-0.5, 0.5)
        }
      }
      if (hasArg(ylim)) {
        ylim <- dots$ylim
      } else {
        if (rescale == TRUE) {
          ylim <- c(-1, 1)
        } else {
          ylim <- c(-0.1, ncol(transM) - 1 + 0.1)
        }
      }
      dots[["xlim"]] <- NULL
      dots[["ylim"]] <- NULL
      dots[["rescale"]] <- NULL
    }
  }


  # Plotting graph
  if (pie == TRUE) {
    pie.values <- lapply(seq_len(nrow(transM)), function(i) x$emission_probs[i,])
    # If slices are combined
    if (combine.slices > 0 &&
        !all(unlist(pie.values)[unlist(pie.values) > 0] > combine.slices)) {
      if (withlegend != FALSE) {
        pie.colors.l <- NULL
        lt <- NULL
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
      pie.colors <- c(pie.colors, combined.slice.color)
      if (withlegend != FALSE) {
        ltext <- c(unique(lt), combined.slice.label)
        pie.colors.l <- c(unique(pie.colors.l), combined.slice.color)
      }
      if (ncol.legend == "auto") {
        if (withlegend == "bottom" || withlegend == "top") {
          ncol.legend <- ceiling(length(pie.colors) / 4)
        } else {
          ncol.legend <- 1
        }
      }
      pie.colors <- c(pie.colors, combined.slice.color)
      # Slices not combined
    } else {
      if (ncol.legend == "auto") {
        if (withlegend == "bottom" || withlegend == "top") {
          ncol.legend <- ceiling(ncol(x$emission_probs) / 4)
        } else {
          ncol.legend <- 1
        }
      }
    }

    if (!is.matrix(layout) && !is.function(layout) &&
         (layout == "horizontal" || layout == "vertical")) {
      if (length(dots) > 0) {
        plotcall <- as.call(c(list(plot.igraph2, g1, layout = glayout,
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
                         xlim = xlim, ylim = ylim, rescale=rescale, main = main), dots))
      } else {
        plotcall <- call("plot.igraph2", g1, layout = glayout,
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
                         xlim = xlim, ylim = ylim, rescale = rescale, main = main)
      }
    } else {
      if (length(dots) > 0) {
        plotcall <- as.call(c(list(plot.igraph2, g1, layout = glayout,
                         vertex.shape = "pie", vertex.pie = pie.values,
                         vertex.pie.color = list(pie.colors),
                         vertex.size = vertex.size,
                         vertex.label = vertex.label, vertex.label.dist = vertex.label.dist,
                         vertex.label.degree = vertex.label.pos,
                         vertex.label.family = vertex.label.family,
                         edge.curved = edge.curved, edge.width = edge.width,
                         edge.label = edge.label,
                         edge.label.family = edge.label.family,
                         edge.arrow.size = edge.arrow.size, main = main), dots))
      } else {
        plotcall <- call("plot.igraph2", g1, layout = glayout,
                         vertex.shape = "pie", vertex.pie = pie.values,
                         vertex.pie.color = list(pie.colors),
                         vertex.size = vertex.size,
                         vertex.label = vertex.label, vertex.label.dist = vertex.label.dist,
                         vertex.label.degree = vertex.label.pos,
                         vertex.label.family = vertex.label.family,
                         edge.curved = edge.curved, edge.width = edge.width,
                         edge.label = edge.label,
                         edge.label.family = edge.label.family,
                         edge.arrow.size = edge.arrow.size, main = main)
      }
    }
  } else {
    if(!is.matrix(layout) && !is.function(layout) &&
         (layout == "horizontal" || layout == "vertical")){
      if (length(dots) > 0) {
        plotcall <- as.call(c(list(plot.igraph2, g1, layout = glayout,
                         vertex.size = vertex.size,
                         vertex.label = vertex.label, vertex.label.dist = vertex.label.dist,
                         vertex.label.degree = vertex.label.pos,
                         vertex.label.family = vertex.label.family,
                         edge.curved = edge.curved, edge.width = edge.width,
                         edge.label = edge.label,
                         edge.label.family = edge.label.family,
                         xlim = xlim, ylim = ylim, rescale = rescale, main = main), dots))
      } else {
        plotcall <- call("plot.igraph2", g1, layout = glayout,
                         vertex.size = vertex.size,
                         vertex.label = vertex.label, vertex.label.dist = vertex.label.dist,
                         vertex.label.degree = vertex.label.pos,
                         vertex.label.family = vertex.label.family,
                         edge.curved = edge.curved, edge.width = edge.width,
                         edge.label = edge.label,
                         edge.label.family = edge.label.family,
                         xlim = xlim, ylim = ylim, rescale = rescale, main = main)
      }
    } else {
      if(length(dots) > 0){
        plotcall <- as.call(c(list(plot.igraph2, g1, layout = glayout,
                         vertex.size = vertex.size,
                         vertex.label = vertex.label, vertex.label.dist = vertex.label.dist,
                         vertex.label.degree = vertex.label.pos,
                         vertex.label.family = vertex.label.family,
                         edge.curved = edge.curved, edge.width = edge.width,
                         edge.label = edge.label,
                         edge.label.family = edge.label.family, main = main), dots))
      } else {
        plotcall <- call("plot.igraph2", g1, layout = glayout,
                         vertex.size = vertex.size,
                         vertex.label = vertex.label, vertex.label.dist = vertex.label.dist,
                         vertex.label.degree = vertex.label.pos,
                         vertex.label.family = vertex.label.family,
                         edge.curved = edge.curved, edge.width = edge.width,
                         edge.label = edge.label,
                         edge.label.family = edge.label.family, main = main)
      }
    }
  }


  # Plotting legend
  if (withlegend != FALSE && pie == TRUE) {
    legendcall <- call("seqlegend", seqdata = x$observations, cpal = pie.colors.l, ltext = ltext,
                       position = legend.pos, fontsize = cex.legend, ncol = ncol.legend,
                       with.missing = FALSE)

  } else {
    legendcall <- NULL
  }


  return(list(plotcall = plotcall, legendcall = legendcall))

  #graphics::layout(1)
}
