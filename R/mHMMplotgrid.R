mHMMplotgrid <- function(x, which.plots = NULL, nrow = NA, ncol = NA, byrow = FALSE,
                         row.prop = "auto", col.prop = "auto", layout = "horizontal", pie = TRUE,
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
                         withlegend = "bottom", legend.pos = "center", ltext = NULL, legend.prop = 0.5,
                         cex.legend = 1, ncol.legend = "auto", cpal = "auto",
                         main = "auto", ...) {

  plot.new()
  opar <- par(no.readonly = TRUE)
  on.exit(opar, add = TRUE)
  on.exit(graphics::layout(1), add = TRUE)

  if (is.null(main)) {
    par(mar = c(0.5, 0.5, 0.5, 0.5))
  } else {
    par(mai = c(0, 0, 1, 0))
  }

  if (!is.null(main)) {
    if (length(main) == 1 && main == "auto") {
      main <- x$cluster_names
    } else if (length(main) != length(x$cluster_names)) {
      warning("The length of the vector provided for the main argument does not match the length of x$cluster_names. Using cluster_names instead.")
      main <- x$cluster_names
    }
  }

  divmodels <- separate_mhmm(x)


  if (!is.na(nrow) && (nrow < 0 || nrow %% 1 != 0)) {
    stop("Argument nrow only accepts positive numeric values.")
  }
  if (!is.na(ncol) && (ncol < 0 || ncol %% 1 != 0)) {
    stop("Argument ncol only accepts positive numeric values.")
  }

  if (!is.logical(byrow)) {
    stop("Argument byrow only accepts values TRUE and FALSE.")
  }

  if (!is.numeric(ncol.legend) && ncol.legend != "auto") {
    stop("Argument ncol.legend only accepts values \"auto\" or a numerical vector.")
  }

  if (!is.numeric(row.prop) && row.prop != "auto") {
    stop("Argument row.prop only accepts values \"auto\" or a numerical vector.")
  } else if (is.numeric(row.prop) && all.equal(sum(row.prop), 1) != TRUE) {
    stop("The elements of the vector provided for row.prop do not sum to 1.")
  }

  if (!is.numeric(col.prop) && col.prop != "auto") {
    stop("Argument col.prop only accepts values \"auto\" or a numerical vector.")
  } else if (is.numeric(col.prop) && all.equal(sum(col.prop), 1) != TRUE) {
    stop("The elements of the vector provided for col.prop do not sum to 1.")
  }


  # Number of plots
  if (is.null(which.plots)) {
    which.plots <- 1:x$n_clusters
  }
  ngridplots <- length(which.plots)


  if (is.na(nrow) && is.na(ncol)) {
    nrow <- ceiling(sqrt(ngridplots))
    ncol <- ceiling(ngridplots/nrow)
  } else if (is.na(nrow)) {
    nrow <- ceiling(ngridplots/ncol)
  } else if (is.na(ncol)) {
    ncol <- ceiling(ngridplots/nrow)
  }

#   # Number of columns in legends
#   if (!is.na(withlegend) && withlegend != FALSE) {
#     if (length(ncol.legend) == 1 && ncol.legend == "auto") {
#       ncol.legend <- rep(2, ngridplots)
#     } else if (length(ncol.legend) == 1 && x$n_clusters > 1) {
#       ncol.legend <- rep(ncol.legend, ngridplots)
#     } else if (length(ncol.legend) != ngridplots) {
#       vertex.label  <- rep(ncol.legend, length.out = ngridplots)
#     }
#   }





  # Cells' proportions
  if (!is.numeric(row.prop) && row.prop == "auto") {
    row.prop <- rep(1/nrow, nrow)
  }
  if (!is.numeric(col.prop) && col.prop == "auto") {
    col.prop <- rep(1/ncol, ncol)
  }
  if (length(row.prop) != nrow) {
    stop("The length of the vector provided for row.prop does not match the number of nrow in the plot.")
  }
  if (length(col.prop) != ncol) {
    stop("The length of the vector provided for col.prop does not match the number of columns in the plot.")
  }

  # Plotting order for layout
  if (!is.na(withlegend) && withlegend != FALSE) {
    if (!byrow) {
      plotlayout <- matrix(c(1:ngridplots,
                             rep(0, nrow * ncol - ngridplots)), nrow = nrow)
      legendlayout <- matrix(c((ngridplots + 1):(2 * ngridplots),
                               rep(0, nrow * ncol - ngridplots)), nrow = nrow)
      if (withlegend == "right") {
        # Matrix for layout
        lmatrix <- cbind(plotlayout[, 1], legendlayout[, 1])
        if (ncol > 1) {
          for (i in 2:ncol) {
            lmatrix <- cbind(lmatrix, plotlayout[, i], legendlayout[, i])
          }
        }
        cprops <- c(col.prop[1] * (1-legend.prop), col.prop[1] * legend.prop)
        if (ncol > 1) {
          for (i in 2:ncol) {
            cprops <- c(cprops, col.prop[i] * (1 - legend.prop),
                        col.prop[i] * legend.prop)
          }
        }
        rprops <- row.prop
      } else if (withlegend == "left") {
        lmatrix <- cbind(legendlayout[, 1], plotlayout[, 1])
        if (ncol > 1) {
          for (i in 2:ncol) {
            lmatrix <- cbind(lmatrix, legendlayout[, i], plotlayout[, i])
          }
        }
        cprops <- c(col.prop[1] * legend.prop, col.prop[1] * (1 - legend.prop))
        if (ncol > 1) {
          for (i in 2:ncol) {
            cprops <- c(cprops, col.prop[i] * legend.prop,
                        col.prop[i] * (1 - legend.prop))
          }
        }
        rprops <- row.prop
      } else if (withlegend == "bottom") {
        lmatrix <- rbind(plotlayout[1, ], legendlayout[1, ])
        if (nrow > 1) {
          for (i in 2:nrow) {
            lmatrix <- rbind(lmatrix, plotlayout[i, ], legendlayout[i, ])
          }
        }
        rprops <- c(row.prop[1] * (1 - legend.prop), row.prop[1] * legend.prop)
        if (nrow > 1) {
          for (i in 2:nrow) {
            rprops <- c(rprops,row.prop[i] * (1 - legend.prop),
                        row.prop[i] * legend.prop)
          }
        }
        cprops <- col.prop
        # withlegend == "top"
      } else {
        lmatrix <- rbind(legendlayout[1, ], plotlayout[1, ])
        if (nrow > 1) {
          for (i in 2:nrow) {
            lmatrix <- rbind(lmatrix, legendlayout[i, ], plotlayout[i, ])
          }
        }
        rprops <- c(row.prop[1] * legend.prop,row.prop[1] * (1 - legend.prop))
        if (nrow > 1) {
          for (i in 2:nrow) {
            rprops <- c(rprops, row.prop[i] * legend.prop,
                        row.prop[i] * (1 - legend.prop))
          }
        }
        cprops <- col.prop
      }
      # byrow = TRUE
    } else {
      plotlayout <- matrix(c(1:ngridplots, rep(0, nrow * ncol - ngridplots)),
                           nrow = nrow, byrow = TRUE)
      legendlayout <- matrix(c((ngridplots + 1):(2 * ngridplots),
                               rep(0,nrow * ncol - ngridplots)),
                             nrow = nrow, byrow = TRUE)
      if (nrow * ncol > ngridplots) {
        plotlayout[plotlayout > ngridplots] <- 0
        legendlayout[legendlayout > (2 * ngridplots)] <- 0
      }
      if (withlegend == "right") {
        # Matrix for layout
        lmatrix <- cbind(plotlayout[, 1], legendlayout[, 1])
        if (ncol > 1) {
          for (i in 2:ncol) {
            lmatrix <- cbind(lmatrix, plotlayout[, i], legendlayout[, i])
          }
        }
        cprops <- c(col.prop[1] * (1 - legend.prop), col.prop[1] * legend.prop)
        if (ncol > 1) {
          for (i in 2:ncol) {
            cprops <- c(cprops,col.prop[i] * (1 - legend.prop),
                        col.prop[i] * legend.prop)
          }
        }
        rprops <- row.prop
      } else if (withlegend == "left") {
        lmatrix <- cbind(legendlayout[,1], plotlayout[,1])
        if (ncol > 1) {
          for (i in 2:ncol) {
            lmatrix <- cbind(lmatrix,legendlayout[,i], plotlayout[,i])
          }
        }
        cprops <- c(col.prop[1] * legend.prop,col.prop[1] * (1-legend.prop))
        if (ncol > 1) {
          for (i in 2:ncol) {
            cprops <- c(cprops,col.prop[i] * legend.prop,col.prop[i] *
                          (1-legend.prop))
          }
        }
        rprops <- row.prop
      } else if (withlegend == "bottom") {
        lmatrix <- rbind(plotlayout[1, ], legendlayout[1, ])
        if (nrow > 1) {
          for (i in 2:nrow) {
            lmatrix <- rbind(lmatrix, plotlayout[i, ], legendlayout[i, ])
          }
        }
        rprops <- c(row.prop[1] * (1-legend.prop), row.prop[1] * legend.prop)
        if (nrow > 1) {
          for (i in 2:nrow) {
            rprops <- c(rprops,row.prop[i] * (1 - legend.prop),
                        row.prop[i] * legend.prop)
          }
        }
        cprops <- col.prop
        # "top"
      } else {
        lmatrix <- rbind(legendlayout[1, ], plotlayout[1, ])
        if (nrow > 1) {
          for (i in 2:nrow) {
            lmatrix <- rbind(lmatrix, legendlayout[i, ], plotlayout[i, ])
          }
        }
        rprops <- c(row.prop[1] * legend.prop, row.prop[1] * (1 - legend.prop))
        if (nrow > 1) {
          for (i in 2:nrow) {
            rprops <- c(rprops,row.prop[i] * legend.prop,row.prop[i] *
                          (1-legend.prop))
          }
        }
        cprops <- col.prop
      }
    }
    # No legends
  } else {
    if (!byrow) {
      lmatrix <- matrix(c(1:ngridplots, rep(0, nrow * ncol-ngridplots)),
                           nrow = nrow)
      cprops <- col.prop
      rprops <- row.prop
      # byrow = TRUE
    } else {
      lmatrix <- matrix(c(1:ngridplots, rep(0, nrow * ncol - ngridplots)),
                           nrow = nrow, byrow = TRUE)
      cprops <- col.prop
      rprops <- row.prop
    }
  }

  graphics::layout(lmatrix, widths = cprops, heights = rprops)


  # Plotting arguments for graphs and legends
  HMMcalls <- list()
  length(HMMcalls) <- ngridplots
  args <- as.list(match.call())[-1]

  args$which.plots <- args$nrow <- args$ncol <- args$byrow <-
    args$row.prop <- args$col.prop <- NULL

  for (p in which.plots) {
    if (length(ncol.legend) > 1) {
      ncolleg <- ncol.legend[p]
    } else {
      ncolleg <- ncol.legend
    }
    args$x <- divmodels[[p]]
    args$main <- main[p]
    args$ncol.legend <- ncolleg
    HMMcalls[[p]] <- do.call(HMMplot, args = args)
  }


  # Plotting graphs
  for (p in which.plots) {
    eval(HMMcalls[[p]]$plotcall)
  }

  # Plotting legends
  if (withlegend != FALSE) {
    for (p in which.plots) {
      eval(HMMcalls[[p]]$legendcall)
    }
  }

  graphics::layout(1)
}
