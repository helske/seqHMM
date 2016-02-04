#' Plot Multidimensional Sequence Plots in a Grid
#'
#' Function \code{gridplot} plots multiple \code{ssp} objects to a
#' grid.
#'
#'
#'
#' @export

#'
#' @param x A list of \code{\link{ssp}} objects.
#'
#' @param nrow,ncol Optional arguments to arrange plots.
#'
#' @param byrow Controls the order of plotting. Defaults to \code{FALSE}, i.e. plots
#'   are arranged column-wise.
#'
#' @param withlegend Defines if and how the legends for the states are plotted.
#'   The default value \code{"auto"} (equivalent to \code{TRUE} and
#'   \code{"many"}) creates separate legends for each requested plot. Other
#'   possibilities are \code{"combined"} (all legends combined) and \code{FALSE}
#'   (no legend).
#'
#' @param legend.pos Defines the positions of the legend boxes relative to the
#'   whole plot. Either one of \code{"bottom"} (equivalent to \code{"auto"}) or
#'   \code{"right"}, or a numerical vector of grid cells (by order) to print the
#'   legends to (the cells must be in one row/column).
#'
#' @param legend.pos2 Defines the positions of the legend boxes relative to the
#'   cell(s). One of \code{"bottomright"}, \code{"bottom"}, \code{"bottomleft"},
#'   \code{"left"}, \code{"topleft"}, \code{"top"} (the default), \code{"topright"},
#'   \code{"right"} and \code{"center"}.
#'
#' @param title.legend The titles for the legend boxes. The default \code{"auto"} takes
#'   the titles from the channel labels provided by the first object in \code{x}.
#'   \code{NA} prints no title.
#'
#' @param ncol.legend (A vector of) the number of columns for the legend(s). The
#'   default \code{"auto"} creates one column for each legend.
#'
#' @param with.missing.legend If set to \code{"auto"} (the default), a legend
#'   for the missing state is added automatically if one or more of the
#'   sequences in data contain missing states. With the value \code{TRUE} a
#'   legend for the missing state is added in any case; equivalently
#'   \code{FALSE} omits the legend for the missing state.
#'
#' @param cex.legend Expansion factor for setting the size of the font for the
#'   labels in the legend. The default value is 1. Values lesser than 1 will
#'   reduce the size of the font, values greater than 1 will increase the size.
#'
#' @param row.prop Sets the proportions of the row heights of the grid. The default
#'   value is \code{"auto"} for even row heights. Takes a vector of values from
#'   0 to 1, with values summing to 1.
#'
#' @param col.prop Sets the proportion of the column heights of the grid. The default
#'   value is \code{"auto"} for even column widths. Takes a vector of values
#'   from 0 to 1, with values summing to 1.
#'
#' @examples
#' \dontrun{
#' data("biofam3c")
#'
#' # Creating sequence objects
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
#' # Preparing plot for state distribution plots of observations for women
#' ssp_f <- ssp(
#'   list(child_seq[biofam3c$covariates$sex == "woman",],
#'        marr_seq[biofam3c$covariates$sex == "woman",],
#'        left_seq[biofam3c$covariates$sex == "woman",]),
#'   type = "d", plots = "obs", title = "Women",
#'   ylab = c("Children", "Married", "Left home"))
#'
#' # Preparing plot for state distribution plots of observations for men
#' # (Updating the previous plot, only arguments that change values)
#' ssp_m <- update(ssp_f, title = "Men",
#'   x = list(child_seq[biofam3c$covariates$sex == "man",],
#'        marr_seq[biofam3c$covariates$sex == "man",],
#'        left_seq[biofam3c$covariates$sex == "man",]))
#'
#' # Plotting state distribution plots of observations for women and men in two columns
#' gridplot(list(ssp_f, ssp_m), ncol = 2, withlegend = FALSE)
#'
#' # Preparing plots for women's state distributions
#' ssp_f2 <- ssp(
#'   list(marr_seq[biofam3c$covariates$sex == "woman",],
#'        child_seq[biofam3c$covariates$sex == "woman",],
#'        left_seq[biofam3c$covariates$sex == "woman",]),
#'   type = "d", border = NA, withlegend = FALSE,
#'   title = "State distributions for women", title.n = FALSE, xtlab = 15:30,
#'   ylab.pos = c(1, 2, 1), ylab = c("Married", "Children", "Left home"))
#'
#' # The same plot with sequences instead of state distributions
#' ssp_f3 <- update(
#'   ssp_f2, type = "I", sortv = "mds.obs", title = "Sequences for women")
#'
#' # State distributions with men's data
#' ssp_m2 <- update(
#'   ssp_f2, title = "State distributions for men",
#'   x = list(marr_seq[biofam3c$covariates$sex == "man",],
#'            child_seq[biofam3c$covariates$sex == "man",],
#'            left_seq[biofam3c$covariates$sex == "man",]))
#'
#' # Men's sequences
#' ssp_m3 <- update(
#'   ssp_m2, type = "I", sortv = "mds.obs", title = "Sequences for men")
#'
#' # Plotting state distributions and index plots of observations
#' # for women and men in two columns (+ one column for legends)
#' gridplot(
#'   list(ssp_f2, ssp_f3, ssp_m2, ssp_m3), ncol = 3, byrow = TRUE,
#'   withlegend = "combined", legend.pos = "right", col.prop = c(0.35, 0.35, 0.3))
#'
#' # The same with different positioning and fixed cells for legends
#' gridplot(
#'   list(ssp_f2, ssp_f3, ssp_m2, ssp_m3), ncol = 2, nrow = 3, byrow = TRUE,
#'   # defining the legend positions by the cell numbers
#'   legend.pos = 3:4)
#'  }
#'
#' @seealso \code{\link{ssp}} for defining the plot before using
#'   \code{gridplot}, and \code{\link{plot.ssp}} for plotting only one ssp object.

gridplot <- function(x, nrow = NA, ncol = NA, byrow = FALSE,
                     withlegend = "auto", legend.pos = "auto",
                     legend.pos2 = "center", title.legend = "auto",
                     ncol.legend = "auto",
                     with.missing.legend = "auto",
                     row.prop = "auto", col.prop = "auto", cex.legend = 1){
  grid.newpage()
  plot.new()
  opar <- par(no.readonly = TRUE)

  if (!is.na(nrow) && (nrow < 0 || nrow %% 1 != 0)) {
    stop("Argument nrow only accepts positive numeric values.")
  }
  if (!is.na(ncol) && (ncol < 0 || ncol %% 1 != 0)) {
    stop("Argument ncol only accepts positive numeric values.")
  }

  if (!is.logical(byrow)) {
    stop("Argument byrow only accepts values TRUE and FALSE.")
  }

  # Checking withlegend
  choices <- c(TRUE, FALSE, "auto", "combined", "many")
  ind <- pmatch(withlegend, choices)
  if (is.na(ind)) {
    stop("Argument withlegend must be one of TRUE, FALSE, \"auto\", \"many\", \"combined\"")
  }
  withlegend <- choices[ind]
  if (withlegend %in% c(TRUE, "auto")){
    withlegend <- "many"
  }

  if(!is.numeric(legend.pos)){
    choices <- c("bottom", "auto", "right")
    ind <- pmatch(legend.pos, choices)
    if (is.na(ind)) {
      stop("Argument legend.pos only accepts values \"bottom\", \"auto\", \"right\", or a numerical vector (one value for each legend).")
    }
    legend.pos <- choices[ind]
    if (legend.pos == "auto") {
      legend.pos <- "bottom"
    }
  }

  legend.pos2 <- match.arg(legend.pos2, c("bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right", "center"))

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
  } else if (is.numeric(col.prop) && all.equal(sum(col.prop),1) != TRUE) {
    stop("The elements of the vector provided for col.prop do not sum to 1.")
  }


  # Number of plots
  ngridplots <- length(x)

  # Checking for classes of x
  for (j in 1:ngridplots) {
    if (!inherits(x[[j]], "ssp")) {
      stop("At least one of your objects in x is not a ssp object. Use ssp to create one.")
    }
  }

  # Checks for automatic legends
  if (!is.na(withlegend) && withlegend != FALSE && ngridplots > 1) {
    nlegend <- x[[1]]$nplots
    for(i in 2:ngridplots){
      if(nlegend  !=  x[[i]]$nplots){
        warning("The number of legends is not the same in all requested plots. Legends could not be printed.")
        withlegend <- FALSE
        break()
      }
      nlegend <- x[[i]]$nplots
    }
  }

  if (x[[1]]$nchannels == 1 && withlegend == "many") {
    withlegend <- "combined"
  }

  # Set numbers of rows and columns if missing
  if (is.na(nrow) && is.na(ncol)) {
    nrow <- ceiling(sqrt(ngridplots))
    ncol <- ceiling(ngridplots / nrow)
    rcfixed <- "none"
  } else if(is.na(nrow)) {
    nrow <- ceiling(ngridplots / ncol)
    rcfixed <- "ncol"
  } else if(is.na(ncol)) {
    ncol <- ceiling(ngridplots / nrow)
    rcfixed <- "nrow"
  } else {
    rcfixed <- "both"
  }


  # New rows and cols for grid, possibly to modify later
  gridnrow <- nrow
  gridncol <- ncol
  emptycells <- gridnrow * gridncol - ngridplots

  # Enough cells without plots?
  if ((is.na(withlegend) && withlegend != FALSE) && length(legend.pos) > 1 && emptycells < ngridplots) {
    warning("There were not enough empty cells for the legends. Legends were positioned automatically.")
    legend.pos <- "bottom"
  }


  # Legend titles
  if (!is.na(withlegend) && withlegend != FALSE && !is.na(title.legend) &&
      title.legend != FALSE && !is.null(title.legend)) {
    # Wrong length for title.legend
    if (length(title.legend) > 1 || (length(title.legend) == 1 && title.legend != "auto")) {
      if (length(title.legend) != x[[1]]$nplots) {
        warning("The length of the vector provided for title.legend does not match the number of legends. Argument title.legend was set to \"auto\".")
        title.legend == "auto"
      }
    }
    # Automatic titles for legends (from the first ssp object)
    if (length(title.legend == 1 && title.legend == "auto")) {
      title.legend <- x[[1]]$ylab
    }
  }

  # Legend positions
  if (!is.na(withlegend) && withlegend != FALSE) {
    # Convert legend positions to "right"/"bottom" for layout
    if (length(legend.pos) > 1){
      if (byrow == TRUE){
        if (max(legend.pos) - min(legend.pos) + 1 > length(legend.pos)) {
          legendp <- "right"
        } else {
          legendp <- "bottom"
        }
        # byrow = FALSE
      } else {
        if (max(legend.pos) - min(legend.pos) + 1 > length(legend.pos)) {
          legendp <- "bottom"
        } else {
          legendp <- "right"
        }
      }
    }

    # Legends at bottom
    if (length(legend.pos) == 1 && legend.pos == "bottom") {
      legendp <- "bottom"
      if (!byrow) {
        if (emptycells == 0) {
          # Make room for legends
          if (rcfixed == "both") {
            warning(paste0("Legends do not fit into the requested grid with ", nrow, " nrows and ", ncol, " columns. A row was added."))
            gridnrow <- nrow + 1
          } else if (rcfixed == "ncol" || rcfixed == "none") {
            gridnrow <- nrow + 1
          } else {
            gridncol <- ncol + 1
          }
          emptycells <- gridnrow * gridncol - ngridplots
        }
        if (emptycells <= gridncol) {
          legend.pos <- c(c(gridncol:1) * gridnrow)[emptycells:1]
        } else {
          legend.pos <- c(gridncol:1) * gridnrow
        }
      # byrow = TRUE
      } else {
        if (emptycells == 0) {
          if (rcfixed == "both") {
            warning(paste0("Legends do not fit into the requested grid with ", nrow, " nrows and ", ncol, " columns. A row was added."))
            gridnrow <- nrow + 1
          } else if (rcfixed == "ncol" || rcfixed == "none") {
            gridnrow <- nrow + 1
          } else {
            gridncol <- ncol + 1
          }
          emptycells <- gridnrow * gridncol - ngridplots
        }
        if (emptycells <= gridncol) {
          legend.pos <- c((gridnrow * gridncol - emptycells+1):(gridnrow * gridncol))
        } else {
          legend.pos <- c((gridnrow * gridncol - floor(emptycells / gridncol) * gridncol + 1):(gridnrow * gridncol))
        }
      }
      if (length(ncol.legend) == 1 && ncol.legend == "auto" && withlegend != TRUE &&
         withlegend == "combined") {
        ncol.legend <- x[[1]]$nplots
      }
      # Legend at right
    } else if (length(legend.pos) == 1 && legend.pos == "right") {
      legendp <- "right"
      if (!byrow) {
        if (emptycells == 0) {
          if (rcfixed == "both") {
            warning(paste0("Legends do not fit into the requested grid with ", nrow, " nrows and ", ncol, " columns. A column was added."))
            gridncol <- ncol + 1
          } else if (rcfixed == "nrow" || rcfixed == "none") {
            gridncol <- ncol + 1
          } else {
            gridnrow <- nrow + 1
          }
          emptycells <- gridnrow * gridncol - ngridplots
        }
        if (emptycells <= gridnrow) {
          legend.pos <- c((gridnrow * gridncol - emptycells+1):(gridnrow * gridncol))
        } else {
          legend.pos <- c((gridnrow * gridncol - floor(emptycells / gridnrow) * gridnrow+1):(gridnrow * gridncol))
        }
        # byrow = TRUE
      } else {
        if (emptycells == 0) {
          if (rcfixed == "both") {
            warning(paste0("Legends do not fit into the requested grid with ", nrow, " nrows and ", ncol, " columns. A column was added."))
            gridncol <- ncol + 1
          } else if (rcfixed == "nrow" || rcfixed == "none") {
            gridncol <- ncol + 1
          } else {
            gridnrow <- nrow + 1
          }
          emptycells <- gridnrow * gridncol - ngridplots
        }
        if (emptycells <= gridnrow) {
          legend.pos <- c(c(gridnrow:1) * gridncol)[emptycells:1]
        } else {
          legend.pos <- c(gridnrow:1) * gridncol
        }
      }
      if (length(ncol.legend) == 1 && ncol.legend == "auto" && withlegend != TRUE &&
         withlegend == "combined") {
        ncol.legend <- 1
      }
    }

    # Legend rows and columns
    if (withlegend == "many") {
      if (length(ncol.legend) == 1 && ncol.legend == "auto") {
        ncol.legend <- rep(1, length.out = x[[1]]$nplots)
        legend.nrow <- sapply(lapply(x[[1]]$obs, "alphabet"), "length")
      } else if (length(ncol.legend) == 1 && x[[1]]$nplots > 1) {
        legend.nrow <- ceiling(sapply(lapply(x[[1]]$obs, "alphabet"), "length") / ncol.legend)
        ncol.legend <- rep(ncol.legend, length.out = x[[1]]$nplots)
      } else if (length(ncol.legend) < x[[1]]$nplots) {
        ncol.legend <- rep(ncol.legend, length.out = x[[1]]$nplots)
        legend.nrow <- ceiling(sapply(lapply(x[[1]]$obs, "alphabet"), "length") / ncol.legend)
      } else if (length(ncol.legend) > x[[1]]$nplots) {
        legend.nrow <- ceiling(sapply(lapply(x[[1]]$obs, "alphabet"), "length") / ncol.legend[1:x[[1]]$nplots])
      } else {
        legend.nrow <- ceiling(x$n_states / ncol.legend)
      }

      if(!is.na(title.legend) && title.legend != FALSE && !is.null(title.legend)){
        legend.nrow <- legend.nrow + 1
      }
    }
  }




  # Cells' proportions
  if (!is.numeric(row.prop) && row.prop == "auto") {
    row.prop <- rep(1 / gridnrow, gridnrow)
  }
  if (!is.numeric(col.prop) && col.prop == "auto") {
    col.prop <- rep(1 / gridncol, gridncol)
  }
  if (length(row.prop) != gridnrow) {
    warning("The length of the vector provided for row.prop does not match the number of nrow in the plot. Argument row.prop was changed to \"auto\".")
    row.prop <- rep(1 / gridnrow, gridnrow)
  }
  if (length(col.prop) != gridncol) {
    warning("The length of the vector provided for col.prop does not match the number of columns in the plot. Argument col.prop was changed to \"auto\".")
    col.prop <- rep(1 / gridncol, gridncol)
  }


  if (!byrow) {
    plotnrow <- rep(c(1:gridnrow), times = gridncol)
    plotncol <- rep(c(1:gridncol), each = gridnrow)
    plotgrid <- cbind(plotnrow, plotncol)
    if (!is.na(withlegend) && withlegend != FALSE) {
      lpos <- matrix(c(1:nrow(plotgrid)), byrow = FALSE, nrow = gridnrow)
      legendplace <- matrix(c(lpos[,] %in% legend.pos), byrow = FALSE, nrow = gridnrow)
      plotplace <- matrix(c(!(lpos[,] %in% legend.pos)), byrow = FALSE, nrow = gridnrow)
      plotcells <- plotgrid[c(plotplace), , drop = FALSE][1:ngridplots, ]
      plotlegend <- plotgrid[c(legendplace), , drop = FALSE]
      if (legendp == "bottom") {
        plotlegend <- plotlegend[order(plotlegend[,2]), , drop = FALSE]
      }
      if (length(ncol.legend) == 1 && ncol.legend == "auto") {
        if (withlegend == "combined") {
          if (length(unique(plotlegend[, 1])) >= length(unique(plotlegend[, 2]))) {
            ncol.legend <- 1
          } else {
            ncol.legend <- x[[1]]$nplots
          }
        } else {
          rep(1, x[[1]]$nplots)
        }
      }
    } else {
      plotcells <- plotgrid[1:ngridplots, , drop = FALSE]
    }
  # byrow = TRUE
  } else {
    plotnrow <- rep(c(1:gridnrow), each = gridncol)
    plotncol <- rep(c(1:gridncol), times = gridnrow)
    plotgrid <- cbind(plotnrow, plotncol)
    if (!is.na(withlegend) && withlegend != FALSE) {
      lpos <- matrix(c(1:nrow(plotgrid)), byrow = TRUE, nrow = gridnrow)
      legendplace <- matrix(c(lpos[,] %in% legend.pos), byrow = FALSE, nrow = gridnrow)
      plotplace <- matrix(c(!(lpos[,] %in% legend.pos)), byrow = FALSE, nrow = gridnrow)
      plotcells <- plotgrid[c(t(plotplace)), , drop = FALSE][1:ngridplots, ]
      plotlegend <- plotgrid[c(t(legendplace)), , drop = FALSE]
      if (legendp == "right") {
        plotlegend <- plotlegend[order(plotlegend[, 2]), , drop = FALSE]
      }
      if (length(ncol.legend) == 1 && ncol.legend == "auto") {
        if (withlegend == "combined") {
          if (length(unique(plotlegend[, 1])) >= length(unique(plotlegend[, 2]))) {
            ncol.legend <- 1
          } else {
            ncol.legend <- x[[1]]$nplots
          }
        } else {
          rep(1, x[[1]]$nplots)
        }
      }
    } else {
      plotcells <- plotgrid[1:ngridplots, , drop = FALSE]
    }
  }


  multitop.vp <- viewport(layout =
                            grid.layout(gridnrow,gridncol,
                                        widths = do.call(unit,
                                                         args = list(col.prop,
                                                                   rep("npc",
                                                                       length(col.prop)))),
                                        heights = do.call(unit,
                                                          args = list(row.prop,
                                                                    rep("npc",
                                                                        length(row.prop))))),
                          width = unit(1, "npc"))
  for (i in 1:ngridplots) {
    assign(paste0("vpplot",i), viewport(layout.pos.row = plotcells[i, 1],
                                        layout.pos.col = plotcells[i, 2],
                                        name = paste0("vpplot", i)))
  }
  if (!is.na(withlegend) && withlegend != FALSE) {
    assign("vplegend", viewport(layout.pos.row = unique(plotlegend[, 1]),
                                layout.pos.col = unique(plotlegend[, 2]),
                                name = "vplegend"))
    vpall <- vpTree(multitop.vp, do.call(vpList,
                                         args = mget(c(paste0("vpplot", 1:ngridplots),
                                                     "vplegend"))))
  } else {
    vpall <- vpTree(multitop.vp, do.call(vpList, args = mget(paste0("vpplot", 1:ngridplots))))
  }


  pushViewport(vpall)

  upViewport()

  # Plots
  for (p in 1:ngridplots) {
    downViewport(paste0("vpplot", p))
    do.call(SSPlotter, args = x[[p]])
    popViewport()
    upViewport()
  }

  for (i in 1:ngridplots) {
    if (x[[i]]$nchannels == 1) {
      x[[i]]$obs <- list(x[[i]]$obs)
    }
  }

  # Legends
  if (!is.na(withlegend) && withlegend != FALSE) {
    # Hidden states (all in one legend)
    if (x[[1]]$plots == "both" || x[[1]]$plots == "hidden.paths") {
      hstext <- NULL
      hscpal <- NULL
      if (x[[1]]$nchannels > 1) {
        for(i in 1:x[[1]]$nchannels){
          hstext <- c(hstext, attr(x[[1]]$hidden.paths.seq, "labels"))
          hscpal <- c(hscpal, attr(x[[1]]$hidden.paths.seq, "cpal"))
        }
      } else {
        hstext <- c(hstext, attr(x[[1]]$hidden.paths.seq, "labels"))
        hscpal <- c(hscpal, attr(x[[1]]$hidden.paths.seq, "cpal"))
      }
    }

    ltext <- NULL
    cpal <- NULL
    if (x[[1]]$plots == "both" || x[[1]]$plots == "obs") {
      ltexts <- rep(list(NULL), x[[1]]$nchannels)
      cpals <- rep(list(NULL), x[[1]]$nchannels)
      anymissing <- FALSE
      for (i in 1:ngridplots) {
        for (j in 1:x[[1]]$nchannels) {
          ltexts[[j]] <- unique(c(c(ltexts[[j]], attr(x[[i]]$obs[[j]], "labels"))))
          cpals[[j]] <- unique(c(cpals[[j]], attr(x[[i]]$obs[[j]], "cpal")))
          if (with.missing.legend == "auto") {
            if (any(x[[i]]$obs[[j]] == "*")) {
              anymissing <- TRUE
            }
          } else {
            anymissing <- with.missing.legend
          }
        }
      }
      for (j in 1:x[[1]]$nchannels) {
        ltext <- c(ltext, unlist(ltexts[[j]]))
        cpal <- c(cpal, unlist(cpals[[j]]))
      }
    }
#     if(x[[1]]$plots == "both" || x[[1]]$plots == "hidden.paths"){
#       ltext <- c(ltext, attr(x[[maxhsplot]]$hidden.paths.seq, "labels"))
#       cpal <- c(cpal, attr(x[[maxhsplot]]$hidden.paths.seq, "cpal"))
#     }

    # Separate legends
    if (withlegend == "many") {
      downViewport("vplegend")
      # Vertical legends
      if(legendp == "right"){
        pushViewport(viewport(layout =
                                grid.layout(nrow = x[[1]]$nplots, ncol = 1,
                                            heights = unit((legend.nrow / sum(legend.nrow)),
                                                         "npc")),
                              width = unit(0.95, "npc")))
        # Legends for channels
        if (x[[1]]$plots == "both" || x[[1]]$plots == "obs") {
          for(i in 1:x[[1]]$nchannels){
            pushViewport(viewport(layout.pos.col = 1, layout.pos.row = i))
            par(plt = gridPLT(), new = TRUE)
            seqlegend(x[[1]]$obs[[i]], fontsize = cex.legend, position = legend.pos2,
                      cpal = cpals[[i]], ltext = ltexts[[i]],
                      ncol = ncol.legend[i], with.missing = with.missing.legend,
                      title = title.legend[i])
            popViewport()
          }
        }
        # Legends for hidden paths
        if (x[[1]]$plots == "both" || x[[1]]$plots == "hidden.paths") {
          pushViewport(viewport(layout.pos.col = 1, layout.pos.row = x[[1]]$nplots))
          par(plt = gridPLT(), new = TRUE)
          seqlegend(x[[1]]$hidden.paths.seq, fontsize = cex.legend,
                    position = legend.pos2, ncol = ncol.legend[length(ncol.legend)],
                    cpal = hscpal, ltext = hstext,
                    with.missing = with.missing.legend,
                    title = title.legend[length(title.legend)])
          popViewport()
        }
        popViewport()
        # Horizontal legends
      } else {
        pushViewport(viewport(layout =
                                grid.layout(ncol = x[[1]]$nplots, nrow = 1,
                                            widths = unit((legend.nrow / sum(legend.nrow)),
                                                        "npc")),
                              width = unit(0.95, "npc")))
        # Legends for channels
        if (x[[1]]$plots == "both" || x[[1]]$plots == "obs") {
          for(i in 1:x[[1]]$nchannels){
            pushViewport(viewport(layout.pos.col = i, layout.pos.row = 1))
            par(plt = gridPLT(), new = TRUE)
            seqlegend(x[[1]]$obs[[i]], fontsize = cex.legend, position = legend.pos2,
                      cpal = cpals[[i]], ltext = ltexts[[i]],
                      ncol = ncol.legend[i], with.missing = with.missing.legend,
                      title = title.legend[i])
            popViewport()
          }
        }
        # Legends for most probable paths
        if(x[[1]]$plots == "both" || x[[1]]$plots == "hidden.paths"){
          pushViewport(viewport(layout.pos.col = x[[1]]$nplots, layout.pos.row = 1))
          par(plt = gridPLT(), new = TRUE)
          seqlegend(x[[1]]$hidden.paths.seq, fontsize = cex.legend,
                    position = legend.pos2, ncol = ncol.legend[length(ncol.legend)],
                    cpal = hscpal, ltext = hstext,
                    with.missing = with.missing.legend,
                    title = title.legend[length(title.legend)])
          popViewport()
        }
        popViewport()
      }
      # Combined legends
    } else if (withlegend == "combined") {
      if (x[[1]]$plots == "both" || x[[1]]$plots == "obs") {
        downViewport("vplegend")
        par(plt = gridPLT(), new = TRUE)
        pushViewport(viewport(width = unit(0.9, "npc")))

        seqlegend(x[[1]]$obs[[1]], fontsize = cex.legend, position = legend.pos2,
                  ncol = ncol.legend, cpal = cpal, ltext = ltext,
                  with.missing = anymissing,
                  missing.color = attr(x[[1]]$obs[[1]],"missing.color"))

        popViewport()

      }
      # Legends for most probable paths
      if(x[[1]]$plots == "both" || x[[1]]$plots == "hidden.paths"){
        ltext <- c(ltext, hstext)
        cpal <- c(cpal, hscpal)

        downViewport("vplegend")
        par(plt = gridPLT(), new = TRUE)
        pushViewport(viewport(width = unit(0.9, "npc")))

        seqlegend(x[[1]]$obs[[1]], fontsize = cex.legend, position = legend.pos2,
                  ncol = ncol.legend, cpal = cpal, ltext = ltext,
                  with.missing = anymissing,
                  missing.color = attr(x[[1]]$obs[[1]],"missing.color"))

        popViewport()
      }


    }
  }

  par(opar)
}

