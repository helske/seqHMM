#' Interactive Stacked Plots of Multichannel Sequences and/or Most Probable
#' Paths for Mixture Hidden Markov Models
#'
#' Function `mssplot` plots stacked sequence plots of observation sequences
#' and/or most probable hidden state paths for each model of the `mhmm`
#' object (model chosen according to the most probable path).
#'
#' @export
#'
#' @param x Mixture hidden Markov model object of class `mhmm`.
#'
#' @param ask If `TRUE` and `which.plots` is NULL, `plot.mhmm` operates in interactive mode, via [menu()]. Defaults to `FALSE`.
#'
#' @param which.plots The number(s) of the requested model(s) as an integer vector. The default `NULL` produces all plots.
#'
#' @param hidden.paths Output from the [hidden_paths()] function. The
#'   default value `NULL` computes hidden paths automatically, if needed.
#'
#' @param plots What to plot. One of `"obs"` for observations (the default),
#'   `"hidden.paths"` for most probable paths of hidden states,
#'   or `"both"` for observations and hidden paths together.
#'
#' @param type The type of the plot. Available types are `"I"` for index
#'   plots and `"d"` for state distribution plots (the default). See
#'   [TraMineR::seqplot()] for details.
#'
#' @param tlim Indexes of the subjects to be plotted (the default is 0,
#' i.e. all subjects are plotted). For example, `tlim = 1:10` plots
#' the first ten subjects in data.
#'
#' @param sortv A sorting variable or a sort method (one of `"from.start"`,
#'   `"from.end"`, `"mds.obs"`, or `"mds.hidden"`) for
#'   `type = "I"`. The value `"mds.hidden"` is only available when
#'   `which = "both"` and `which = "hidden.paths"`. Options `"mds.obs"` and
#'   `"mds.hidden"` automatically arrange the sequences according to the
#'   scores of multidimensional scaling (using [stats::cmdscale()]) for the
#'   observed data or hidden states paths.
#'   MDS scores are computed from distances/dissimilarities using a metric
#'   defined in argument `dist.method`. See [TraMineR::plot.stslist()] for
#'   more details on `"from.start"` and `"from.end"`.
#'
#' @param sort.channel The number of the channel according to which the
#'   `"from.start"` or `"from.end"` sorting is done. Sorting according
#'   to hidden states is called with value 0. The default value is 1 (the first
#'   channel).
#'
#' @param dist.method The metric to be used for computing the distances of the
#'   sequences if multidimensional scaling is used for sorting. One of "OM"
#'   (optimal matching, the default), "LCP" (longest common prefix), "RLCP"
#'   (reversed LCP, i.e. longest common suffix), "LCS" (longest common
#'   subsequence), "HAM" (Hamming distance), and "DHD" (dynamic Hamming distance).
#'   Transition rates are used for defining substitution costs if needed. See
#'   [TraMineR::seqdef()] for more information on the metrics.
#'
#' @param with.missing Controls whether missing states are included in state
#'   distribution plots (`type = "d"`). The default is `FALSE`.
#'
#' @param missing.color Alternative color for representing missing values
#'   in the sequences. By default, this color is taken from the `missing.color`
#'   attribute of the sequence object.
#'
#' @param title A vector of main titles for the graphics. The default is `NA`: if
#'   `title.n = TRUE`, the name of the cluster and the number of subjects is plotted.
#'   `FALSE` prints no titles, even when `title.n = TRUE`.
#'
#' @param title.n Controls whether the number of subjects is printed in the main
#'   titles of the plots. The default is `TRUE`: n is plotted if `title`
#'   is anything but `FALSE`.
#'
#' @param cex.title Expansion factor for setting the size of the font for the main
#'   titles. The default value is 1. Values lesser than 1 will reduce the size of
#'   the font, values greater than 1 will increase the size.
#'
#' @param title.pos Controls the position of the main titles of the plots. The
#'   default value is 1. Values greater than 1 will place the title higher.
#'
#' @param with.legend Defines if and where the legend for the states is plotted.
#'   The default value `"auto"` (equivalent to `TRUE` and
#'   `"right"`) creates separate legends for each requested plot and
#'   positions them on the right-hand side of the plot. Other possible values
#'   are `"bottom"`,
#'   `"right.combined"`, and `"bottom.combined"`, of which the last
#'   two create a combined legend in the selected position. `FALSE` prints no legend.
#'
#' @param ncol.legend (A vector of) the number of columns for the legend(s). The
#'   default `"auto"` creates one column for each legend.
#'
#' @param with.missing.legend If set to `"auto"` (the default), a legend
#'   for the missing state is added automatically if one or more of the
#'   sequences in the data/channel contains missing states and `type = "I"`.
#'   If `type = "d"` missing states are omitted from the legends unless
#'   `with.missing = TRUE`. With the value `TRUE` a
#'   legend for the missing state is added in any case; equivalently
#'   `FALSE` omits the legend for the missing state.
#'
#' @param legend.prop Sets the proportion of the graphic area used for plotting
#'   the legend when `with.legend` is not `FALSE`. The default value is
#'   0.3. Takes values from 0 to 1.
#'
#' @param cex.legend Expansion factor for setting the size of the font for the
#'   labels in the legend. The default value is 1. Values lesser than 1 will
#'   reduce the size of the font, values greater than 1 will increase the size.
#'
#' @param hidden.states.colors A vector of colors assigned to hidden states. The default
#'   value `"auto"` uses the colors assigned to the `stslist` object (created
#'   with [TraMineR::seqdef()]) if `hidden.paths` is given; otherwise colors from
#'   [colorpalette()] are automatically used.
#'
#' @param hidden.states.labels Labels for the hidden states. The default value
#'   `"auto"` uses the names provided in `x$state_names` if `x` is
#'   an `hmm` object; otherwise the number of the hidden state.
#'
#' @param xaxis Controls whether an x-axis is plotted below the plot at the
#'   bottom. The default value is `TRUE`.
#'
#' @param xlab An optional label for the x-axis. If set to `NA`, no label
#'   is drawn.
#'
#' @param xtlab Optional labels for the x-axis tick labels.  If unspecified, the
#'   column names of the `seqdata` sequence object are used (see
#'   [TraMineR::seqdef()]).
#'
#' @param xlab.pos Controls the position of the x-axis label. The default value
#'   is 1. Values greater than 1 will place the label further away from the plot.
#'
#' @param yaxis Controls whether or not to plot the y-axis. The default is `FALSE`.
#'
#' @param ylab Labels for the channels shown as labels for y-axes.
#'   A vector of names for each channel
#'   (observations). The default value `"auto"` uses the names provided in
#'   `x$channel_names` if `x` is an `hmm` object; otherwise the
#'   names of the list in `x` if given, or the
#'   number of the channel if names are not given. `FALSE` prints no labels.
#'
#' @param hidden.states.title Optional label for the hidden state plot (in the
#'   y-axis). The default is `"Hidden states"`.
#'
#' @param ylab.pos Controls the position of the y axis labels (labels for
#'   channels and/or hidden states). Either `"auto"` or a numerical vector
#'   indicating how far away from the plots the titles are positioned. The
#'   default value `"auto"` positions all titles on line 1.
#'   Shorter vectors are recycled.
#'
#' @param cex.lab Expansion factor for setting the size of the font for the axis
#'   labels. The default value is 1. Values lesser than 1 will reduce the size
#'   of the font, values greater than 1 will increase the size.
#'
#' @param cex.axis Expansion factor for setting the size of the font for the x-axis
#'   tick labels. The default value is 1. Values lesser than 1 will reduce the size of
#'   the font, values greater than 1 will increase the size.
#' @param respect_void  If `TRUE` (default), states at the time points
#' corresponding to TraMineR's void in the observed sequences are set to void
#' in the hidden state sequences as well.
#' @param ... Other arguments to be passed on to
#'   [TraMineR::seqplot()].
#'
#' @seealso [build_mhmm()] and [fit_model()] for building and
#'   fitting mixture hidden Markov models, [hidden_paths()] for
#'   computing the most probable paths (Viterbi paths) of hidden states,
#'   [plot.mhmm()] for plotting `mhmm` objects as directed graphs, and
#'   [colorpalette()] for default colors.
#'
mssplot <- function(
    x, ask = FALSE, which.plots = NULL, hidden.paths = NULL,
    plots = "obs", type = "d", tlim = 0,
    sortv = NULL, sort.channel = 1, dist.method = "OM",
    with.missing = FALSE, missing.color = NULL,
    title = NA, title.n = TRUE, cex.title = 1, title.pos = 1,
    with.legend = "auto", ncol.legend = "auto",
    with.missing.legend = "auto",
    legend.prop = 0.3, cex.legend = 1,
    hidden.states.colors = "auto", hidden.states.labels = "auto",
    xaxis = TRUE, xlab = NA, xtlab = NULL, xlab.pos = 1,
    ylab = "auto", hidden.states.title = "Hidden states",
    yaxis = FALSE, ylab.pos = "auto",
    cex.lab = 1, cex.axis = 1, respect_void = TRUE, ...) {
  .Deprecated("stacked_sequence_plot")

  # Checking for class of x
  if (!inherits(x, "mhmm")) {
    stop("Your object x is not a mhmm object. Use build_mhmm to create one.")
  }

  oldPar <- par(no.readonly = TRUE)
  on.exit(par(oldPar), add = TRUE)

  oldWarn <- options("warn")
  options(warn = 1)
  on.exit(options(oldWarn), add = TRUE)

  # ssp arguments (besides mhmm object and hidden.paths)
  args <- as.list(match.call())[-(1:2)]
  if ("ask" %in% names(args)) {
    args <- args[-which(names(args) == "ask")]
  }
  if ("which.plots" %in% names(args)) {
    args <- args[-which(names(args) == "which.plots")]
  }
  if ("hidden.paths" %in% names(args)) {
    args <- args[-which(names(args) == "hidden.paths")]
  }
  if (!("title" %in% names(args))) {
    titles <- x$cluster_names
  } else {
    if (length(title) == 1 && (is.na(title) || !title)) {
      titles <- rep(title, x$n_clusters)
    } else if (length(title) != x$n_clusters) {
      warning("The length of the vector provided for the title argument does not match the number of clusters. Automatic titles were used instead.")
      titles <- x$cluster_names
    } else {
      titles <- eval(args$title)
    }
    args <- args[-which(names(args) == "title")]
  }
  if (length(ylab) == 1 && ylab == "auto") {
    args$ylab <- x$channel_names
  }

  if (is.null(hidden.paths)) {
    hidden.paths <- suppressWarnings(suppressMessages(hidden_paths(x, respect_void = respect_void)))
  }

  if (!("hidden.states.labels" %in% names(args))) {
    hidden.states.labels <- unlist(x$state_names)
  }
  hidden.pathslabs <- list()
  k <- 0
  for (i in 1:x$n_clusters) {
    hidden.pathslabs[[i]] <- hidden.states.labels[(k + 1):(k + x$n_states[i])]
    k <- k + x$n_states[i]
  }
  n_alphabet <- length(alphabet(hidden.paths))
  if (!("hidden.states.colors" %in% names(args))) {
    if (n_alphabet <= 200) {
      hidden.states.colors <- seqHMM::colorpalette[[n_alphabet]]
    } else {
      stop(
        "Model contains ", n_alphabet, " hidden states, which ",
        " is more than supported by the default color palette. Specify your ",
        " own color palette with the argument 'hidden.states.colors'."
      )
    }
  }
  if (n_alphabet != length(hidden.states.colors)) {
    stop(
      "The number of hidden states is ", n_alphabet,
      " but the supplied color palette contains only ",
      length(hidden.states.colors), "colours."
    )
  }

  hidden.pathscols <- list()
  k <- 0
  for (i in 1:x$n_clusters) {
    hidden.pathscols[[i]] <- hidden.states.colors[(k + 1):(k + x$n_states[i])]
    k <- k + unname(x$n_states[i])
  }


  # Clusters determined by hidden_paths
  hp_by_cluster_logic <- rep(list(NULL), x$n_clusters)
  hp_by_cluster <- rep(list(NULL), x$n_clusters)
  mm <- NULL
  for (i in 1:x$n_clusters) {
    # Find matching cluster names from the first hidden state of each individual
    if (length(unique(unlist(x$state_names))) == length(unlist(x$state_names))) {
      hp_by_cluster_logic[[i]] <- hidden.paths[, 1] %in% x$state_names[[i]]
    } else {
      hp_by_cluster_logic[[i]] <- grepl(paste0(x$cluster_names[i], ":"), hidden.paths[, 1])
    }
    hp_by_cluster[[i]] <- hidden.paths[hp_by_cluster_logic[[i]], ]
    # Give a warning, if no subjects assigned to cluster
    if (sum(hp_by_cluster_logic[[i]]) == 0) {
      mm <- c(mm, i)
    }
  }
  if (length(mm) > 0) {
    warning(paste("When computing the most probable paths, no subjects were assigned to following clusters:", paste(x$cluster_names[mm], collapse = ", ")))
  }


  if (!is.null(which.plots)) {
    if (any(!is.numeric(which.plots)) || any(!(which.plots %in% 1:x$n_clusters))) {
      stop(paste0("The which.plot argument only accepts numerical values between 1 and ", x$n_clusters, "."))
    } else if (any(which.plots %in% mm)) {
      warning("You requested cluster(s) with no subjects. Plotting only relevant clusters.")
      which.plots <- setdiff(which.plots, mm)
    }
  } else if (!ask && is.null(which.plots)) {
    which.plots <- 1:x$n_clusters
    # removing clusters with no subjects (according to hidden.paths)
    which.plots <- setdiff(which.plots, mm)
  }


  if (x$n_channels == 1) {
    x$observations <- list(x$observations)
  }
  if (ask && is.null(which.plots)) {
    tmenu <- 1:x$n_clusters
    tmenu <- setdiff(tmenu, mm)
    tmenunames <- x$cluster_names[tmenu]
    plot.new()
    repeat {
      pick <- utils::menu(tmenunames, title = "\n Select cluster (or 0 to exit):\n")
      if (pick == 0) {
        return(invisible())
      } else {
        args$x <- lapply(x$observations, \(y) y[hp_by_cluster_logic[[pick]], ])
        if (plots != "obs") {
          args$hidden.states.labels <- hidden.pathslabs[[pick]]
          args$hidden.paths <- hp_by_cluster[[pick]]
          states <- paste(
            x$cluster_names[pick],
            x$state_names[[pick]],
            sep = ":"
          )
          attr(args$hidden.paths, "alphabet") <- states
          if (attr(args$hidden.paths, "nr") %in% levels(args$hidden.paths[[1]])) {
            states <- c(states, attr(args$hidden.paths, "nr"))
          }
          if (attr(args$hidden.paths, "void") %in% levels(args$hidden.paths[[1]])) {
            states <- c(states, attr(args$hidden.paths, "void"))
          }
          args$hidden.paths[] <- lapply(
            args$hidden.paths, factor, levels = states
          )
          attr(args$hidden.paths, "labels") <- args$hidden.states.labels
          attr(args$hidden.paths, "cpal") <- hidden.pathscols[[pick]]
          args$hidden.states.colors <- hidden.pathscols[[pick]]
          if (!is.null(sortv) && sortv == "mds.hidden") {
            if (length(args$hidden.states.labels) == 1) {
              args$sortv <- "mds.obs"
            } else {
              args$sortv <- "mds.hidden"
            }
          }
        }
        args$title <- titles[tmenu[pick]]
        do.call(ssplotM, args = args)
      }
    }
  } else if (ask && !is.null(which.plots)) {
    tmenu <- which.plots
    tmenunames <- x$cluster_names[which.plots]
    plot.new()
    repeat {
      pick <- utils::menu(tmenunames, title = "\n Select cluster (or 0 to exit):\n")
      if (pick == 0) {
        return(invisible())
      } else {
        args$x <- lapply(x$observations, \(y) y[hp_by_cluster_logic[[pick]], ])
        if (plots != "obs") {
          args$hidden.states.labels <- hidden.pathslabs[[pick]]
          args$hidden.paths <- hp_by_cluster[[pick]]
          states <- paste(
            x$cluster_names[pick],
            x$state_names[[pick]],
            sep = ":"
          )
          attr(args$hidden.paths, "alphabet") <- states
          if (attr(args$hidden.paths, "nr") %in% levels(args$hidden.paths[[1]])) {
            states <- c(states, attr(args$hidden.paths, "nr"))
          }
          if (attr(args$hidden.paths, "void") %in% levels(args$hidden.paths[[1]])) {
            states <- c(states, attr(args$hidden.paths, "void"))
          }
          args$hidden.paths[] <- lapply(
            args$hidden.paths, factor, levels = states
          )
          attr(args$hidden.paths, "labels") <- args$hidden.states.labels
          attr(args$hidden.paths, "cpal") <- hidden.pathscols[[pick]]
          args$hidden.states.colors <- hidden.pathscols[[pick]]
          if (!is.null(sortv) && sortv == "mds.hidden") {
            if (length(args$hidden.states.labels) == 1) {
              args$sortv <- "mds.obs"
            } else {
              args$sortv <- "mds.hidden"
            }
          }
        }
        args$title <- titles[tmenu[pick]]
        do.call(ssplotM, args = args)
      }
    }
  } else {
    ask <- length(which.plots) > 1
    plot.new()
    for (i in which.plots) {
      args$x <- lapply(x$observations, \(y) y[hp_by_cluster_logic[[i]], ])
      if (plots != "obs") {
        args$hidden.paths <- hp_by_cluster[[i]]
        states <- paste(
          x$cluster_names[i],
          x$state_names[[i]],
          sep = ":"
        )
        attr(args$hidden.paths, "alphabet") <- states
        if (attr(args$hidden.paths, "nr") %in% levels(args$hidden.paths[[1]])) {
          states <- c(states, attr(args$hidden.paths, "nr"))
        }
        if (attr(args$hidden.paths, "void") %in% levels(args$hidden.paths[[1]])) {
          states <- c(states, attr(args$hidden.paths, "void"))
        }
        args$hidden.paths[] <- lapply(
          args$hidden.paths, factor, levels = states
        )
        args$hidden.states.labels <- hidden.pathslabs[[i]]
        attr(args$hidden.paths, "labels") <- args$hidden.states.labels
        attr(args$hidden.paths, "cpal") <- hidden.pathscols[[i]]
        args$hidden.states.colors <- hidden.pathscols[[i]]
        if (!is.null(sortv) && sortv == "mds.hidden") {
          if (length(args$hidden.states.labels) == 1) {
            args$sortv <- "mds.obs"
          } else {
            args$sortv <- "mds.hidden"
          }
        }
      }
      args$title <- titles[i]
      do.call(ssplotM, args = args)
      if (ask) {
        op <- par(ask = TRUE)
      }
    }
    # par(ask = FALSE)
  }
  invisible()
}
