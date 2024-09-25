#' Stacked Sequence Plots of Multichannel Sequences and/or Most Probable
#' Paths from Hidden Markov Models
#'
#' Function `stacked_sequence_plot` draws stacked sequence plots of sequence 
#' object created with the [TraMineR::seqdef] function or observations and/or most
#' probable paths of model objects of `seqHMM` (e.g., `hmm` and `mhmm`).
#'
#' @param x Either a hidden Markov model object of class `hmm`, `mhmm`, `nhmm`, 
#' or `mnhmm`, a sequence object of class `stslist` (created with the 
#' [TraMineR::seqdef()] function) or a list of `stslist` objects.
#' @param plots What to plot. One of `"obs"` for observations (the default),
#' `"hidden_paths"` for most probable paths of hidden states,
#' or `"both"` for observations and hidden paths together. Latter two options 
#' are only possible for model objects.
#' @param type The type of the plot. Available types are `"index"` for sequence 
#' index plots and `"distribution"` for state distribution plots (the default). 
#' See [ggseqplot::ggseqiplot()] and [ggseqplot::ggseqdplot()] for details.
#' @param ids Indexes of the subjects to be plotted (the default is all). For 
#' example, `ids = c(1:10, 15) plots the first ten subjects and subject 15 in 
#' the data.
#' @param sort_by A sorting variable or a sort method (one of `"none`, `"start"`,
#' `"end"`, or `"mds"` for `type = "index"`. Option `"mds"` arranges the 
#' sequences according to the scores of multidimensional scaling (using 
#' [stats::cmdscale()]). Default is `"none"`, i.e., no sorting. Numeric vectors 
#' are passed to `sortv` argument of [ggseqplot::ggseqiplot()].
#' @param sort_channel Name of the channel which should be used for the 
#' sorting. Alternatively value `"Hidden states"` uses the hidden state 
#' sequences for sorting. Default is to sort by the first channel in the data. 
#' If `sort_by = "mds"`, all channels are used for defining the sorting.
#' @param dist_method The metric to be used for computing the distances of the
#' sequences if multidimensional scaling is used for sorting. One of `"OM"`
#' (optimal matching, the default), `"LCP"` (longest common prefix), `"RLCP"`
#' (reversed LCP, i.e. longest common suffix), `"LCS"` (longest common
#' subsequence), `"HAM"` (Hamming distance), and `"DHD"` (dynamic Hamming 
#' distance). Transition rates are used for defining substitution costs if 
#' needed. See [TraMineR::seqdef()] for more information on the metrics.
#' @param group Variable used for grouping the sequences in each channel, which 
#' is passed to [ggseqplot::ggseqiplot()] and [ggseqplot::ggseqdplot()]. By 
#' default, no grouping is done, except for mixture models where the grouping 
#' is based on most probable clusters (defined by the most probable hidden 
#' paths). Grouping by clusters can be overloaded by supplying variable for 
#' `group` or by setting `group = NA`.
#' @param legend_position Position of legend for each channel, 
#' passed to `legend.position` argument of [ggplot2::theme()]. Either a vector 
#' of length 1, or of length matching the number of channels to be plotted.
#' @param ... Other arguments to [ggseqplot::ggseqiplot()] or 
#' [ggseqplot::ggseqdplot()].
#' @export
#' @examples
#' p <- stacked_sequence_plot(
#'   mhmm_biofam, 
#'   plots = "both", 
#'   type = "d", 
#'   legend_position = c("right", "right", "right", "none")
#' )
#' library(ggplot2)
#' p & theme(plot.margin = unit(c(1, 1, 0, 2), "mm"))
#' 
stacked_sequence_plot <- function(
    x, plots = "obs", type = "d", ids,
    sort_by = "none", sort_channel, dist_method = "OM", group = NULL, 
    legend_position = "right", ...) {
  
  plots <- match.arg(plots, c("obs", "hidden_paths", "both"))
  
  if (inherits(x, c("hmm", "nhmm", "mhmm", "mnhmm"))) {
    if (is.null(group) && inherits(x, c("mhmm", "mnhmm"))) {
      hp <- hidden_paths(x)
      group <- factor(
        most_probable_cluster(x, type = "viterbi", hp = hp),
        levels = x$cluster_names
      )
    } else {
      if (plots != "obs") {
        hp <- hidden_paths(x)
      }
    }
    if (plots == "both") {
      if (x$n_channels == 1) {
        y <- list(x$observations, hp)
      } else {
        y <- c(x$observations, list(hp))
      }
      channel_names <- c(x$channel_names, "Hidden states")
      n_channels <- x$n_channels + 1
    }
    if (plots == "hidden_paths") {
      y <- hp
      channel_names <- "Hidden states"
      n_channels <- 1
    }
    if (plots == "obs") {
      y <- x$observations
      channel_names <- x$channel_names
      n_channels <- x$n_channels
    }
  } else {
    stopifnot_(
      plots == "obs",
      "Cannot draw most probable hidden paths as {.arg x} is not a model object."
    )
    y <- x
    if (TraMineR::is.stslist(y)) {
      n_channels <- 1
      channel_names <- "Observations"
    } else {
      stopifnot_(
        all(unlist(lapply(y, inherits, "stslist"))),
        "{.arg x} is not an {.cls stslist} object or a list of such objects."
      )
      n_channels <- length(y)
      if (is.null(channel_names <- names(y))) {
        channel_names <- paste("Channel", seq_len(n_channels))
      }
    }
  }
  if (!missing(ids)) {
    if (n_channels == 1) {
      y <- y[ids, ]
    } else {
      for (i in seq_len(n_channels)) {
        y[[i]] <- y[[i]][ids, ]
      }
    }
  }
  if (!identical(sort_by, "none")) {
    if (length(sort_by) == 1) {
      sort_by <- match.arg(sort_by, c("none", "start", "end", "mds"))
      if (missing(sort_channel)) sort_channel <- channel_names[1]
      stopifnot_(
        sort_channel %in% channel_names || sort_channel %in% seq_len(n_channels),
        paste0("{.arg sort_channel} should be either ",
        "{cli::cli_vec(channel_names, style = list('vec-last' = ' or '))}, ",
        "or integer between 1 and {n_channels}."
        )
      )
      y <- sort_sequences(y, sort_by, sort_channel, dist_method)
      sort_by <- NULL
    } else {
      stopifnot_(
        length(sort_by) == nrow(y),
        "Length of {.arg sort_by} does not match the number of sequences 
        {nrow(y)}."
      )
    }
  }
  if (identical(group, NA)) group <- NULL
  
  type <- match.arg(type, c("distribution", "index"))
  
  if (n_channels == 1) {
    if (type == "distribution") {
      cpal_y <- setNames(attr(y, "cpal"), attr(y, "labels"))
      p <- ggseqplot::ggseqdplot(y, group = group, sortv = sort_by, ...) + 
        ggplot2::theme(legend.position = legend_position) +
        ggplot2::ylab("Proportion") +
        ggplot2::xlab("Time")
      suppressMessages(
        p <- p +
          ggplot2::scale_fill_manual(values = cpal_y)
      )
    }
    if (type == "index") {
      cpal_y <- setNames(attr(y, "cpal"), attr(y, "labels"))
      p <- ggseqplot::ggseqiplot(y, group = group, sortv = sort_by, ...) + 
        ggplot2::theme(legend.position = legend_position) +
        ggplot2::ylab("Sequence") +
        ggplot2::xlab("Time")
      suppressMessages(
        p <- p +
          ggplot2::scale_fill_manual(values = cpal_y)
      )
      if (!isTRUE(list(...)$border)) {
        suppressMessages(
          p <- p +
            ggplot2::scale_colour_manual(values = cpal_y)
        )
      }
    }
  } else {
    if (length(legend_position) == 1) {
      legend_position <- rep(legend_position, n_channels)
    }
    stopifnot_(
      length(legend_position) == n_channels,
      "Length of {.arg legend_position} should be 1 or 
      {n_channels} (number of channels to be drawn)."
    )
    p <- vector("list", n_channels)
    if (type == "distribution") {
      for (i in seq_len(n_channels)) {
        cpal_y <- setNames(attr(y[[i]], "cpal"), attr(y[[i]], "labels"))
        p[[i]] <- ggseqplot::ggseqdplot(y[[i]], group = group, ...) + 
          ggplot2::theme(legend.position = legend_position[i]) +
          ggplot2::ggtitle(channel_names[i]) +
          ggplot2::ylab("Proportion") +
          ggplot2::xlab("Time")
        suppressMessages(
          p[[i]] <- p[[i]] + 
            ggplot2::scale_fill_manual(values = cpal_y)
        )
      }
    }
    if (type == "index") {
      for (i in seq_len(n_channels)) {
        cpal_y <- setNames(attr(y[[i]], "cpal"), attr(y[[i]], "labels"))
        p[[i]] <- ggseqplot::ggseqiplot(y[[i]], group = group, ...) + 
          ggplot2::theme(legend.position = legend_position[i]) +
          ggplot2::ggtitle(channel_names[i]) +
          ggplot2::ylab("Sequence") +
          ggplot2::xlab("Time")
        suppressMessages(
          p[[i]] <- p[[i]] + 
            ggplot2::scale_fill_manual(values = cpal_y)
        )
        if (!isTRUE(list(...)$border)) {
          suppressMessages(
            p[[i]] <- p[[i]] +
              ggplot2::scale_colour_manual(values = cpal_y)
          )
        }
      }
    }
    p <- patchwork::wrap_plots(p, ncol = 1, ...)
  }
  p
}
