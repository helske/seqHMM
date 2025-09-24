#' Stacked Plots of Multichannel Sequences and/or Most Probable
#' Paths from Hidden Markov Models
#'
#' Function `ssplot` plots stacked sequence plots of sequence object
#' created with the [seqdef()] function or observations and/or most
#' probable paths of `hmm` objects.
#' 
#' This function is deprecated and will be removed in future versions of `seqHMM`.
#'
#' @export
#'
#' @param x Either a hidden Markov model object of class `hmm` or a state
#'   sequence object of class `stslist` (created with the [TraMineR::seqdef()])
#'   function) or a list of state sequence objects.
#'
#' @param hidden.paths Output from [hidden_paths()] function. Optional, if
#'   `x` is a `hmm` object or if `type = "obs"`.
#'
#' @param plots What to plot. One of `"obs"` for observations (the default),
#'   `"hidden.paths"` for most probable paths of hidden states,
#'   or `"both"` for observations and hidden paths together.
#'
#' @param type The type of the plot. Available types are `"I"` for sequence index
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
#'   hidden paths are available. Options `"mds.obs"` and
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
#' @param title Main title for the graphic. The default is `NA`: if
#'   `title.n = TRUE`, only the number of subjects is plotted. `FALSE`
#'   prints no title, even when `title.n = TRUE`.
#'
#' @param title.n Controls whether the number of subjects (in the
#'   first channel) is printed in the title of the plot. The default is
#'   `TRUE`: n is plotted if `title` is anything but `FALSE`.
#'
#' @param cex.title Expansion factor for setting the size of the font for the
#'   title. The default value is 1. Values lesser than 1 will reduce the size of
#'   the font, values greater than 1 will increase the size.
#'
#' @param title.pos Controls the position of the main title of the plot. The
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
#'   default `"auto"` determines number of columns depending on the position of
#'   the legend.
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
#' @param yaxis Controls whether or not to plot the y-axis. The default is `FALSE`.
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
#'
#' @param ... Other arguments to be passed on to
#'   [TraMineR::seqplot()].
ssplot <- function(x, hidden.paths = NULL,
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
  
  args <- as.list(match.call())[-1]
  args[[1]] <- eval(args[[1]], envir = parent.frame())
  sspargs <- do.call(ssp, args = args)
  plot.new()
  grid.newpage()
  savepar <- par(no.readonly = TRUE)
  on.exit(savepar, add = TRUE)
  do.call(SSPlotter, args = sspargs)
  par(savepar)
}
