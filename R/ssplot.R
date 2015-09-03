#' Stacked Plots of Multichannel Sequences and/or Most Probable 
#' Paths from Hidden Markov Models
#' 
#' Function \code{ssplot} plots stacked sequence plots of sequence object created with the \code{\link{seqdef}} function or observations and/or most probable paths of \code{hmm} objects.
#' 
#' 
#' 
#' @export
#' 
#' @param x Either hidden Markov model object of class \code{hmm} or a 
#'   sequence object created with the \code{\link{seqdef}} function or a list of
#'   sequence objects.
#'   
#'   
#' @param mpp Output from \code{\link{most_probable_path}} function. Optional, if 
#'   \code{x} is a hmm object or if \code{type=="obs"}.
#'   
#' @param plots What to plot. One of \code{"obs"} for observations (the default), 
#'   \code{"mpp"} for most probable paths, or \code{"both"} for observations 
#'   and most probable paths.
#'   
#' @param type The type of the plot. Available types are \code{"I"} for index 
#'   plots and \code{"d"} for state distribution plots. See 
#'   \code{\link{seqplot}} for details.
#'   
#' @param sortv A sorting variable or a sort method (one of \code{"from.start"},
#'   \code{"from.end"}, \code{"mds.obs"}, or \code{"mds.mpp"}) for 
#'   \code{type=="I"}. The value \code{"mds.mpp"} is only available for 
#'   \code{which="both"} and \code{which="mpp"}. Options \code{"mds.obs"} and 
#'   \code{"mds.mpp"} automatically arrange the sequences according to the 
#'   scores of multidimensional scaling (using \code{\link{cmdscale}}) for the 
#'   observed or hidden states path data from \code{\link{most_probable_path}}. 
#'   MDS scores are computed from distances/dissimilarities using a metric 
#'   defined in argument \code{dist.method}. See \code{\link{plot.stslist}} for 
#'   more details on \code{"from.start"} and \code{"from.end"}.
#'   
#' @param sort.channel The number of the channel according to which the 
#'   \code{"from.start"} or \code{"from.end"} sorting is done. Sorting according
#'   to hidden states is called with value 0. The default value is 1 (the first 
#'   channel).
#'   
#' @param dist.method The metric to be used for computing the distances of the 
#'   sequences if multidimensional scaling is used for sorting. One of "OM" 
#'   (optimal Matching, the default), "LCP" (longest common prefix), "RLCP" 
#'   (reversed LCP, i.e. longest common suffix), "LCS" (longest common 
#'   subsequence), "HAM" (Hamming distance), "DHD" (dynamic Hamming distance). 
#'   Transition rates are used for defining substitution costs if needed. See
#'   \code{\link{seqdef}} for more information on the metrics.
#'   
#' @param with.missing Controls whether missing states are included in state 
#'   distribution plots (\code{type="d"}). The default is \code{FALSE}.
#'   
#' @param title Title for the graphic. The default is \code{NA}: if 
#'   \code{title.n=TRUE}, only the number of subjects is plotted. \code{FALSE} 
#'   prints no title, even when \code{title.n=TRUE}.
#'   
#' @param title.n Controls whether the number of subjects is printed in the 
#'   title of the plot. The default is \code{TRUE}: n is plotted if \code{title}
#'   is anything but \code{FALSE}.
#'   
#' @param cex.title Expansion factor for setting the size of the font for the 
#'   title. The default value is 1. Values lesser than 1 will reduce the size of
#'   the font, values greater than 1 will increase the size.
#'   
#' @param title.pos Controls the position of the main title of the plot. The 
#'   default value is 1. Values greater than 1 will place the title higher.
#'   
#' @param withlegend Defines if and where the legend for the states is plotted. 
#'   The default value \code{"auto"} (equivalent to \code{TRUE} and 
#'   \code{right.many}) creates separate legends for each requested plot and 
#'   sets the positions automatically. Other possible values are \code{"right"},
#'   \code{"bottom"} and \code{"bottom.many"}, of which the first two create a 
#'   combined legend in the selected position and the last one creates separate 
#'   legends for each requested plot at the bottom of the graph. Value 
#'   \code{FALSE} prints no legend.
#'   
#' @param ncol.legend (A vector of) the number of columns for the legend(s). The
#'   default \code{"auto"} creates one column for each legend.
#'   
#' @param with.missing.legend If set to \code{"auto"} (the default), a legend 
#'   for the missing state is added automatically if one or more of the 
#'   sequences in the data/channel contains missing states and \code{type="I"}. 
#'   If \code{type="d"} missing states are omitted from the legends unless 
#'   \code{with.missing=TRUE}. With the value \code{TRUE} a 
#'   legend for the missing state is added in any case; equivalently 
#'   \code{FALSE} omits the legend for the missing state.
#'   
#' @param legend.prop Sets the proportion of the graphic area used for plotting 
#'   the legend when \code{withlegend} is not \code{FALSE}. The default value is
#'   0.3. Takes values from 0 to 1.
#'   
#' @param cex.legend Expansion factor for setting the size of the font for the 
#'   labels in the legend. The default value is 1. Values lesser than 1 will 
#'   reduce the size of the font, values greater than 1 will increase the size.
#'   
#' @param mpp.color A vector of colors assigned to hidden states. The default 
#'   value \code{"auto"} uses the colors assigned to the stslist object created 
#'   with \code{seqdef} if \code{mpp} is given; otherwise otherwise colors from 
#'   \code{\link{colorpalette}} are automatically used. 
#'   
#' @param mpp.labels Labels for the hidden states. The default value 
#'   \code{"auto"} uses the names provided in \code{x$state_names} if \code{x} is
#'   an hmm object; otherwise the number of the hidden state.
#'   
#' @param xaxis Controls whether an x-axis is plotted below the plot at the 
#'   bottom. The default value is \code{TRUE}.
#'   
#' @param xlab An optional label for the x-axis. If set to \code{NA}, no label 
#'   is drawn.
#'   
#' @param xtlab Optional labels for the x-axis tick labels.  If unspecified, the
#'   column names of the \code{seqdata} sequence object are used (see 
#'   \code{\link{seqdef}}).
#'   
#' @param xlab.pos Controls the position of the x axis label. The default value 
#'   is 1. Values greater to 1 will place the label further away from the plot.
#'   
#' @param ylab Labels for the channels. A vector of names for each channel 
#'   (observations). The default value \code{"auto"} uses the names provided in 
#'   \code{x$channel_names} if \code{x} is an hmm object; otherwise the 
#'   number of the channel. \code{FALSE} prints no labels.
#'   
#' @param hidden.states.title Optional label for the hidden state plot (in the 
#'   y-axis). The default is \code{"Hidden states"}.
#'   
#' @param ylab.pos Controls the position of the y axis labels (labels for 
#'   channels and/or hidden states). Either \code{"auto"} or a numerical vector 
#'   indicating on how far away from the plots the titles are positioned. The 
#'   default value \code{"auto"} positions all titles on line 1.
#'   
#' @param cex.lab Expansion factor for setting the size of the font for the axis
#'   labels. The default value is 1. Values lesser than 1 will reduce the size 
#'   of the font, values greater than 1 will increase the size.
#'   
#' @param cex.axis Expansion factor for setting the size of the font for the 
#'   axis. The default value is 1. Values lesser than 1 will reduce the size of 
#'   the font, values greater than 1 will increase the size.
#'   
#' @param ... Other arguments to be passed to \code{\link{seqplot}} to produce 
#'   the appropriate plot method.
#'   
#' @examples 
#' require(TraMineR)
#' 
#' data(biofam)
#' biofam <- biofam[1:500,]
#' 
#' #' # Single-channel data
#' 
#' biofam.seq <- seqdef(
#'   biofam[, 10:25], 
#'   states = c("Parent", "Left", "Married", "Left+Marr",
#'              "Left+Child", "Left+Marr+Child", "Divorced"),
#'   start = 15
#'   )
#' 
#' # Starting values for the emission matrix
#' B <- matrix(NA, nrow = 4, ncol = 7)
#' B[1,] <- seqstatf(biofam.seq[, 1:4])[, 2] + 0.1
#' B[2,] <- seqstatf(biofam.seq[, 5:8])[, 2] + 0.1
#' B[3,] <- seqstatf(biofam.seq[, 9:12])[, 2] + 0.1
#' B[4,] <- seqstatf(biofam.seq[, 13:15])[, 2] + 0.1
#' B <- B / rowSums(B)
#' 
#' # Starting values for the transition matrix
#' A <- matrix(c(0.80, 0.10, 0.05, 0.05,
#'               0.05, 0.80, 0.10, 0.05,
#'               0.05, 0.05, 0.80, 0.10,
#'               0.05, 0.05, 0.10, 0.80), nrow=4, ncol=4, byrow=TRUE)
#' 
#' # Starting values for initial state probabilities
#' initial_probs <- c(0.9, 0.07, 0.02, 0.01)
#' 
#' # Building a hidden Markov model with starting values
#' bHMM <- build_hmm(
#'   observations = biofam.seq, transition_matrix = A, 
#'   emission_matrix = B, initial_probs = initial_probs
#' )
#' 
#' #########################################
#' 
#' # Multichannel data
#' 
#' # Building one channel per type of event left, children or married
#' bf <- as.matrix(biofam[, 10:25])
#' children <-  bf == 4 | bf == 5 | bf == 6
#' married <- bf == 2 | bf == 3 | bf == 6
#' left <- bf == 1 | bf == 3 | bf == 5 | bf == 6
#' 
#' children[children == TRUE] <- "Children"
#' children[children == FALSE] <- "Childless"
#' 
#' married[married == TRUE] <- "Married"
#' married[married == FALSE] <- "Single"
#' 
#' left[left == TRUE] <- "Left home"
#' left[left == FALSE] <- "With parents"
#' 
#' # Building sequence objects
#' child.seq <- seqdef(children)
#' marr.seq <- seqdef(married)
#' left.seq <- seqdef(left)
#' 
#' ## Choosing colors
#' attr(child.seq, "cpal") <- c("#66C2A5", "#FC8D62")
#' attr(marr.seq, "cpal") <- c("#E7298A", "#E6AB02")
#' attr(left.seq, "cpal") <- c("#A6CEE3", "#E31A1C")
#' 
#' 
#' # Plotting state distribution plots of observations
#' ssplot(list(child.seq, marr.seq, left.seq), type = "d", plots = "obs")
#' 
#' # Plotting sequence index plots of observations
#' ssplot(
#'   list(child.seq, marr.seq, left.seq), type = "I", plots = "obs",
#'   # Sorting subjects according to the beginning of the 2nd channel (marr.seq)
#'   sortv = "from.start", sort.channel = 2, 
#'   # Controlling the size, positions, and names for channel labels
#'   ylab.pos = c(1, 2, 1), cex.lab = 1, ylab = c("Children", "Married", "Left home"), 
#'   # Plotting without legend
#'   withlegend = FALSE
#'   )
#' 
#' # Fitting hidden Markov model
#' 
#' # Initial values for emission matrices
#' B_child <- matrix(NA, nrow = 3, ncol = 2)
#' B_child[1,] <- seqstatf(child.seq[, 1:5])[, 2] + 0.1
#' B_child[2,] <- seqstatf(child.seq[, 6:10])[, 2] + 0.1
#' B_child[3,] <- seqstatf(child.seq[, 11:15])[, 2] + 0.1
#' B_child <- B_child / rowSums(B_child)
#' 
#' B_marr <- matrix(NA, nrow = 3, ncol = 2)
#' B_marr[1,] <- seqstatf(marr.seq[, 1:5])[,2] + 0.1
#' B_marr[2,] <- seqstatf(marr.seq[, 6:10])[,2] + 0.1
#' B_marr[3,] <- seqstatf(marr.seq[, 11:15])[,2] + 0.1
#' B_marr <- B_marr / rowSums(B_marr)
#' 
#' B_left <- matrix(NA, nrow = 3, ncol = 2)
#' B_left[1,] <- seqstatf(left.seq[, 1:5])[, 2] + 0.1
#' B_left[2,] <- seqstatf(left.seq[, 6:10])[, 2] + 0.1
#' B_left[3,] <- seqstatf(left.seq[, 11:15])[, 2] + 0.1
#' B_left <- B_left / rowSums(B_left)
#' 
#' # Initial values for transition matrix
#' A <- matrix(c(0.9, 0.07, 0.03,
#'                 0,  0.9,  0.1,
#'                 0,    0,    1), nrow = 3, ncol = 3, byrow = TRUE)
#' 
#' # Initial values for initial state probabilities
#' initial_probs <- c(0.9, 0.09, 0.01)
#' 
#' # Building hidden Markov model with initial parameter values
#' bHMM <- build_hmm(
#'   observations = list(child.seq, marr.seq, left.seq), 
#'   transition_matrix = A,
#'   emission_matrix = list(B_child, B_marr, B_left), 
#'   initial_probs = initial_probs
#'   )
#' 
#' # Fitting hidden Markov model
#' HMM <- fit_hmm(bHMM)
#' 
#' # Plotting observations and  the most probable paths of hidden states
#' ssplot(
#'   HMM$model, type = "I", plots = "both", 
#'   # Sorting according to multidimensional scaling of hidden states paths
#'   sortv = "mds.mpp", 
#'   ylab = c("Children", "Married", "Left home"), 
#'   # Controlling title
#'   title = "Biofam", cex.title = 1.5,
#'   # Labels for x axis and tick marks
#'   xtlab = 15:30, xlab = "Age"
#'   )
#' 
#' # Computing the most probable hidden state paths
#' mpp <- most_probable_path(HMM$model)$mpp
#' mpp.seq <- seqdef(
#'   mpp, labels = c("Hidden state 1", "Hidden state 2", "Hidden state 3")
#'   )
#' 
#' # Plotting observations and hidden state paths
#' ssplot(
#'   HMM$model, type = "I", plots = "mpp", 
#'   # Sequence object of most probable paths
#'   mpp = mpp.seq,
#'   # Sorting according to the end of hidden state paths
#'   sortv = "from.end", sort.channel = 0,
#'   # Contolling legend position, type, and proportion
#'   withlegend = "bottom", legend.prop = 0.15,
#'   # Plotting without title and y label
#'   title = FALSE, ylab = FALSE
#'   )
#'   
#' @seealso \code{\link{gridplot}} for plotting multiple ssp 
#'   objects created with \code{\link{ssp}}, \code{\link{build_hmm}} and \code{\link{fit_hmm}} for building and 
#'   fitting Hidden Markov models, and \code{\link{most_probable_path}} for 
#'   computing the most probable paths (Viterbi paths) of hidden states.


ssplot <- function(x, mpp=NULL,
                plots="obs", type="I", 
                sortv=NULL, sort.channel=1, dist.method="OM",
                with.missing=FALSE,
                title=NA, title.n=TRUE, cex.title=1, title.pos=1,
                withlegend="auto", ncol.legend="auto", 
                with.missing.legend="auto",                         
                legend.prop=0.3, cex.legend=1,
                mpp.color="auto", mpp.labels="auto",
                xaxis=TRUE, xlab=NA, xtlab=NULL, xlab.pos=1,
                ylab="auto", hidden.states.title="Hidden states", 
                ylab.pos="auto", 
                cex.lab=1, cex.axis=1, ...){
  
  sspargs <- do.call(ssp,args=as.list(match.call())[-1])
  plot.new()  
  grid.newpage()
  savepar <- par(no.readonly = TRUE)
  do.call(SSPlotter,args=sspargs)
  par(savepar)
}