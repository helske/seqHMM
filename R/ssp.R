#' Define Arguments for Plotting Multichannel Sequences and/or Most Probable 
#' Paths from Hidden Markov Models
#' 
#' Function \code{ssp} defines the arguments for plotting with 
#' \code{\link{ssplot}} or \code{\link{gridplot}}.
#' 
#' 
#' 
#' @export
#' 
#' @param x Either hidden Markov model object of class \code{HMModel} or a 
#'   sequence object created with the \code{\link{seqdef}} function or a list of
#'   sequence objects.
#'   
#'   
#' @param mpp Output from \code{\link{mostProbablePath}} function. Optional, if 
#'   \code{x} is a HMModel object or if \code{type=="obs"}.
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
#'   observed or hidden states path data from \code{\link{mostProbablePath}}. 
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
#'   \code{"auto"} uses the names provided in \code{x$stateNames} if \code{x} is
#'   an HMModel object; otherwise the number of the hidden state.
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
#'   \code{x$channelNames} if \code{x} is an HMModel object; otherwise the 
#'   number of the channel. \code{FALSE} prints no labels.
#'   
#' @param hiddenStates.title Optional label for the hidden state plot (in the 
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
#' ## Building one channel per type of event left, children or married
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
#' ## Building sequence objects
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
#' # Defining the plot for state distribution plots of observations
#' ssp1 <- ssp(list(child.seq, marr.seq, left.seq), type = "d", plots = "obs")
#' # Plotting ssp1
#' plot(ssp1)
#' 
#' # Defining the plot for sequence index plots of observations
#' ssp2 <- ssp(
#'   list(child.seq, marr.seq, left.seq), type = "I", plots = "obs", 
#'   # Sorting subjects according to the beginning of the 2nd channel (marr.seq)
#'   sortv = "from.start", sort.channel = 2, 
#'   # Controlling the size, positions, and names for channel labels
#'   ylab.pos = c(1, 2, 1), cex.lab = 1, ylab = c("Children", "Married", "Left home"), 
#'   # Plotting without legend
#'   withlegend = FALSE
#'   )
#' plot(ssp2)
#' 
#' # Computing hidden Markov model
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
#' initialProbs <- c(0.9, 0.09, 0.01)
#' 
#' # Building hidden Markov model with initial parameter values
#' bHMM <- buildHMM(
#'   observations = list(child.seq, marr.seq, left.seq), 
#'   transitionMatrix = A,
#'   emissionMatrix = list(B_child, B_marr, B_left), 
#'   initialProbs = initialProbs
#'   )
#' 
#' # Fitting hidden Markov model
#' HMM <- fitHMM(bHMM)
#' 
#' # Plotting observations and hidden states (most probable) paths
#' ssp3 <- ssp(
#'   HMM$model, type = "I", plots = "both", 
#'   # Sorting according to multidimensional scaling of hidden states paths
#'   sortv = "mds.mpp", 
#'   ylab = c("Children", "Married", "Left home"), 
#'   # Controlling title
#'   title = "Biofam", cex.title = 1.5,
#'   # Labels for x axis and tick marks
#'   xtlab = 15:30, xlab = "Age"
#'   )
#' plot(ssp3)
#' 
#' # Computing the most probable paths
#' mpp <- mostProbablePath(HMM$model)$mpp
#' mpp.seq <- seqdef(
#'   mpp, labels=c("Hidden state 1", "Hidden state 2", "Hidden state 3")
#'   )
#' 
#' # Plotting observations and hidden state paths
#' ssp4 <- ssp(
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
#' plot(ssp4)
#' 
#' @return Object of class \code{ssp}.
#'   
#' @seealso \code{\link{plot.ssp}} for plotting objects created with 
#'   \code{ssp}, \code{\link{gridplot}} for plotting multiple ssp 
#'   objects, \code{\link{buildHMM}} and \code{\link{fitHMM}} for building and 
#'   fitting Hidden Markov models, and \code{\link{mostProbablePath}} for 
#'   computing the most probable paths (Viterbi paths) of hidden states.


ssp <- function(x, mpp=NULL,
                       plots="obs", type="I", 
                       sortv=NULL, sort.channel=1, dist.method="OM",
                       with.missing=FALSE,
                       title=NA, title.n=TRUE, cex.title=1, title.pos=1,
                       withlegend="auto", ncol.legend="auto", 
                       with.missing.legend="auto",                         
                       legend.prop=0.3, cex.legend=1,
                       mpp.color="auto", mpp.labels="auto",
                       xaxis=TRUE, xlab=NA, xtlab=NULL, xlab.pos=1,
                       ylab="auto", hiddenStates.title="Hidden states", 
                       ylab.pos="auto", 
                       cex.lab=1, cex.axis=1, ...){
  
  arguments <- list()
  
  if(!inherits(x,"HMModel") && (plots=="both" || plots=="mpp") && is.null(mpp)){
    stop(paste("For plotting the most probable paths, you need to add the argument mpp or give an object of class HMModel to x."))
  }
  
  if(!is.null(mpp) && !inherits(mpp,"stslist")){
    stop(paste("Object for argument mpp is not a state sequence object. Use seqdef to create one."))
  }
  
  plots <- match.arg(plots, c("both", "obs", "mpp"))
  if(withlegend!=FALSE && withlegend!=TRUE){
    withlegend <- match.arg(withlegend, c("auto", "right.many", "right",
                                          "bottom", "bottom.many"))
  }
  if(type=="I" && !is.numeric(sortv) && !is.null(sortv) && (sortv %in% c("from.start", "from.end", "mds.obs", "mds.mpp"))==FALSE){
    warning("Argument sortv only accepts values \"from.start\", \"from.end\", \"mds.mpp\" or a numerical vector (one value for each sequence).")
    sortv <- NULL
  }
  if(!is.numeric(ncol.legend) && ncol.legend!="auto"){
    warning("Argument ncol.legend only accepts values \"auto\" or a numerical vector.")
    ncol.legend <- "auto"
  }
  if(!is.numeric(ylab.pos) && ylab.pos!="auto"){
    warning("Argument ylab.pos only accepts values \"auto\" or a numerical vector.")
    ylab.pos <- "auto"
  }
  if(type!="I" && type!="d"){
    stop("Argument \"type\" should be one of \"I\" or \"d\".")
  }
  if(legend.prop<0 || legend.prop>1){
    warning("Argument legend.prop only accepts values between 0 and 1. Proportion was set to 0.3.")
    legend.prop <- 0.3
  }
  
  dist.method <- match.arg(dist.method, c("OM", "LCP", "RLCP", "LCS", "HAM", "DHD"))
  
  if(inherits(x,"HMModel")){
    obs <- x$observations
    channelNames <- x$channelNames
    if(length(ylab)>1 || (!is.na(ylab) && ylab!=FALSE)){
      if(length(ylab)==1 && ylab=="auto"){
        ylab <- x$channelNames
      }else if(length(ylab)==1 && 
                 x$numberOfChannels>1 && ylab!="auto"){
        warning("The length of ylab does not match the number of channels.")
        ylab <- rep(ylab, x$numberOfChannels)
        channelNames <- ylab
      }else if(length(ylab) < x$numberOfChannels && !is.na(ylab)){
        warning("The length of ylab does not match the number of channels.")
        ylab <- rep(ylab, x$numberOfChannels)
        channelNames <- ylab
      }else if(length(ylab) > x$numberOfChannels){
        warning("The length of ylab does not match the number of channels.")
        ylab <- ylab[1:x$numberOfChannels]
        channelNames <- ylab
      }
    }
    # Single channel stslist
  }else if(inherits(x, "stslist")){
    obs <- x
    channelNames <- 1
    if(length(ylab)>1 || (!is.na(ylab) && ylab!=FALSE)){
      if(length(ylab)==1 && ylab=="auto"){
        ylab <- 1
      }else if(length(ylab) > 1){
        warning("The length of ylab does not match the number of channels (1).")
        ylab <- ylab[1]
        channelNames <- ylab
      }
    }
    # List of stslists
  }else{
    for(i in 1:length(x)){
      if(!inherits(x[[i]], "stslist")){
        stop("At least one of the list members is not an stslist object. Use seqdef to create one or provide an object of class HMModel.")
      }
    }
    obs <- x
    channelNames <- 1:length(obs)
    if(length(ylab)>1 || (!is.na(ylab) && ylab!=FALSE)){
      if(length(ylab)==1 && ylab=="auto"){
        ylab <- 1:length(obs)
      }else if(length(ylab)==1 && 
                 length(obs)>1 && ylab!="auto"){
        warning("The length of ylab does not match the number of channels.")
        ylab <- rep(ylab, length(obs))
        channelNames <- ylab
      }else if(length(ylab) < length(obs)){
        warning("The length of ylab does not match the number of channels.")
        ylab <- rep(ylab, length(obs))
        channelNames <- ylab
      }else if(length(ylab) > length(obs)){
        warning("The length of ylab does not match the number of channels.")
        ylab <- ylab[1:length(obs)]
        channelNames <- ylab
      }
    }
  }  
  
  # Number of channels
  nchannels <- length(channelNames)
  
  # Check the number of sequences
  ncheck <- NULL
  if(plots=="both" || plots=="obs"){
    if(nchannels>1){
      for(i in 1:nchannels){
        ncheck[i] <- dim(obs[[i]])[1]
      }
    }else{
      ncheck[1] <- dim(obs)[1]
    }
  }
  if((plots=="both" || plots=="mpp") && !is.null(mpp)){
    ncheck[(length(ncheck)+1)] <- dim(mpp)[1]
  }
  if(length(unique(ncheck))>1){
    warning("The number of sequences is not the same in all channels.")
  }
  
  legend.c.prop <- 0
  legend.r.prop <- 0
  
  if(withlegend==TRUE || withlegend=="auto" || withlegend=="right"
     || withlegend=="right.many"){
    legend.c.prop <- legend.prop
  }else if(withlegend=="bottom" || withlegend=="bottom.many"){
    legend.r.prop <- legend.prop
  }
  
  # Number of plots and positions of y labels
  if(plots=="both"){
    nplots <- nchannels+1 
    if(is.character(ylab.pos)){
      if(ylab.pos=="auto"){
        ylab.pos <- rep(1,(nchannels+1))
      }else{
        stop(paste("Argument ylab.pos only accepts the value \"auto\" or a numeric vector."))
      }
    }else if(length(ylab.pos) == 1){
      ylab.pos <- rep(ylab.pos, (nchannels+1))
    }else if(length(ylab.pos)!=(nchannels+1)){
      warning(paste0("The vector provided for ylab.pos does not match the number of requested plots (",
                     (nchannels+1), ")"))
    }
  }else if(plots=="obs"){
    nplots <- nchannels
    if(is.character(ylab.pos)){
      if(ylab.pos=="auto"){
        ylab.pos <- rep(1,nchannels)
      }else{
        stop(paste("Argument ylab.pos only accepts the value \"auto\" or a numeric vector."))
      }
    }else if(length(ylab.pos) == 1){
      ylab.pos <- rep(ylab.pos, nchannels)
    }else if(length(ylab.pos)!=nchannels){
      warning(paste0("The vector provided for ylab.pos does not match the number of requested plots (",
                     nchannels, ")"))
    }
  }else if(plots=="mpp"){
    nplots <- 1
    if(is.character(ylab.pos)){
      if(ylab.pos=="auto"){
        ylab.pos <- 1
      }else{
        stop(paste("Argument ylab.pos only accepts the value \"auto\" or a numeric vector."))
      }
    }else if(length(ylab.pos)!=1){
      warning(paste0("The vector provided for ylab.pos does not match the number of requested plots (",
                     1, ")"))
    }
  }
  
  if(type == "I"){
    ylab.pos <- ylab.pos + 0.5
  }
  
  # Space for viewports
  if((is.na(title) && title.n==FALSE) || (!is.na(title) && title==FALSE)){
    title.pos <- 0
  }
  if(length(ylab)==1 && (ylab==FALSE || is.na(ylab))){
    ylab.space <- 0
  }else if(!is.na(ylab.pos) && ylab.pos!=FALSE && max(ylab.pos)<(-1)){
    ylab.space <- -1
  }else{
    ylab.space <- max(ylab.pos)*cex.lab
  }
  if(is.na(xlab) || xlab==FALSE){
    xlab.pos <- 0
  }
  if(xaxis==FALSE){
    xaxis.space <- 0
  }else{
    xaxis.space <- 1
  }
  if(length(xtlab)==1 && (is.null(xtlab) || is.na(xtlab) || xtlab==FALSE)){
    xt.space <- 0
  }else{
    xt.space <- 1
  }
  
  arguments <- list(obs=obs, nchannels=nchannels, channelNames=channelNames, nplots=nplots, 
                    legend.c.prop=legend.c.prop, legend.r.prop=legend.r.prop,
                    ylab.space=ylab.space, xaxis.space=xaxis.space, xt.space=xt.space)
  
  # Columns for legends
  if(plots=="both"){
    if(length(ncol.legend)==1 && ncol.legend=="auto"){
      if(withlegend==TRUE || withlegend=="auto" || withlegend=="right.many" ||
           withlegend=="bottom.many"){
        ncol.legend <- rep(1, (nchannels+1))
      }else if(withlegend=="right"){
        ncol.legend <- 1
      }else if(withlegend=="bottom"){
        ncol.legend <- nchannels+1
      }
    }else if((withlegend==TRUE || withlegend=="auto") && length(ncol.legend)>(nchannels+1)){
      warning(paste0("The length of ncol.legend does not match the number of requested plots. Only the first ", nchannels+1, " arguments of \"ncol.legend\" were used."))
    }else if((withlegend==TRUE || withlegend=="auto") && length(ncol.legend)<(nchannels+1)){
      warning(paste0("The length of ncol.legend does not match the number of requested plots. The last were arranged in 1 column."))
      ncol.legend[(nchannels+1-length(ncol.legend)+1):(nchannels+1)] <- 1
      ncol.legend <- ncol.legend
    }else if((withlegend=="right" || withlegend=="bottom") && 
               length(ncol.legend)>1){
      warning(paste0("The length of ncol.legend does not match the number of requested legends (1). Only the first argument of \"ncol.legend\" was used."))
      ncol.legend <- ncol.legend[1]
    }
  }else if(plots=="obs"){
    if(length(ncol.legend)==1 && ncol.legend=="auto"){
      if(withlegend==TRUE || withlegend=="auto" || withlegend=="right.many" ||
           withlegend=="bottom.many"){
        ncol.legend <- rep(1, nchannels)
      }else if(withlegend=="right"){
        ncol.legend <- 1
      }else if(withlegend=="bottom"){
        ncol.legend <- nchannels
      }
    }else if((withlegend==TRUE || withlegend=="auto") && length(ncol.legend)>nchannels){
      warning(paste0("The length of ncol.legend does not match the number of requested plots. Only the first ", nchannels, " arguments of \"ncol.legend\" were used."))
    }else if((withlegend==TRUE || withlegend=="auto") && length(ncol.legend)<nchannels){
      warning(paste0("The length of ncol.legend does not match the number of requested plots. The last were arranged in 1 column."))
      ncol.legend[(nchannels-length(ncol.legend)+1):nchannels] <- 1
      ncol.legend <- ncol.legend
    }else if((withlegend=="right" || withlegend=="bottom") && 
               length(ncol.legend)>1){
      warning(paste0("The length of ncol.legend does not match the number of requested legends (1). Only the first argument of \"ncol.legend\" was used."))
      ncol.legend <- ncol.legend[1]
    }
  }else if(plots=="mpp"){
    if(length(ncol.legend)==1 && ncol.legend=="auto"){
      ncol.legend <- 1
    }else if((withlegend==TRUE || withlegend=="auto") && length(ncol.legend)>1){
      warning(paste0("The length of ncol.legend does not match the number of requested plots. Only the first argument of \"ncol.legend\" was used."))
    }else if((withlegend=="right" || withlegend=="bottom") && 
               length(ncol.legend)>1){
      warning(paste0("The length of ncol.legend does not match the number of requested legends (1). Only the first argument of \"ncol.legend\" was used."))
      ncol.legend <- ncol.legend[1]
    }
  }  
  
  # Most probable paths
  if(plots=="both" || plots=="mpp" || (plots=="obs" && !is.null(mpp))){
    if(!is.null(mpp)){
      if(!is.null(mpp.labels) && length(mpp.labels)==1 && mpp.labels=="auto"){
        mpp.labels <- attr(mpp, "labels")
      }else if(length(mpp.labels)==1 && is.null(mpp.labels)){
        mpp.labels <- rep("", length(alphabet(mpp)))
      }else if(!is.null(mpp.labels) && length(mpp.labels)!=length(alphabet(mpp))){
        warning("The number of labels for hidden states does not match the number of hidden states. Given labels were not used.")
        mpp.labels <- attr(mpp, "labels")
      }
    }
    # Computing mpp
    if(is.null(mpp)){
      mpp <- suppressMessages(mostProbablePath(x)$mpp)
      if(length(mpp.labels)==1 && is.null(mpp.labels)){
        mpp.labels <- rep("", length(alphabet(mpp)))
      }
      if(length(mpp.labels)==1 && mpp.labels=="auto"){
        mpp.labels <- paste("State", alphabet(mpp))
      }else if(!is.null(mpp.labels) && length(mpp.labels)!=length(alphabet(mpp))){
        warning("The number of labels for hidden states does not match the number of hidden states. Labels were not used.")
        mpp.labels <- paste("State", alphabet(mpp))
      }
    }
    if(nchannels>1){
      mpp.seq <- suppressWarnings(suppressMessages(seqdef(mpp, 
                                         missing.color=attr(obs[[1]],"missing.color"),
                                         labels=mpp.labels,
                                         start=attr(obs[[1]],"start"),
                                         xtstep=attr(obs[[1]],"xtstep"))))
    }else{
      mpp.seq <- suppressWarnings(suppressMessages(seqdef(mpp, 
                                         missing.color=attr(obs,"missing.color"),
                                         labels=mpp.labels,
                                         start=attr(obs,"start"),
                                         xtstep=attr(obs,"xtstep"))))
    }
    # Color palette for mpp
    if(length(mpp.color)==1 && mpp.color=="auto" && length(mpp)>1){
      if(is.null(attr(mpp, "cpal"))){
        attr(mpp.seq, "cpal") <- seqHMM::colorpalette[[length(alphabet(mpp.seq))]]
      }else{
        attr(mpp.seq, "cpal") <- attr(mpp, "cpal")
      }
    }else if(!is.null(mpp.labels) && length(mpp.labels)==1 && mpp.color=="auto" && length(mpp)==1 && is.null(mpp)){
      attr(mpp.seq, "cpal") <- seqHMM::colorpalette[[length(alphabet(mpp.seq))]]
    }else if(all(isColor(mpp.color))){
      if(length(mpp.color)!=length(alphabet(mpp.seq))){
        warning(paste0("Number of colors assigned to mpp.color does not match the number of hidden states. \n
                       There are ", length(alphabet(mpp.seq)), 
                       " hidden states but ", length(mpp.color), " color(s)."))
      }
      attr(mpp.seq, "cpal") <- rep(mpp.color, length(alphabet(mpp.seq)))[1:length(alphabet(mpp.seq))]
    }else{
      stop(paste("Please provide a vector of colors for argument mpp.color or use value \"auto\" to automatically determine gray scale colors."))
    }
    arguments <- c(arguments, list(mpp.seq=mpp.seq))
    # Sorting sequences according to multidimensional scaling score of mpp
    if(!is.null(sortv)){
      if(length(sortv)==1 && sortv=="mds.mpp"){
        dist.mpp <- suppressMessages(seqdist(mpp.seq, method=dist.method, 
                                             sm="TRATE", with.missing=TRUE))
        sortv <- cmdscale(dist.mpp, k=1)
        #         arguments <- c(arguments, list(sortv=sortv))
      }
    } 
    }
  
  # Ordering sequences for sortv
  if(type=="I" && !is.null(sortv)){
    if(nchannels>1){
      # Multidimensional scaling on observations
      if(length(sortv)==1 && sortv=="mds.obs"){
        dist.obs <- suppressMessages(seqdistmc(obs, method=dist.method, 
                                               sm="TRATE", with.missing=TRUE))
        sortv <- cmdscale(dist.obs, k=1)
        #         arguments <- c(arguments, list(sortv=sortv))
      }
      # Sorting from start or end
      if(length(sortv)==1 && (sortv=="from.start" || sortv=="from.end")){
        end <- if (sortv == "from.end") {
          if(sort.channel>0 && sort.channel<=length(obs)){
            max(seqlength(obs[[sort.channel]]))
          }else if(sort.channel==0){
            if(plots=="both" || plots=="mpp"){
              max(seqlength(mpp.seq))
            }else{
              warning("Most probable paths are only computed automatically for argument plots=\"both\" or plots=\"mpp\". Sequences were not sorted.")
              sortv <- NULL
            }
          }else{
            warning(paste0("For data with ", length(obs), " channels, the value for sort.channel must be a non-negative integer smaller or equal to ", length(obs), ". Sequences were not sorted."))
            sortv <- NULL
          }
        }else {
          1
        }
        beg <- if (sortv == "from.end") {
          1
        }else{
          if(sort.channel>0 && sort.channel<=length(obs)){
            max(seqlength(obs[[sort.channel]]))
          }else if(sort.channel==0){
            if(plots=="both" || plots=="mpp"){
              max(seqlength(mpp.seq))
            }else{
              warning("Most probable paths are only computed automatically for argument plots=\"both\" or plots=\"mpp\". Sequences were not sorted.")
              sortv <- NULL
            }
          }else{
            warning(paste0("For data with ", length(obs), " channels, the value for sort.channel must be a non-negative integer smaller or equal to ", length(obs), ". Sequences were not sorted."))
            sortv <- NULL
          }
        }
        if(sort.channel>0 && sort.channel<=length(obs)){
          orderv <- do.call(order, as.data.frame(obs[[sort.channel]])[, end:beg])
          arguments <- c(arguments, list(orderv=orderv))
        }else if(sort.channel==0 && plots!="obs"){
          orderv <- do.call(order, as.data.frame(mpp.seq)[, end:beg])
          arguments <- c(arguments, list(orderv=orderv))
        }
      }
    }else{
      if(length(sortv)==1 && sortv=="mds.obs"){
        dist.obs <- suppressMessages(seqdist(obs, method=dist.method, 
                                             sm="TRATE", with.missing=TRUE))
        sortv <- cmdscale(dist.obs, k=1)
        #         arguments <- c(arguments, list(sortv=sortv))
      }else if(length(sortv)==1 && (sortv=="from.start" || sortv=="from.end")){
        end <- if (sortv == "from.end") {
          if(sort.channel==1){
            max(seqlength(obs))
          }else if(sort.channel==0){
            if(plots=="both" || plots=="mpp"){
              max(seqlength(mpp.seq))
            }else{
              warning("Most probable paths are only computed automatically for argument plots=\"both\" or plots=\"mpp\". Sequences were not sorted.")
              sortv <- NULL
            }
          }else{
            warning(paste0("For data with 1 channel, the value for sort.channel must be 0 (for most probable paths) or 1 (for observations). Sequences were not sorted."))
            sortv <- NULL
          }
        }else {
          1
        }
        beg <- if (sortv == "from.end") {
          1
        }else{
          if(sort.channel==1){
            max(seqlength(obs))
          }else if(sort.channel==0){
            if(plots=="both" || plots=="mpp"){
              max(seqlength(mpp.seq))
            }else{
              warning("Most probable paths are only computed automatically for argument plots=\"both\" or plots=\"mpp\". Sequences were not sorted.")
              sortv <- NULL
            }
          }else{
            warning(paste0("For data with 1 channel, the value for sort.channel must be 0 (for most probable paths) or 1 (for observations). Sequences were not sorted."))
            sortv <- NULL
          }
        }
        if(sort.channel==1){
          orderv <- do.call(order, as.data.frame(obs[[sort.channel]])[, end:beg])
          arguments <- c(arguments, list(orderv=orderv))
        }else if(sort.channel==0 && plots!="obs"){
          orderv <- do.call(order, as.data.frame(mpp.seq)[, end:beg])
          arguments <- c(arguments, list(orderv=orderv))
        }
      }
    }
  }
  
  
  # Plotting observations
  if(plots=="obs" && xaxis==TRUE){
    plotxaxis <- TRUE
  }else{
    plotxaxis <- FALSE
  }  
  arguments <- c(arguments, list(plotxaxis=plotxaxis))
  if(plots=="both" || plots=="obs"){ 
    if(type=="I"){ 
      if(length(sortv)==1 && sortv=="mds.mpp" && plots=="obs" && is.null(mpp)){
        warning("Most probable paths are only computed automatically for argument plots=\"both\" or plots=\"mpp\". Sequences were not sorted.")
        sortv <- NULL
      }
    }
  }
  
  # Legends
  if(type=="d"){
    if(with.missing==FALSE && with.missing.legend!=FALSE && with.missing.legend!=TRUE &&
                                 with.missing.legend=="auto"){
      with.missing.legend <- FALSE
    }
  }
  
  if(length(list(...))==0){
    arguments <- c(arguments, list(mpp=mpp, plots=plots, type=type, 
                                   sortv=sortv, sort.channel=sort.channel, 
                                   with.missing=with.missing,
                                   title=title, title.n=title.n, cex.title=cex.title, title.pos=title.pos,
                                   withlegend=withlegend, ncol.legend=ncol.legend, 
                                   with.missing.legend=with.missing.legend,                         
                                   legend.prop=legend.prop, cex.legend=cex.legend,
                                   mpp.color=mpp.color, mpp.labels=mpp.labels,
                                   xaxis=xaxis, xlab=xlab, xtlab=xtlab, xlab.pos=xlab.pos,
                                   ylab=ylab, hiddenStates.title=hiddenStates.title, 
                                   ylab.pos=ylab.pos, 
                                   cex.lab=cex.lab, cex.axis=cex.axis, call=match.call()))
  }else{
    arguments <- c(arguments, list(mpp=mpp, plots=plots, type=type, 
                                   sortv=sortv, sort.channel=sort.channel, 
                                   with.missing=with.missing,
                                   title=title, title.n=title.n, cex.title=cex.title, title.pos=title.pos,
                                   withlegend=withlegend, ncol.legend=ncol.legend, 
                                   with.missing.legend=with.missing.legend,                         
                                   legend.prop=legend.prop, cex.legend=cex.legend,
                                   mpp.color=mpp.color, mpp.labels=mpp.labels,
                                   xaxis=xaxis, xlab=xlab, xtlab=xtlab, xlab.pos=xlab.pos,
                                   ylab=ylab, hiddenStates.title=hiddenStates.title, 
                                   ylab.pos=ylab.pos, 
                                   cex.lab=cex.lab, cex.axis=cex.axis, call=match.call()),
                   list(...))    
  }
  
  class(arguments) <- "ssp"
  arguments
  
}
