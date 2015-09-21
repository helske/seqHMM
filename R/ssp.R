#' Define Arguments for Plotting Multichannel Sequences and/or Most Probable 
#' Paths from Hidden Markov Models
#' 
#' Function \code{ssp} defines the arguments for plotting with 
#' \code{\link{ssplot}} or \code{\link{gridplot}}.
#' 
#' 
#' 
#' @export
#' @param x Either hidden Markov model object of class \code{hmm} or a 
#'   sequence object created with the \code{\link{seqdef}} function or a list of
#'   sequence objects.
#'
#' @param hidden.paths Output from \code{\link{hidden_paths}} function. Optional, if 
#'   \code{x} is a hmm object or if \code{type=="obs"}.
#'   
#' @param plots What to plot. One of \code{"obs"} for observations (the default), 
#'   \code{"hidden.paths"} for most probable paths, or \code{"both"} for observations 
#'   and most probable paths.
#'   
#' @param type The type of the plot. Available types are \code{"I"} for index 
#'   plots and \code{"d"} for state distribution plots (the default). See 
#'   \code{\link{seqplot}} for details.
#'   
#' @param sortv A sorting variable or a sort method (one of \code{"from.start"},
#'   \code{"from.end"}, \code{"mds.obs"}, or \code{"mds."}) for 
#'   \code{type=="I"}. The value \code{"mds.hidden"} is only available for 
#'   \code{which="both"} and \code{which="hidden.paths"}. Options \code{"mds.obs"} and 
#'   \code{"mds.hidden"} automatically arrange the sequences according to the 
#'   scores of multidimensional scaling (using \code{\link{cmdscale}}) for the 
#'   observed or hidden states path data from \code{\link{hidden_paths}}. 
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
#' @param hidden.states.colors A vector of colors assigned to hidden states. The default 
#'   value \code{"auto"} uses the colors assigned to the stslist object created 
#'   with \code{seqdef} if \code{hidden.paths} is given; otherwise otherwise colors from 
#'   \code{\link{colorpalette}} are automatically used. 
#'   
#' @param hidden.states.labels Labels for the hidden states. The default value 
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
#'   \code{x$channel_names} if \code{x} is an \code{hmm} object; otherwise the 
#'   names of the list in \code{x} if given, or the
#'   number of the channel if names are not given. \code{FALSE} prints no labels.
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
#' 
#' \dontrun{
#' 
#' data(biofam3c)
#' 
#' ## Building sequence objects
#' child.seq <- seqdef(biofam3c$children, start = 15)
#' marr.seq <- seqdef(biofam3c$married, start = 15)
#' left.seq <- seqdef(biofam3c$left, start = 15)
#' 
#' ## Choosing colors
#' attr(child.seq, "cpal") <- c("#66C2A5", "#FC8D62")
#' attr(marr.seq, "cpal") <- c("#AB82FF", "#E6AB02", "#E7298A")
#' attr(left.seq, "cpal") <- c("#A6CEE3", "#E31A1C")
#' 
#' 
#' # Defining the plot for state distribution plots of observations
#' ssp1 <- ssp(list(child.seq, marr.seq, left.seq))
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
#' # Plotting hidden Markov models
#' 
#' # Loading data
#' data(hmm_biofam)
#' 
#' # Plotting observations and hidden states (most probable) paths
#' ssp3 <- ssp(
#'   hmm_biofam, type = "I", plots = "both", 
#'   # Sorting according to multidimensional scaling of hidden states paths
#'   sortv = "mds.hidden", 
#'   ylab = c("Children", "Married", "Left home"), 
#'   # Controlling title
#'   title = "Biofam", cex.title = 1.5,
#'   # Labels for x axis and tick marks
#'   xtlab = 15:30, xlab = "Age"
#'   )
#' plot(ssp3)
#' 
#' # Computing the most probable paths of hidden states
#' hid <- hidden_paths(hmm_biofam)
#' alphabet(hid) <- paste("Hidden state", 1:4)
#' 
#' # Plotting observations and hidden state paths
#' ssp4 <- ssp(
#'   hmm_biofam, type = "I", plots = "hidden.paths", 
#'   # Sequence object of most probable paths
#'   hidden.paths = hid,
#'   # Sorting according to the end of hidden state paths
#'   sortv = "from.end", sort.channel = 0,
#'   # Contolling legend position, type, and proportion
#'   withlegend = "bottom", legend.prop = 0.15,
#'   # Plotting without title and y label
#'   title = FALSE, ylab = FALSE
#'   )
#' plot(ssp4)
#' }
#' @return Object of class \code{ssp}.
#'   
#' @seealso \code{\link{plot.ssp}} for plotting objects created with 
#'   the \code{ssp} function; \code{\link{gridplot}} for plotting multiple \code{ssp} 
#'   objects; \code{\link{build_hmm}} and \code{\link{fit_hmm}} for building and 
#'   fitting hidden Markov models; \code{\link{hidden_paths}} for 
#'   computing the most probable paths of hidden states; and \code{\link{biofam3c}}
#'   \code{\link{hmm_biofam}} for information on the data and model used in the example.


ssp <- function(x, hidden.paths = NULL,
                       plots = "obs", type = "d", 
                       sortv = NULL, sort.channel = 1, dist.method = "OM",
                       with.missing = FALSE,
                       title = NA, title.n = TRUE, cex.title = 1, title.pos = 1,
                       withlegend = "auto", ncol.legend = "auto", 
                       with.missing.legend = "auto",                         
                       legend.prop = 0.3, cex.legend = 1,
                       hidden.states.colors = "auto", hidden.states.labels = "auto",
                       xaxis = TRUE, xlab = NA, xtlab = NULL, xlab.pos = 1,
                       ylab = "auto", hidden.states.title = "Hidden states", 
                       ylab.pos = "auto", 
                       cex.lab = 1, cex.axis = 1, ...){
  
  arguments <- list()
  
  if(!inherits(x,"hmm") && (plots=="both" || plots=="hidden.paths") && is.null(hidden.paths)){
    stop(paste("For plotting the most probable paths, you need to add the argument hidden.paths or give an object of class hmm to x."))
  }
  
  if(!is.null(hidden.paths) && !inherits(hidden.paths,"stslist")){
    stop(paste("Object for argument hidden.paths is not a state sequence object. Use seqdef to create one."))
  }
  
  plots <- match.arg(plots, c("both", "obs", "hidden.paths"))
  if(withlegend!=FALSE && withlegend!=TRUE){
    withlegend <- match.arg(withlegend, c("auto", "right.many", "right",
                                          "bottom", "bottom.many"))
  }
  if(type=="I" && !is.numeric(sortv) && !is.null(sortv) && (sortv %in% c("from.start", "from.end", "mds.obs", "mds.hidden"))==FALSE){
    warning("Argument sortv only accepts values \"from.start\", \"from.end\", \"mds.hidden\" or a numerical vector (one value for each sequence).")
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
  
  if (inherits(x,"hmm")) {
    obs <- x$observations
    channel_names <- x$channel_names
    if (length(ylab) > 1 || (!is.na(ylab) && ylab != FALSE)) {
      if (plots != "hidden.paths") {
        if (length(ylab)==1 && ylab=="auto") {
          ylab <- x$channel_names
        } else {
          ylab <- rep(ylab, length.out = x$n_channels)
          channel_names <- ylab
        }
      } else {
        if(!is.null(ylab) && hidden.states.title == "Hidden states"){
          warning("Argument ylab only modifies channel titles (observations). Did you mean to change hidden.states.title?")
        }
      }
    }
    # Single channel stslist
  } else if(inherits(x, "stslist")) {
    obs <- x
    channel_names <- "Observations"
    if (length(ylab) > 1 || (!is.na(ylab) && ylab != FALSE)) {
      if (length(ylab) == 1 && ylab == "auto") {
        ylab <- "Observations"
      } else if(length(ylab) > 1) {
        ylab <- ylab[1]
        channel_names <- ylab
      }
    }
    # List of stslists
  } else {
    for (i in 1:length(x)) {
      if(!inherits(x[[i]], "stslist")){
        stop("At least one of the list members is not an stslist object. Use seqdef to create one or provide an object of class hmm.")
      }
    }
    obs <- x
    if (!is.null(names(x))) {
      channel_names <- names(x)
    } else {
      channel_names <- 1:length(obs)
    }
    if (length(ylab) > 1 || (!is.na(ylab) && ylab != FALSE)) {
      if (length(ylab)==1 && ylab=="auto") {
        ylab <- channel_names
      } else {
        ylab <- rep(ylab, length.out = length(obs))
        channel_names <- ylab
      }
    }
  }  
  
  # Number of channels
  nchannels <- length(channel_names)
  
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
  if((plots=="both" || plots=="hidden.paths") && !is.null(hidden.paths)){
    ncheck[(length(ncheck)+1)] <- dim(hidden.paths)[1]
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
  if (plots=="both") {
    nplots <- nchannels+1 
    if (is.character(ylab.pos)) {
      if (ylab.pos=="auto") {
        ylab.pos <- rep(1, (nchannels + 1) )
      } else {
        stop(paste("Argument ylab.pos only accepts the value \"auto\" or a numeric vector."))
      }
    }else if (length(ylab.pos) != (nchannels+1) ){
      ylab.pos <- rep(ylab.pos, length.out = (nchannels+1))
    }
  }else if (plots == "obs") {
    nplots <- nchannels
    if (is.character(ylab.pos)) {
      if (ylab.pos == "auto") {
        ylab.pos <- rep(1, nchannels)
      } else {
        stop(paste("Argument ylab.pos only accepts the value \"auto\" or a numeric vector."))
      }
    } else if(length(ylab.pos) != nchannels) {
      ylab.pos <- rep(ylab.pos, length.out = nchannels)
    }
  }
  
  if (type == "I") {
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
  xaxis.space <- ifelse(xaxis, 1, 0)
  if(length(xtlab)==1 && (is.null(xtlab) || is.na(xtlab) || xtlab==FALSE)){
    xt.space <- 0
  }else{
    xt.space <- 1
  }
  
  arguments <- list(obs=obs, nchannels=nchannels, channel_names=channel_names, nplots=nplots, 
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
      ncol.legend <- ncol.legend[1:(nchannels + 1)]
    }else if((withlegend==TRUE || withlegend=="auto") && length(ncol.legend)<(nchannels+1)){
      ncol.legend <- rep(ncol.legend, length.out = (nchannels + 1))
    }else if((withlegend=="right" || withlegend=="bottom") && 
               length(ncol.legend)>1){
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
      ncol.legend <- ncol.legend[1:nchannels]
    }else if((withlegend==TRUE || withlegend=="auto") && length(ncol.legend)<nchannels){
      ncol.legend <- rep(ncol.legend, length.out = nchannels)
      ncol.legend <- ncol.legend
    }else if((withlegend=="right" || withlegend=="bottom") && 
               length(ncol.legend)>1){
      ncol.legend <- ncol.legend[1]
    }
  }else if(plots=="hidden.paths"){
    if(length(ncol.legend)==1 && ncol.legend=="auto"){
      ncol.legend <- 1
    }else if((withlegend==TRUE || withlegend=="auto" || 
              withlegend=="right" || withlegend=="bottom") && length(ncol.legend)>1){
      ncol.legend <- ncol.legend[1]
    }
  }  
  
  # Most probable paths
  if(plots=="both" || plots=="hidden.paths" || (plots=="obs" && !is.null(hidden.paths))){
    if(!is.null(hidden.paths)){
      if(!is.null(hidden.states.labels) && length(hidden.states.labels)==1 && hidden.states.labels=="auto"){
        hidden.states.labels <- attr(hidden.paths, "labels")
      }else if(length(hidden.states.labels)==1 && is.null(hidden.states.labels)){
        hidden.states.labels <- rep("", length(alphabet(hidden.paths)))
      }else if(!is.null(hidden.states.labels) && length(hidden.states.labels)!=length(alphabet(hidden.paths))){
        warning("The number of labels for hidden states does not match the number of hidden states. Given labels were not used.")
        hidden.states.labels <- attr(hidden.paths, "labels")
      }
    }
    # Computing hidden.paths
    if(is.null(hidden.paths)){
      hidden.paths <- suppressMessages(hidden_paths(x))
      if(length(hidden.states.labels)==1 && is.null(hidden.states.labels)){
        hidden.states.labels <- rep("", length(alphabet(hidden.paths)))
      }
      if(length(hidden.states.labels)==1 && hidden.states.labels=="auto"){
        hidden.states.labels <- paste("State", alphabet(hidden.paths))
      }else if(!is.null(hidden.states.labels) && length(hidden.states.labels)!=length(alphabet(hidden.paths))){
        warning("The number of labels for hidden states does not match the number of hidden states. Labels were not used.")
        hidden.states.labels <- paste("State", alphabet(hidden.paths))
      }
    }
    if(nchannels>1){
      hidden.paths.seq <- suppressWarnings(suppressMessages(seqdef(hidden.paths, 
                                         missing.color=attr(obs[[1]],"missing.color"),
                                         labels=hidden.states.labels,
                                         start=attr(obs[[1]],"start"),
                                         xtstep=attr(obs[[1]],"xtstep"))))
    }else{
      hidden.paths.seq <- suppressWarnings(suppressMessages(seqdef(hidden.paths, 
                                         missing.color=attr(obs,"missing.color"),
                                         labels=hidden.states.labels,
                                         start=attr(obs,"start"),
                                         xtstep=attr(obs,"xtstep"))))
    }
    # Color palette for hidden.paths
    if(length(hidden.states.colors)==1 && hidden.states.colors=="auto" && length(hidden.paths)>1){
      if(is.null(attr(hidden.paths, "cpal"))){
        attr(hidden.paths.seq, "cpal") <- seqHMM::colorpalette[[length(alphabet(hidden.paths.seq))]]
      }else{
        attr(hidden.paths.seq, "cpal") <- attr(hidden.paths, "cpal")
      }
    }else if(!is.null(hidden.states.labels) && length(hidden.states.labels)==1 && hidden.states.colors=="auto" && length(hidden.paths)==1 && is.null(hidden.paths)){
      attr(hidden.paths.seq, "cpal") <- seqHMM::colorpalette[[length(alphabet(hidden.paths.seq))]]
    }else if(all(isColor(hidden.states.colors))){
      if(length(hidden.states.colors)!=length(alphabet(hidden.paths.seq))){
        warning(paste0("Number of colors assigned to hidden.states.colors does not match the number of hidden states. \n
                       There are ", length(alphabet(hidden.paths.seq)), 
                       " hidden states but ", length(hidden.states.colors), " color(s)."))
      }
      attr(hidden.paths.seq, "cpal") <- rep(hidden.states.colors, length(alphabet(hidden.paths.seq)))[1:length(alphabet(hidden.paths.seq))]
    }else{
      stop(paste("Please provide a vector of colors for argument hidden.states.colors or use value \"auto\" to automatically determine gray scale colors."))
    }
    arguments <- c(arguments, list(hidden.paths.seq=hidden.paths.seq))
    # Sorting sequences according to multidimensional scaling score of hidden.paths
    if(!is.null(sortv)){
      if(length(sortv)==1 && sortv=="mds.hidden"){
        dist.hidden.paths <- suppressMessages(seqdist(hidden.paths.seq, method=dist.method, 
                                             sm="TRATE", with.missing=TRUE))
        sortv <- cmdscale(dist.hidden.paths, k=1)
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
      }
      # Sorting from start or end
      if(length(sortv)==1 && (sortv=="from.start" || sortv=="from.end")){
        end <- if (sortv == "from.end") {
          if(sort.channel>0 && sort.channel<=length(obs)){
            max(seqlength(obs[[sort.channel]]))
          }else if(sort.channel==0){
            if(plots=="both" || plots=="hidden.paths"){
              max(seqlength(hidden.paths.seq))
            }else{
              warning("Most probable paths are only computed automatically for argument plots=\"both\" or plots=\"hidden.paths\". Sequences were not sorted.")
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
            if(plots=="both" || plots=="hidden.paths"){
              max(seqlength(hidden.paths.seq))
            }else{
              warning("Most probable paths are only computed automatically for argument plots=\"both\" or plots=\"hidden.paths\". Sequences were not sorted.")
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
          orderv <- do.call(order, as.data.frame(hidden.paths.seq)[, end:beg])
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
            if(plots=="both" || plots=="hidden.paths"){
              max(seqlength(hidden.paths.seq))
            }else{
              warning("Most probable paths are only computed automatically for argument plots=\"both\" or plots=\"hidden.paths\". Sequences were not sorted.")
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
            if(plots=="both" || plots=="hidden.paths"){
              max(seqlength(hidden.paths.seq))
            }else{
              warning("Most probable paths are only computed automatically for argument plots=\"both\" or plots=\"hidden.paths\". Sequences were not sorted.")
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
          orderv <- do.call(order, as.data.frame(hidden.paths.seq)[, end:beg])
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
      if(length(sortv)==1 && sortv=="mds.hidden" && plots=="obs" && is.null(hidden.paths)){
        warning("Most probable paths are only computed automatically for argument plots=\"both\" or plots=\"hidden.paths\". Sequences were not sorted.")
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
    arguments <- c(arguments, list(hidden.paths=hidden.paths, plots=plots, type=type, 
                                   sortv=sortv, sort.channel=sort.channel, 
                                   with.missing=with.missing,
                                   title=title, title.n=title.n, cex.title=cex.title, title.pos=title.pos,
                                   withlegend=withlegend, ncol.legend=ncol.legend, 
                                   with.missing.legend=with.missing.legend,                         
                                   legend.prop=legend.prop, cex.legend=cex.legend,
                                   hidden.states.colors=hidden.states.colors, hidden.states.labels=hidden.states.labels,
                                   xaxis=xaxis, xlab=xlab, xtlab=xtlab, xlab.pos=xlab.pos,
                                   ylab=ylab, hidden.states.title=hidden.states.title, 
                                   ylab.pos=ylab.pos, 
                                   cex.lab=cex.lab, cex.axis=cex.axis, call=match.call()))
  }else{
    arguments <- c(arguments, list(hidden.paths=hidden.paths, plots=plots, type=type, 
                                   sortv=sortv, sort.channel=sort.channel, 
                                   with.missing=with.missing,
                                   title=title, title.n=title.n, cex.title=cex.title, title.pos=title.pos,
                                   withlegend=withlegend, ncol.legend=ncol.legend, 
                                   with.missing.legend=with.missing.legend,                         
                                   legend.prop=legend.prop, cex.legend=cex.legend,
                                   hidden.states.colors=hidden.states.colors, hidden.states.labels=hidden.states.labels,
                                   xaxis=xaxis, xlab=xlab, xtlab=xtlab, xlab.pos=xlab.pos,
                                   ylab=ylab, hidden.states.title=hidden.states.title, 
                                   ylab.pos=ylab.pos, 
                                   cex.lab=cex.lab, cex.axis=cex.axis, call=match.call()),
                   list(...))    
  }
  
  class(arguments) <- "ssp"
  arguments
  
}
