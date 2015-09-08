#' Interactive Stacked Plots of Multichannel Sequences and/or Most Probable 
#' Paths for Mixture Hidden Markov Models
#' 
#' Function \code{mssplot} plots stacked sequence plots of observation sequences 
#' and/or most probable hidden state paths for each model of the \code{mhmm} 
#' object (model chosen according to the most probable path).
#' 
#' 
#' 
#' @export
#' 
#' @param x Mixture hidden Markov model object of class \code{mhmm}.
#'
#' @param ask If true and \code{which.plots} is NULL, \code{plot.mhmm} operates in interactive mode, via \code{\link{menu}}. Defaults to \code{FALSE}.
#' 
#' @param which.plots The number(s) of the requested model as an integer vector. The default \code{NULL} produces all plots.
#'   
#' @param mpp Output from \code{\link{hidden_paths}} function.
#'   
#' @param plots What to plot. One of \code{"obs"} for observations, \code{"mpp"}
#'   for most probable paths, or \code{"both"} for observations 
#'   and most probable paths (the default).
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
#' @param title A vector of titles for the graphics. The default is \code{NA}: if 
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
#' @param mpp.color A vector of colors assigned to hidden states (as ordered by 
#'   the \code{\link{hidden_paths}} function). The default value \code{"auto"} uses 
#'   the colors assigned to the stslist object created with \code{seqdef} if \code{mpp} 
#'   is given; otherwise colors from \code{\link{colorpalette}} are automatically used. 
#'   
#' @param mpp.labels Labels for the hidden states. The default value 
#'   \code{"auto"} uses the labels of the \code{mpp} argument if given; otherwise the number
#'   of the hidden state.
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
#' ## Building one channel per type of event left, children or married
#' bf <- as.matrix(biofam[, 10:25])
#' children <-  bf == 4 | bf == 5 | bf == 6
#' married <- bf == 2 | bf == 3 | bf == 6
#' left <- bf == 1 | bf == 3 | bf == 5 | bf == 6 | bf == 7
#' 
#' children[children == TRUE] <- "Children"
#' children[children == FALSE] <- "Childless"
#' # Divorced parents
#' div <- bf[(rowSums(bf == 7) > 0 & rowSums(bf == 5) > 0) | 
#'             (rowSums(bf == 7) > 0 & rowSums(bf == 6) > 0),]
#' children[rownames(bf) %in% rownames(div) & bf == 7] <- "Children"
#' 
#' married[married == TRUE] <- "Married"
#' married[married == FALSE] <- "Single"
#' married[bf == 7] <- "Divorced"
#' 
#' left[left == TRUE] <- "Left home"
#' left[left == FALSE] <- "With parents"
#' # Divorced living with parents (before divorce)
#' wp <- bf[(rowSums(bf == 7) > 0 & rowSums(bf == 2) > 0 & rowSums(bf == 3) == 0 &  
#'           rowSums(bf == 5) == 0 & rowSums(bf == 6) == 0) | 
#'          (rowSums(bf == 7) > 0 & rowSums(bf == 4) > 0 & rowSums(bf == 3) == 0 &  
#'          rowSums(bf == 5) == 0 & rowSums(bf == 6) == 0),]
#' left[rownames(bf) %in% rownames(wp) & bf == 7] <- "With parents"
#' 
#' ## Building sequence objects
#' child.seq <- seqdef(children, start = 15)
#' marr.seq <- seqdef(married, start = 15)
#' left.seq <- seqdef(left, start = 15)
#' 
#' ## Starting values for emission probabilities
#' 
#' # Cluster 1
#' alphabet(child.seq) # Checking for the order of observed states
#' B1_child <- matrix(c(0.99, 0.01, # High probability for childless
#'                      0.99, 0.01,
#'                      0.99, 0.01,
#'                      0.99, 0.01), nrow = 4, ncol = 2, byrow = TRUE)
#' 
#' alphabet(marr.seq)                      
#' B1_marr <- matrix(c(0.01, 0.01, 0.98, # High probability for single
#'                     0.01, 0.01, 0.98,
#'                     0.01, 0.98, 0.01, # High probability for married
#'                     0.98, 0.01, 0.01), # High probability for divorced
#'                     nrow = 4, ncol = 3, byrow = TRUE)                   
#' 
#' alphabet(left.seq)
#' B1_left <- matrix(c(0.01, 0.99, # High probability for living with parents
#'                     0.99, 0.01, # High probability for having left home
#'                     0.99, 0.01,
#'                     0.99, 0.01), nrow = 4, ncol = 2, byrow = TRUE)
#' 
#' # Cluster 2
#' B2_child <- matrix(c(0.99, 0.01, # High probability for childless
#'                      0.99, 0.01,
#'                      0.99, 0.01,
#'                      0.01, 0.99), nrow = 4, ncol = 2, byrow = TRUE)
#'                      
#' B2_marr <- matrix(c(0.01, 0.01, 0.98, # High probability for single
#'                     0.01, 0.01, 0.98,
#'                     0.01, 0.98, 0.01, # High probability for married
#'                     0.29, 0.7, 0.01),
#'                    nrow = 4, ncol = 3, byrow = TRUE)                   
#' 
#' B2_left <- matrix(c(0.01, 0.99, # High probability for living with parents
#'                     0.99, 0.01,
#'                     0.99, 0.01,
#'                     0.99, 0.01), nrow = 4, ncol = 2, byrow = TRUE) 
#' 
#' # Cluster 3
#' B3_child <- matrix(c(0.99, 0.01, # High probability for childless
#'                      0.99, 0.01,
#'                      0.01, 0.99,
#'                      0.99, 0.01,
#'                      0.01, 0.99,
#'                      0.01, 0.99), nrow = 6, ncol = 2, byrow = TRUE)
#' 
#' B3_marr <- matrix(c(0.01, 0.01, 0.98, # High probability for single
#'                     0.01, 0.01, 0.98,
#'                     0.01, 0.01, 0.98,
#'                     0.01, 0.98, 0.01,
#'                     0.01, 0.98, 0.01, # High probability for married
#'                     0.98, 0.01, 0.01), # High probability for divorced
#'                    nrow = 6, ncol = 3, byrow = TRUE)                   
#' 
#' B3_left <- matrix(c(0.01, 0.99, # High probability for living with parents
#'                     0.99, 0.01,
#'                     0.50, 0.50,
#'                     0.01, 0.99,
#'                     0.99, 0.01,
#'                     0.99, 0.01), nrow = 6, ncol = 2, byrow = TRUE) 
#' 
#' # Initial values for transition matrices
#' A1 <- matrix(c(0.8,   0.16, 0.03, 0.01,
#'                  0,    0.9, 0.07, 0.03, 
#'                  0,      0,  0.9,  0.1, 
#'                  0,      0,    0,    1), 
#'              nrow = 4, ncol = 4, byrow = TRUE)
#' 
#' A2 <- matrix(c(0.8, 0.10, 0.05,  0.03, 0.01, 0.01,
#'                  0,  0.7,  0.1,   0.1, 0.05, 0.05,
#'                  0,    0, 0.85,  0.01,  0.1, 0.04,
#'                  0,    0,    0,   0.9, 0.05, 0.05,
#'                  0,    0,    0,     0,  0.9,  0.1,
#'                  0,    0,    0,     0,    0,    1), 
#'              nrow = 6, ncol = 6, byrow = TRUE)
#' 
#' # Initial values for initial state probabilities 
#' initial_probs1 <- c(0.9, 0.07, 0.02, 0.01)
#' initial_probs2 <- c(0.9, 0.04, 0.03, 0.01, 0.01, 0.01)
#' 
#' # Birth cohort
#' biofam$cohort <- cut(biofam$birthyr, c(1908, 1935, 1945, 1957))
#' biofam$cohort <- factor(
#'   biofam$cohort, labels=c("1909-1935", "1936-1945", "1946-1957")
#' )
#' 
#' # Build mixture HMM
#' bMHMM <- buildMixHMM(
#'   observations = list(child.seq, marr.seq, left.seq),
#'   transition_matrix = list(A1,A1,A2),
#'   emission_matrix = list(list(B1_child, B1_marr, B1_left),
#'                         list(B2_child, B2_marr, B2_left), 
#'                         list(B3_child, B3_marr, B3_left)),
#'   initial_probs = list(initial_probs1, initial_probs1, initial_probs2),
#'   formula = ~ sex + cohort, data = biofam,
#'   cluster_names = c("Cluster 1", "Cluster 2", "Cluster 3"),
#'   channel_names = c("Parenthood", "Marriage", "Left home")
#'   )
#' 
#' 
#' \dontrun{
#' # Interactive plot
#'  mssplot(bMHMM)
#' }
#' # Plotting the first cluster only
#' mssplot(bMHMM, which.plots=1)
#'   
#' @seealso \code{\link{build_mhmm}} and \code{\link{fit_mhmm}} for building and 
#'   fitting Hidden Markov models, \code{\link{hidden_paths}} for 
#'   computing the most probable paths (Viterbi paths) of hidden states, and
#'   \code{\link{plot.mhmm}} for plotting parameters of \code{mhmm} objects.


mssplot <- function(x, ask = FALSE, which.plots = NULL, mpp=NULL,
                    plots="both", type="I", 
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
  
  # Checking for class of x
  if(!inherits(x, "mhmm")){
    stop("Your object x is not a mhmm object. Use build_mhmm and fit_mhmm to create one.")
  }
  
  oldPar <- par(no.readonly=TRUE)
  on.exit(par(oldPar))
  
  oldWarn <- options("warn")
  options(warn=1)
  on.exit(options(oldWarn), add=TRUE)
  
  allobs <- x$obs
  
  # ssp arguments (besides mhmm object and mpp)
  args <- as.list(match.call())[-(1:2)]
  if("ask" %in% names(args)){
    args <- args[-which(names(args)=="ask")]
  }
  if("which.plots" %in% names(args)){
    args <- args[-which(names(args)=="which.plots")]
  }
  if("mpp" %in% names(args)){
    args <- args[-which(names(args)=="mpp")]
  }
  if(!("title" %in% names(args))){
    titles <- x$cluster_names
  }else{
    if(length(title)!=x$number_of_clusters){
      warning("The length of the vector provided for the title argument does not match the number of clusters. Automatic titles were used instead.")
      titles <- x$cluster_names
    }else{
      titles <- args$title
    }
    args <- args[-which(names(args)=="title")]
  }
  if(length(ylab)==1 && ylab=="auto"){
    args$ylab <- x$channel_names
  }
  
  if(is.null(mpp)){
    mpp <- suppressWarnings(suppressMessages(hidden_paths(x)))
  }
  
  if(!("mpp.labels" %in% names(args))){
    mpp.labels <- NULL
    for(i in 1:x$number_of_clusters){
      mpp.labels <- c(mpp.labels, paste("State", 1:x$numberOfStates[i]))
    }
  }
  mpplabs <- list()
  k <- 0
  for(i in 1:x$number_of_clusters){
    mpplabs[[i]] <- mpp.labels[(k+1):(k+x$numberOfStates[i])]
    k <- k+x$numberOfStates[i]
  }
  
  if(!("mpp.color" %in% names(args))){
    mpp.color <- seqHMM::colorpalette[[length(alphabet(mpp$mpp))]]
  }
  mppcols <- list()
  k <- 0
  for(i in 1:x$number_of_clusters){
    mppcols[[i]] <- mpp.color[(k+1):(k+x$numberOfStates[i])]
    k <- k+x$numberOfStates[i]
  }
  
  mppm <- unique(mpp$cluster)
  mm <- NULL
  if(length(mppm)<x$number_of_clusters){
    mm <- which(!(x$cluster_names%in%mppm))
    warning(paste("When computing the most probable paths, no subjects were assigned to following clusters:", paste(x$cluster_names[mm], collapse=", ")))
  }
  
  if(!is.null(which.plots)){
    if(any(!is.numeric(which.plots)) || any(!(which.plots %in% 1:x$number_of_clusters))){
      stop(paste0("The which.plot argument only accepts numerical values between 1 and ", x$number_of_clusters, "."))
    }else if(any(which.plots %in% mm)){
      warning("You requested cluster(s) with no subjects. Plotting only relevant clusters.")
      which.plots <- setdiff(which.plots, mm)
    }
  }else if(!ask && is.null(which.plots)){
    which.plots <- 1:x$number_of_clusters
    # removing clusters with no subjects (according to mpp)
    which.plots <- setdiff(which.plots, mm)
  }
  
  
  
  if (ask && is.null(which.plots)) {
    tmenu <- 1:x$number_of_clusters
    tmenu <- setdiff(tmenu, mm)
    tmenunames <- x$cluster_names[tmenu]
    plot.new()
    repeat {
      pick <- menu(tmenunames, title = "\n Select cluster (or 0 to exit):\n")
      if(pick==0){
        return(invisible())
      }else{
        args$x <- lapply(allobs, function(y) y[mpp$cluster==x$cluster_names[[tmenu[pick]]],])
        args$mpp.labels <- mpplabs[[pick]]
        args$mpp <- suppressWarnings(suppressMessages(
          seqdef(mpp$mpp[mpp$cluster==x$cluster_names[[tmenu[pick]]],], 
                 labels=args$mpp.labels)))
        args$mpp.color <- mppcols[[pick]]
        args$title <- titles[tmenu[pick]]
        do.call(ssplotM,args=args)
      }
    }
  }else if (ask && !is.null(which.plots)) {
    tmenu <- which.plots
    tmenunames <- x$cluster_names[which.plots]
    plot.new()
    repeat {
      pick <- menu(tmenunames, title = "\n Select cluster (or 0 to exit):\n")
      if(pick==0){
        return(invisible())
      }else{
        args$x <- lapply(allobs, function(y) y[mpp$cluster==x$cluster_names[[tmenu[pick]]],])
        args$mpp.labels <- mpplabs[[pick]]
        args$mpp <- suppressWarnings(suppressMessages(
          seqdef(mpp$mpp[mpp$cluster==x$cluster_names[[tmenu[pick]]],], 
                 labels=args$mpp.labels)))
        args$mpp.color <- mppcols[[pick]]
        args$title <- titles[tmenu[pick]]
        do.call(ssplotM,args=args)
      }
    }
  }else{
    ask <- length(which.plots) > 1
    plot.new()
    for (i in which.plots) {
      args$x <- lapply(allobs, function(y) y[mpp$cluster==x$cluster_names[[i]],])
      args$mpp.labels <- mpplabs[[i]]
      args$mpp <- suppressWarnings(suppressMessages(
        seqdef(mpp$mpp[mpp$cluster==x$cluster_names[[i]],], labels=args$mpp.labels)))
      args$mpp.color <- mppcols[[i]]
      args$title <- titles[i]
      do.call(ssplotM,args=args)
      if (ask) {
        op <- par(ask = TRUE)
      }
    }
    # par(ask = FALSE)
  }
  invisible()
}