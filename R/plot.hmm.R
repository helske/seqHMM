
#' Plot hidden Markov models
#' 
#' Function \code{plot.hmm} plots a directed graph with pie charts of 
#' emission probabilities as vertices/nodes.
#' 
#' @import igraph
#' @export
#' @importFrom TraMineR seqlegend
#' @param x A hidden Markov model object of class hmm created with 
#'   \code{\link{build_hmm}} and \code{\link{fit_hmm}}. Multichannel 
#'   hmm objects are automatically transformed to single channel objects. 
#'   See function \code{\link{mc_to_sc}} for more information on the 
#'   transformation.
#' @param layout specifies the layout of the vertices (nodes). Accepts a 
#'   numerical matrix, a layout function, or either of \code{"horizontal"} (the 
#'   default) and \code{"vertical"}. Options \code{"horizontal"} and 
#'   \code{"vertical"} position vertices at the same horizontal or vertical 
#'   line. A two-column matrix can be used to give the x and y coordinates of 
#'   the vertices, or alternatively an igraph \code{\link[igraph]{layout}} 
#'   function can be called.
#' @param pie Are vertices plotted as pie charts of emission probabilities? 
#'   Defaulting to TRUE.
#' @param vertex.size The size of the vertex, given as a scalar or numerical 
#'   vector. The default value is 40.
#' @param vertex.label Labels for the vertices. Possible options include 
#'   \code{"initial.probs"}, \code{"names"}, \code{NA}, and a character or 
#'   numerical vector. The default \code{"initial.probs"} prints the initial 
#'   probabilities of the model and \code{"names"} prints the names of the 
#'   hidden states as labels. \code{NA} prints no labels.
#' @param vertex.label.dist The distance of the label of the vertex from its 
#'   center. The default value \code{"auto"} places the label outside the 
#'   vertex.
#' @param vertex.label.pos The position of the label of the vertex, relative to 
#'   the center of the vertices. A scalar or numerical vector giving the 
#'   position(s) as radians or one of \code{"bottom"} (pi/2 as radians), 
#'   \code{"top"} (-pi/2), \code{"left"} (pi), or \code{"right"} (0).
#' @param vertex.label.family,edge.label.family The font family to be used for
#'   vertex/edge labels. See argument \code{family} in \code{\link{par}} for
#'   more information.
#' @param loops Defines whether transitions to the same state are plotted.
#' @param edge.curved Defines whether to plot curved edges (arcs, arrows) 
#'   between the vertices. A logical or numerical vector or scalar. A numerical 
#'   value specifies the curvature of the edge. The default value \code{TRUE} 
#'   gives curvature of 0.5 to all edges. See \code{\link{igraph.plotting}} for 
#'   more information.
#' @param edge.label Labels for the edges. Possible options include 
#'   \code{"auto"}, \code{"NA"}, and a character or numerical vector. The 
#'   default \code{"auto"} prints the transition probabilities as edge labels. 
#'   \code{NA} prints no labels.
#' @param edge.width The width of the edges. The default \code{"auto"} plots the
#'   widths according to the transition probabilities between the hidden states.
#'   Other possibilities are a single value or a numerical vector giving the 
#'   widths.
#' @param cex.edge.width An expansion factor for the edge widths. Defaults to 1.
#' @param edge.arrow.size The size of the arrows in edges (constant). Defaults to 1.5.
#' @param label.signif Rounds labels of model parameters to the specified number
#'   of significant digits, 2 by default. Ignored for user-given labels.
#' @param label.scientific Defines if scientific notation is to be used to 
#'   describe small numbers. Defaults to \code{FALSE}, e.g. 0.0001 instead of 
#'   1e-04. Ignored for user-given labels.
#' @param label.max.length Maximum number of digits in labels of model 
#'   parameters. Ignored for user-given labels.
#' @param trim A scalar between 0 and 1 giving the highest probability of 
#'   transitions that are plotted as edges, defaults to 1e-15.
#' @param combine.slices A scalar between 0 and 1 giving the highest probability
#'   of emission probabilities that are combined into one state. The dafault 
#'   value is 0.05.
#' @param combined.slice.color The color of the slice including the smallest 
#'   emission probabilitis that user wants to combine (only if argument 
#'   \code{"combine.slices"} is greater than 0). The default color is white.
#' @param combined.slice.label The label for the combined states (when argument 
#'   \code{"combine.slices"} is greater than 0) to appear in the legend.
#' @param withlegend defines if and where the legend of the state colors is 
#'   plotted. Possible values include \code{"bottom"} (the default), 
#'   \code{"top"}, \code{"left"}, and \code{"right"}. \code{FALSE} omits the 
#'   legend.
#' @param ltext Optional description of the (combined) observed states to appear
#'   in the legend. A vector of character strings. See \code{\link{seqplot}} for
#'   more information.
#' @param legend.prop The proportion of the graphic area used for plotting the 
#'   legend. A scalar between 0 and 1, defaults to 0.5.
#' @param cex.legend Expansion factor for setting the size of the font for the 
#'   labels in the legend. The default value is 1. Values lesser than 1 will 
#'   reduce the size of the font, values greater than 1 will increase the size.
#' @param ncol.legend The number of columns for the legend. The default value 
#'   \code{"auto"} sets the number of columns automatically.
#' @param cpal Optional color palette for the (combinations of) observed states.
#'   The default value \code{"auto"} uses automatic color palette. Otherwise a 
#'   vector of length \code{x$n_symbols} is given, i.e. requires a color 
#'   specified for all (combinations of) observed states even if they are not 
#'   plotted (if the probability is less than combine.slices).
#' @param ... Other parameters passed on to \code{\link{plot.igraph}} such as 
#'   \code{vertex.color}, \code{vertex.label.cex}, \code{edge.lty}, 
#'   \code{margin}, or \code{main}.
#'   
#' @seealso \code{\link{build_hmm}} and \code{\link{fit_hmm}} for building and 
#'   fitting Hidden Markov models, \code{\link{mc_to_sc}} for transforming 
#'   multistate hmm objects to single channel objects, and 
#'   \code{\link{plot.igraph}} for the general plotting function of directed graphs.
#'   
#' @examples 
#' require(TraMineR)
#'   
#' data(biofam) 
#' 
#' biofam <- biofam[1:500,]
#' 
#' #' #' # Single-channel data
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
#' # Fitting HMM
#' HMM <- fit_hmm(bHMM)
#' 
#' # Plotting HMM
#' plot(HMM$model)
#' 
#' # Modifying edge curvature
#' plot(HMM$model, edge.curved = c(0, -0.6, 0.6, 0.3, -0.6, 0.3, 0))
#' 
#' #########################################
#' 
#' # Multichannel data
#'   
#' # Building one channel per type of event left, children or married 
#' bf <- as.matrix(biofam[, 10:25]) 
#' children <-  bf==4 | bf==5 | bf==6 
#' married <- bf == 2 | bf== 3 | bf==6 
#' left <- bf==1 | bf==3 | bf==5 | bf==6
#'   
#' children[children==TRUE] <- "Children" 
#' children[children==FALSE] <- "Childless"
#'   
#' married[married==TRUE] <- "Married" 
#' married[married==FALSE] <- "Single"
#'   
#' left[left==TRUE] <- "Left home" 
#' left[left==FALSE] <- "With parents"
#'   
#' # Building sequence objects 
#' child.seq <- seqdef(children) 
#' marr.seq <- seqdef(married) 
#' left.seq <- seqdef(left)
#'   
#' # Initial values for emission matrices 
#' B_child <- matrix(NA, nrow=4, ncol=2) 
#' B_child[1,] <- seqstatf(child.seq[,1:4])[,2]+0.1 
#' B_child[2,] <- seqstatf(child.seq[,5:8])[,2]+0.1 
#' B_child[3,] <- seqstatf(child.seq[,9:11])[,2]+0.1 
#' B_child[4,] <- seqstatf(child.seq[,12:15])[,2]+0.1 
#' B_child <- B_child/rowSums(B_child)
#'   
#' B_marr <- matrix(NA, nrow=4, ncol=2) 
#' B_marr[1,] <- seqstatf(marr.seq[,1:4])[,2]+0.1 
#' B_marr[2,] <- seqstatf(marr.seq[,5:8])[,2]+0.1 
#' B_marr[3,] <- seqstatf(marr.seq[,9:11])[,2]+0.1 
#' B_marr[4,] <- seqstatf(marr.seq[,12:15])[,2]+0.1 
#' B_marr <- B_marr/rowSums(B_marr)
#'   
#' B_left <- matrix(NA, nrow=4, ncol=2) 
#' B_left[1,] <- seqstatf(left.seq[,1:4])[,2]+0.1 
#' B_left[2,] <- seqstatf(left.seq[,5:8])[,2]+0.1 
#' B_left[3,] <- seqstatf(left.seq[,9:11])[,2]+0.1 
#' B_left[4,] <- seqstatf(left.seq[,12:15])[,2]+0.1 
#' B_left <- B_left/rowSums(B_left)
#'   
#' # Initial values for transition matrix 
#' A <- matrix(c(0.9,   0.06, 0.03, 0.01,
#'                 0,    0.9, 0.07, 0.03, 
#'                 0,      0,  0.9,  0.1, 
#'                 0,      0,    0,    1), 
#'             nrow=4, ncol=4, byrow=TRUE)
#'   
#' # Initial values for initial state probabilities 
#' initial_probs <- c(0.9, 0.07, 0.02, 0.01)
#'   
#' # Building hidden Markov model with initial parameter values 
#' bHMM <- build_hmm(observations=list(child.seq, marr.seq, left.seq), 
#'                  transition_matrix=A, 
#'                  emission_matrix=list(B_child, B_marr, B_left),
#'                  initial_probs=initial_probs)
#'   
#' # Fitting hidden Markov model 
#' HMM <- fit_hmm(bHMM)
#'   
#' # Plotting hmm object 
#' plot(HMM$model)
#' 
#' # Plotting HMM with
#' plot(HMM$model,
#'      # larger vertices
#'      vertex.size=50,
#'      # varying curvature of edges
#'      edge.curved=c(0,-0.7,0.6,0,-0.7,0),
#'      # legend with two columns and less space
#'      ncol.legend=2, legend.prop=0.4,
#'      # new label for combined slice
#'      combined.slice.label="States with probability < 0.05")
#' 
#' # Plotting HMM with given coordinates
#' plot(HMM$model,
#'      # layout given in 2x4 matrix
#'      # x coordinates in the first column
#'      # y coordinates in the second column
#'      layout=matrix(c(1,3,3,5,
#'                      0,-1,1,0), ncol=2),
#'      # larger vertices
#'      vertex.size=50,
#'      # straight edges
#'      edge.curved=FALSE,
#'      # thinner edges and arrows
#'      cex.edge.width=0.5, edge.arrow.size=1,
#'      # varying positions for vertex labels (initial probabilities)
#'      vertex.label.pos=c(pi, pi/2, -pi/2, 0),
#'        # different legend properties
#'        withlegend="top", legend.prop=0.3, cex.legend=1.1,
#'        # different axes
#'        xlim=c(0.5,5.5), ylim=c(-1.5,1.5), rescale=FALSE,
#'        # all states (not combining states with small probabilities)
#'        combine.slices=0,
#'        # legend with two columns
#'        ncol.legend=2)
#' 
#' # Plotting HMM with own color palette 
#' ## States with emission probability less than 0.2 removed
#' plot(HMM$model, cpal=1:10, combine.slices=0.2)
#'   
#' # Plotting HMM without pie graph and with a different layout 
#' require(igraph)
#' # Setting the seed for a random layout
#' set.seed(321)
#' plot(HMM$model,
#'      # Without pie graph
#'      pie=FALSE,
#'      # Fruchterman-Reingold layout
#'      layout=layout.fruchterman.reingold,
#'      # Straight edges and probabilities of moving to the same state
#'      edge.curved=FALSE, loops=TRUE,
#'      # Labels with three significant digits
#'      label.signif=3,
#'      # Fixed edge width
#'      edge.width=1,
#'      # Remove edges with probability less than 0.02
#'      trim=0.02,
#'      # Hidden state names as vertex labels
#'      vertex.label="names",
#'      # Labels insidde vertices
#'      vertex.label.dist=0)


plot.hmm <- function(x, layout="horizontal", pie=TRUE, 
                         vertex.size=40, vertex.label="initial.probs", 
                         vertex.label.dist="auto", vertex.label.pos="bottom",
                         vertex.label.family="sans",
                         loops=FALSE, edge.curved=TRUE, edge.label="auto", 
                         edge.width="auto", cex.edge.width=1, 
                         edge.arrow.size=1.5, edge.label.family="sans",
                         label.signif=2, label.scientific=FALSE, label.max.length=6,
                         trim=1e-15, 
                         combine.slices=0.05, combined.slice.color="white", 
                         combined.slice.label="others",
                         withlegend="bottom", ltext=NULL, legend.prop=0.5, 
                         cex.legend=1, ncol.legend="auto", cpal="auto", ...){
  
  
  # Saving and changing marginals
  oldPar <- par(no.readonly = TRUE)
  par(mar=c(0.5,0.5,0.5,0.5))
  on.exit(par(oldPar))
  on.exit(par(mfrow=c(1,1)))
  
  dots <- list(...)
  
  labelprint <- function(z, labs){
    if(labs==TRUE && (z > 0.001 || z==0)){
      labs <- FALSE
    }
    if(z < 10^-(label.max.length)){
      z <- prettyNum(signif(round(z, digits=label.max.length), digits=label.signif), scientific=labs)
    }else{
      z <- prettyNum(signif(z, digits=label.signif), scientific=labs)
    }
  }
  
  if(!is.matrix(layout) && !is.function(layout)){
    layout <- match.arg(layout, c("horizontal", "vertical"))
  }
  
  if(!is.numeric(vertex.label.pos)){
    vertex.label.pos <- match.arg(vertex.label.pos, c("bottom", "top", "left", "right"))
  }
  
  if(!is.logical(withlegend)){
    withlegend <- match.arg(withlegend, c("bottom", "top", "left", "right"))
  }
  
  if(x$n_channels>1){
    x <- mc_to_sc(x)
  }
  
  if(pie==FALSE && withlegend!=FALSE){
    withlegend <- FALSE
  }
  
  
  # Positions of vertex labels
  if(!is.numeric(vertex.label.pos)){
    if(vertex.label.pos=="bottom"){
      vertex.label.pos <- pi/2
    }else if(vertex.label.pos=="top"){
      vertex.label.pos <- -pi/2
    }else if(vertex.label.pos=="left"){
      vertex.label.pos <- pi
    }else{
      vertex.label.pos <- 0
    }
  }
  
  # Vertex labels
  if(length(vertex.label)==1){
    if(!is.na(vertex.label) && vertex.label!=FALSE){
      if(vertex.label=="initial.probs"){
        vertex.label <- sapply(x$initial_probs, labelprint, labs=label.scientific)
      }else if(vertex.label=="names"){
        vertex.label <- x$state_names
      }
    }
  }else if(length(vertex.label)<length(x$state_names)){
    warning("The length of the vector provided for the argument \"vertex.label\" is less than the number of hidden states. The vertor was repeated to archieve the correct length.")
    vertex.label <- rep(vertex.label, length.out=length(x$state_names))
  }else if(length(vertex.label)<length(x$state_names)){
    warning(paste("The length of the vector provided for the argument \"vertex.label\" is more than the number of number of hidden states. Only the first", length(x$state_names), "labels were used."))
    vertex.label <- vertex.label[1:length(x$state_names)]
  }
  
  # Vertex label distances
  if(is.character(vertex.label.dist)){
    match.arg(vertex.label.dist, c("auto"))
    vertex.label.dist <- vertex.size*0.4/10
  }else if(length(vertex.label.dist)>1 && length(vertex.label.dist)<x$n_states){
    warning("The length of the vector provided for the argument \"vertex.label.dist\" is less than the number of edges. The vector was repeated to archieve the correct length.")
    vertex.label.dist <- rep(vertex.label.dist, length.out=length(x$n_states))
  }else if(length(vertex.label.dist)>1 && length(vertex.label.dist)>x$n_states){
    warning(paste("The length of the vector provided for the argument \"vertex.label.dist\" is more than the number of edges. Only the first", length(x$n_states), "labels were used."))
    vertex.label.dist <- vertex.label.dist[1:length(x$n_states)]
  }
  
  
  # Trimming
  transM <- x$transition_matrix
  transM[transM<trim] <- 0
  
  # Adjacency matrix
  edges <- transM
  edges[edges>0] <- 1
  if(loops==FALSE){
    diag(edges) <- 0
  }
  
  # Vector of non-zero transition probabilities
  transitions <- transM
  if(loops==FALSE){
    diag(transitions) <- 0
  }
  transitions <- t(transitions)[t(transitions)>0]
  
  # Edge labels
  if(!is.na(edge.label) && edge.label!=FALSE){
    if(length(edge.label)==1 && is.character(edge.label)){
      match.arg(edge.label, c("auto"))
      edge.label <- sapply(transitions, labelprint, labs=label.scientific)
    }else if(length(edge.label)==1 && edge.label==TRUE){
      edge.label <- sapply(transitions, labelprint, labs=label.scientific)
    }else if(length(edge.label)>1 && length(edge.label)<length(transitions)){
      warning("The length of the vector provided for the argument \"edge.label\" is less than the number of edges. The vector was repeated to archieve the correct length.")
      edge.label <- rep(edge.label, length.out=length(transitions))
    }else if(length(edge.label)>1 && length(edge.label)>length(transitions)){
      warning(paste("The length of the vector provided for the argument \"edge.label\" is more than the number of edges. Only the first", length(transitions), "labels were used."))
      edge.label <- edge.label[1:length(transitions)]
    }
  }
  
  # Edge widths
  if(is.character(edge.width)){
    match.arg(edge.width, c("auto"))
    edge.width <- transitions*(7/max(transitions))*cex.edge.width
  }else if(length(edge.width)>1 && edge.width<length(transitions)){
    warning("The length of the vector provided for the argument \"edge.width\" is less than the number of edges. The vector was repeated to archieve the correct length.")
    edge.width <- rep(edge.width, length.out=length(transitions))
  }else if(length(edge.width)>1 && length(edge.width)>length(transitions)){
    warning(paste("The length of the vector provided for the argument \"edge.width\" is more than the number of edges. Only the first", length(transitions), "labels were used."))
    edge.width <- edge.width[1:length(transitions)]
  }

  # Defining the graph
  g1 <- graph.adjacency(edges, mode="directed")
  
  # Layout of the graph
  if(is.function(layout)){
    glayout <- layout(g1)
  }else if(is.matrix(layout)){
    glayout <- layout
  }else{
    if(layout=="horizontal"){
      glayout <- layout_on_grid(g1, width=x$n_states)
    }else if(layout=="vertical"){
      glayout <- layout_on_grid(g1, width=1)
    }
  }
  
  
  # Colors for the (combinations of) observed states
  if(length(cpal)==1 && cpal=="auto"){
    pie.colors <- attr(x$observations, "cpal")
  }else if(length(cpal)!=ncol(x$emiss)){
    warning("The length of the vector provided for argument cpal does not match the number of observed states. Automatic color palette was used.")
    pie.colors <- attr(x$observations, "cpal")
  }else if(!all(isColor(cpal))){
    stop(paste("Please provide a vector of colors for argument cpal or use value \"auto\" to use automatic color palette."))
  }else{
    pie.colors <- cpal
  }
  pie.colors.l <- pie.colors
  
  
  # Legend position and number of columns
  if(withlegend!=FALSE && pie==TRUE){
    if(!is.null(ltext)){
      if(length(ltext)!=x$n_symbols){
        warning("The length of the argument ltext does not match the number of observed states.")
      }
    }else{
      ltext <- x$symbol_names
    }
    if(withlegend=="bottom" || withlegend==TRUE){
      graphics::layout(matrix(1:2, nrow=2), heights=c(1-legend.prop, legend.prop))
    }else if(withlegend=="top"){
      graphics::layout(matrix(2:1, nrow=2), heights=c(legend.prop, 1-legend.prop))
    }else if(withlegend=="right"){
      graphics::layout(matrix(1:2, ncol=2), widths=c(1-legend.prop, legend.prop))
    }else{
      graphics::layout(matrix(2:1, ncol=2), widths=c(legend.prop, 1-legend.prop))
    }
    par(cex=1)
  }
  
  
  if(!is.matrix(layout) && !is.function(layout)){
    if(layout=="horizontal"){
      if(hasArg(rescale)){
        rescale <- dots$rescale
      }else{
        rescale=FALSE
      }
      if(hasArg(xlim)){
        xlim <- dots$xlim
      }else{
        if(rescale==TRUE){
          xlim <- c(-1,1)
        }else{
          xlim <- c(-0.1,ncol(transM)-1+0.1)
        }
      }
      if(hasArg(ylim)){
        ylim <- dots$ylim
      }else{
        if(rescale==TRUE){
          ylim <- c(-1,1)
        }else{
          ylim <- c(-0.5,0.5)
        }
      }
      dots[["xlim"]] <- NULL
      dots[["ylim"]] <- NULL
      dots[["rescale"]] <- NULL
    }else if(layout=="vertical"){
      if(hasArg(xlim)){
        xlim <- dots$xlim
      }else{
        if(rescale==TRUE){
          xlim <- c(-1,1)
        }else{
          xlim <- c(-0.5,0.5)
        }
      }
      if(hasArg(ylim)){
        ylim <- dots$ylim
      }else{
        if(rescale==TRUE){
          ylim <- c(-1,1)
        }else{
          ylim <- c(-0.1,ncol(transM)-1+0.1)
        }
      }
      dots[["xlim"]] <- NULL
      dots[["ylim"]] <- NULL
      dots[["rescale"]] <- NULL
    }
  }
  
  
  # Plotting graph
  if(pie==TRUE){
    pie.values <- lapply(seq_len(nrow(transM)), function(i) x$emission_matrix[i,])
    if(combine.slices>0){
      pie.colors.l <- NULL
      if(withlegend!=FALSE){
        lt <- NULL
        for(i in 1:x$n_states){
          cs.prob <- sum(pie.values[[i]][pie.values[[i]]<combine.slices])
          pie.values[[i]][pie.values[[i]]<combine.slices] <- 0
          pie.colors.l <- c(pie.colors.l,pie.colors[pie.values[[i]]>=combine.slices])
          lt <- c(lt, ltext[pie.values[[i]]>=combine.slices])
          pie.values[[i]] <- c(pie.values[[i]], cs.prob)
        }
        ltext <- c(unique(lt), combined.slice.label)
      }else{
        for(i in 1:x$n_states){
          cs.prob <- sum(pie.values[[i]][pie.values[[i]]<combine.slices])
          pie.values[[i]][pie.values[[i]]<combine.slices] <- 0
          pie.colors.l <- c(pie.colors.l,pie.colors[pie.values[[i]]>=combine.slices])
          pie.values[[i]] <- c(pie.values[[i]], cs.prob)
        }
      }
      pie.colors <- c(pie.colors, combined.slice.color)
      pie.colors.l <- c(unique(pie.colors.l), combined.slice.color)
      if(ncol.legend=="auto"){
        if(withlegend=="bottom" || withlegend==TRUE || withlegend=="top"){
          ncol.legend <- ceiling(length(pie.colors.l)/4)
        }else{
          ncol.legend <- 1
        }
      }
    }else{
      if(ncol.legend=="auto"){
        if(withlegend=="bottom" || withlegend==TRUE || withlegend=="top"){
          ncol.legend <- ceiling(ncol(x$emission_matrix)/4)
        }else{
          ncol.legend <- 1
        }
      }
    }
    
    if(!is.matrix(layout) && !is.function(layout) && (layout=="horizontal" || layout=="vertical")){
      do.call(plot.igraph2, c(list(g1, layout=glayout, 
                                  vertex.shape="pie", vertex.pie=pie.values,
                                  vertex.pie.color=list(pie.colors),
                                  vertex.size=vertex.size, 
                                  vertex.label=vertex.label, vertex.label.dist=vertex.label.dist, 
                                  vertex.label.degree=vertex.label.pos,
                                  vertex.label.family=vertex.label.family,
                                  edge.curved=edge.curved, edge.width=edge.width, 
                                  edge.label=edge.label, 
                                  edge.label.family=edge.label.family, 
                                  edge.arrow.size=edge.arrow.size,
                                  xlim=xlim, ylim=ylim, rescale=rescale), dots))
    }else{
      do.call(plot.igraph2, c(list(g1, layout=glayout, 
                                  vertex.shape="pie", vertex.pie=pie.values,
                                  vertex.pie.color=list(pie.colors),
                                  vertex.size=vertex.size, 
                                  vertex.label=vertex.label, vertex.label.dist=vertex.label.dist, 
                                  vertex.label.degree=vertex.label.pos,
                                  vertex.label.family=vertex.label.family,
                                  edge.curved=edge.curved, edge.width=edge.width, 
                                  edge.label=edge.label, 
                                  edge.label.family=edge.label.family,
                                  edge.arrow.size=edge.arrow.size), dots))
    }
  }else{
    if(!is.matrix(layout) && !is.function(layout) && (layout=="horizontal" || layout=="vertical")){
      do.call(plot.igraph2, c(list(g1, layout=glayout, 
                                  vertex.size=vertex.size, 
                                  vertex.label=vertex.label, vertex.label.dist=vertex.label.dist, 
                                  vertex.label.degree=vertex.label.pos,
                                  vertex.label.family=vertex.label.family,
                                  edge.curved=edge.curved, edge.width=edge.width, 
                                  edge.label=edge.label, 
                                  edge.label.family=edge.label.family, 
                                  xlim=xlim, ylim=ylim, rescale=rescale), dots))
    }else{
      do.call(plot.igraph2, c(list(g1, layout=glayout, 
                                  vertex.size=vertex.size, 
                                  vertex.label=vertex.label, vertex.label.dist=vertex.label.dist, 
                                  vertex.label.degree=vertex.label.pos,
                                  vertex.label.family=vertex.label.family,
                                  edge.curved=edge.curved, edge.width=edge.width, 
                                  edge.label=edge.label, 
                                  edge.label.family=edge.label.family), dots))
    }
  }  
  
  
  # Plotting legend
  if(withlegend!=FALSE && pie==TRUE){
    seqlegend(x$observations, cpal=pie.colors.l, ltext=ltext, 
              position="center", fontsize=cex.legend, ncol=ncol.legend,
              with.missing=FALSE)
    
  }
  
  
  par(oldPar)
  
}