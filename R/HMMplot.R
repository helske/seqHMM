# Used in plot.mhmm through mHMMplotgrid

HMMplot <- function(x, layout="horizontal", pie=TRUE, 
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
                    cex.legend=1, ncol.legend="auto", cpal="auto", 
                    legend.pos="center", main = "auto", ...){
  
  
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
    if (!(layout %in% c("horizontal", "vertical"))) {
      stop("Argument layout only accepts numerical matrices, igraph layout functions, or strings \"horizontal\" and \"vertical\".")
    }
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
  transM <- x$transition_probs
  transM[transM<trim] <- 0
  
  # Adjacency matrix
  edges <- transM
  edges[edges>0] <- 1
  if(!loops){
    diag(edges) <- 0
  }
  
  # Vector of non-zero transition probabilities
  transitions <- transM
  if(!loops && length(transitions) > 1){
    diag(transitions)  <- 0
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
      glayout <- layout.grid(g1, width=x$n_states)
    }else if(layout=="vertical"){
      glayout <- layout.grid(g1, width=1)
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
  }
  
  
  if(!is.matrix(layout) && !is.function(layout)){
    if(layout=="horizontal"){
      if(hasArg(rescale)){
        rescale <- dots$rescale
      }else{
        rescale <- FALSE
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
    pie.values <- lapply(seq_len(nrow(transM)), function(i) x$emission_probs[i,])
    if(combine.slices>0 && !all(unlist(pie.values)[unlist(pie.values) > 0] > combine.slices)){
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
          ncol.legend <- ceiling(ncol(x$emission_probs)/4)
        }else{
          ncol.legend <- 1
        }
      }
    }
    
    if(!is.matrix(layout) && !is.function(layout) && 
         (layout=="horizontal" || layout=="vertical")){
      if(length(dots)>0){
        plotcall <- as.call(c(list(plot.igraph2, g1, layout=glayout, 
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
                         xlim=xlim, ylim=ylim, rescale=rescale, main = main), dots))
      }else{
        plotcall <- call("plot.igraph2", g1, layout=glayout, 
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
                         xlim=xlim, ylim=ylim, rescale=rescale, main = main)
      }
    }else{
      if(length(dots)>0){
        plotcall <- as.call(c(list(plot.igraph2, g1, layout=glayout, 
                         vertex.shape="pie", vertex.pie=pie.values,
                         vertex.pie.color=list(pie.colors),
                         vertex.size=vertex.size, 
                         vertex.label=vertex.label, vertex.label.dist=vertex.label.dist, 
                         vertex.label.degree=vertex.label.pos,
                         vertex.label.family=vertex.label.family,
                         edge.curved=edge.curved, edge.width=edge.width, 
                         edge.label=edge.label, 
                         edge.label.family=edge.label.family,
                         edge.arrow.size=edge.arrow.size, main = main), dots))
      }else{
        plotcall <- call("plot.igraph2", g1, layout=glayout, 
                         vertex.shape="pie", vertex.pie=pie.values,
                         vertex.pie.color=list(pie.colors),
                         vertex.size=vertex.size, 
                         vertex.label=vertex.label, vertex.label.dist=vertex.label.dist, 
                         vertex.label.degree=vertex.label.pos,
                         vertex.label.family=vertex.label.family,
                         edge.curved=edge.curved, edge.width=edge.width, 
                         edge.label=edge.label, 
                         edge.label.family=edge.label.family,
                         edge.arrow.size=edge.arrow.size, main = main)
      }
    }
  }else{
    if(!is.matrix(layout) && !is.function(layout) && 
         (layout=="horizontal" || layout=="vertical")){
      if(length(dots)>0){
        plotcall <- as.call(c(list(plot.igraph2, g1, layout=glayout, 
                         vertex.size=vertex.size, 
                         vertex.label=vertex.label, vertex.label.dist=vertex.label.dist, 
                         vertex.label.degree=vertex.label.pos,
                         vertex.label.family=vertex.label.family,
                         edge.curved=edge.curved, edge.width=edge.width, 
                         edge.label=edge.label, 
                         edge.label.family=edge.label.family, 
                         xlim=xlim, ylim=ylim, rescale=rescale, main = main), dots))
      }else{
        plotcall <- call("plot.igraph2", g1, layout=glayout, 
                         vertex.size=vertex.size, 
                         vertex.label=vertex.label, vertex.label.dist=vertex.label.dist, 
                         vertex.label.degree=vertex.label.pos,
                         vertex.label.family=vertex.label.family,
                         edge.curved=edge.curved, edge.width=edge.width, 
                         edge.label=edge.label, 
                         edge.label.family=edge.label.family, 
                         xlim=xlim, ylim=ylim, rescale=rescale, main = main)
      }
    }else{
      if(length(dots)>0){
        plotcall <- as.call(c(list(plot.igraph2, g1, layout=glayout, 
                         vertex.size=vertex.size, 
                         vertex.label=vertex.label, vertex.label.dist=vertex.label.dist, 
                         vertex.label.degree=vertex.label.pos,
                         vertex.label.family=vertex.label.family,
                         edge.curved=edge.curved, edge.width=edge.width, 
                         edge.label=edge.label, 
                         edge.label.family=edge.label.family, main = main), dots))
      }else{
        plotcall <- call("plot.igraph2", g1, layout=glayout, 
                         vertex.size=vertex.size, 
                         vertex.label=vertex.label, vertex.label.dist=vertex.label.dist, 
                         vertex.label.degree=vertex.label.pos,
                         vertex.label.family=vertex.label.family,
                         edge.curved=edge.curved, edge.width=edge.width, 
                         edge.label=edge.label, 
                         edge.label.family=edge.label.family, main = main)
      }
    }
  }  
  
  
  # Plotting legend
  if(withlegend!=FALSE && pie==TRUE){
    legendcall <- call("seqlegend", seqdata=x$observations, cpal=pie.colors.l, ltext=ltext, 
                       position=legend.pos, fontsize=cex.legend, ncol=ncol.legend,
                       with.missing=FALSE)
    
  }else{
    legendcall <- NULL
  }
  
  
  return(list(plotcall=plotcall, legendcall=legendcall))
  
  #graphics::layout(1)
}