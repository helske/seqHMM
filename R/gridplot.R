#' Plot Multidimensional Sequence Plots in a Grid
#' 
#' Function \code{gridplot} plots multiple \code{ssp} objects to a
#' grid.
#' 
#' 
#' 
#' @export

#'   
#' @param x A list of \code{ssp} objects.
#'   
#' @param nrow,ncol Optional arguments to arrange plots.
#'   
#' @param byrow Controls the order of plotting. Defaults to \code{FALSE}, i.e. plots
#'   are arranged columnwise.
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
#'   sequences in data contains missing states. With the value \code{TRUE} a
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
#' require(TraMineR)
#' data(biofam3c)
#' 
#' # Creating sequence objects
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
#' # Preparing plot for state distribution plots of observations for women
#' ssp_f <- ssp(
#'   list(child.seq[biofam3c$covariates$sex == "woman",], 
#'        marr.seq[biofam3c$covariates$sex == "woman",], 
#'        left.seq[biofam3c$covariates$sex == "woman",]),
#'   type = "d", plots = "obs", title = "Women", 
#'   ylab = c("Children", "Married", "Left home")
#'   )
#' 
#' # Preparing plot for state distribution plots of observations for men
#' ssp_m <- ssp(
#'   list(child.seq[biofam3c$covariates$sex == "man",], 
#'        marr.seq[biofam3c$covariates$sex == "man",], 
#'        left.seq[biofam3c$covariates$sex == "man",]), 
#'   type = "d", plots = "obs", title = "Men", 
#'   ylab = c("Children", "Married", "Left home")
#'   )
#' 
#' # Plotting state distribution plots of observations for women and men in two columns 
#' gridplot(list(ssp_f, ssp_m), ncol = 2, withlegend = FALSE)
#' 
#' # Preparing plots for women's state distributions
#' ssp_f2 <- ssp(
#'   list(marr.seq[biofam3c$covariates$sex == "woman",], 
#'        child.seq[biofam3c$covariates$sex == "woman",],
#'        left.seq[biofam3c$covariates$sex == "woman",]),
#'   type = "d", border = NA, withlegend = FALSE, 
#'   title = "State distributions for women", title.n = FALSE, xtlab = 15:30,
#'   ylab.pos = c(1, 2, 1), ylab = c("Married", "Children", "Left home")
#'   )
#' 
#' # The same plot with sequences instead of state distributions
#' ssp_f3 <- update(
#'   ssp_f2, type = "I", sortv="mds.obs", title = "Sequences for women"
#'   )
#' 
#' # State distributions with men's data
#' ssp_m2 <- update(
#'   ssp_f2, title = "State distributions for men", 
#'   x = list(marr.seq[biofam3c$covariates$sex == "man",], 
#'            child.seq[biofam3c$covariates$sex == "man",],
#'            left.seq[biofam3c$covariates$sex == "man",]),
#'   )
#' 
#' # Men's sequences
#' ssp_m3 <- update(
#'   ssp_m2, type = "I", sortv="mds.obs", title = "Sequences for women"
#'   )
#'   
#' # Plotting state distributions and index plots of observations 
#' # for women and men in two columns 
#' gridplot(
#'   list(ssp_f2, ssp_f3, ssp_m2, ssp_m3), ncol=3, byrow=TRUE, 
#'   withlegend="combined", legend.pos="right", col.prop=c(0.35,0.35,0.3)
#'   )
#'   
#' # The same with different positioning
#' gridplot(
#'   list(ssp_f2, ssp_f3, ssp_m2, ssp_m3), ncol = 2, nrow = 3, byrow=TRUE, 
#'   # defining the legend positions by the cell numbers
#'   legend.pos = 3:4
#'   )
#'                   
#'                   
#' @seealso \code{\link{ssp}} for defining the plot before using
#'   \code{gridplot}, and \code{\link{plot.ssp}} for plotting only one ssp object.

gridplot <- function(x, nrow=NA, ncol=NA, byrow=FALSE,
                     withlegend="auto", legend.pos="auto", 
                     legend.pos2="center", title.legend="auto",
                     ncol.legend="auto", 
                     with.missing.legend="auto",                       
                     row.prop="auto", col.prop="auto", cex.legend=1){
  grid.newpage()
  plot.new()
  opar <- par(no.readonly=TRUE)
  
  
  if(withlegend!=FALSE && withlegend!=TRUE){
    withlegend <- match.arg(withlegend, c("auto", "combined", "many"))
  }
  
  if(!is.numeric(ncol.legend) && ncol.legend!="auto"){
    warning("Argument ncol.legend only accepts values \"auto\" or a numerical vector.")
    ncol.legend <- "auto"
  }
  
  if(!is.numeric(row.prop) && row.prop!="auto"){
    warning("Argument row.prop only accepts values \"auto\" or a numerical vector.")
    row.prop <- "auto"
  }else if(is.numeric(row.prop) && all.equal(sum(row.prop),1)!=TRUE){
    warning("The elements of the vector provided for row.prop do not sum to 1. Argument row.prop was changed to \"auto\".")
    row.prop <- "auto"
  }
  
  if(!is.numeric(col.prop) && col.prop!="auto"){
    warning("Argument col.prop only accepts values \"auto\" or a numerical vector.")
    col.prop <- "auto"
  }else if(is.numeric(col.prop) && all.equal(sum(col.prop),1)!=TRUE){
    warning("The elements of the vector provided for col.prop do not sum to 1. Argument col.prop was changed to \"auto\".")
    col.prop <- "auto"
  }
  
  
  # Number of plots
  ngridplots <- length(x)
  
  # Checking for classes of x
  for(j in 1:ngridplots){
    if(!inherits(x[[j]], "ssp")){
      stop("At least one of your objects in x is not a ssp object. Use ssp to create one.")
    }
  }
  
  if(!is.na(withlegend) && withlegend!=FALSE && ngridplots > 1){
    nlegend <- x[[1]]$nplots
    for(i in 2:ngridplots){
      if(nlegend != x[[i]]$nplots){
        warning("The number of requested plots is not the same in all requested plots. Legends could not be printed.")
        withlegend <- FALSE
        break()
      }
      nlegend <- x[[i]]$nplots
    }  
  }
  
  if(x[[1]]$nchannels == 1 && (withlegend == "many" || 
                               withlegend == TRUE || 
                               withlegend == "auto")){
    withlegend <- "combined"
  }
  
  if(is.na(nrow) && is.na(ncol)){
    nrow <- ceiling(sqrt(ngridplots))
    ncol <- ceiling(ngridplots/nrow)
    rcfixed <- "none"
  }else if(is.na(nrow)){
    nrow <- ceiling(ngridplots/ncol)
    rcfixed <- "ncol"
  }else if(is.na(ncol)){
    ncol <- ceiling(ngridplots/nrow)
    rcfixed <- "nrow"
  }else{
    rcfixed <- "both"
  }
  
  
  
  emptycells <- nrow*ncol-ngridplots
  
  if((is.na(withlegend) && withlegend!=FALSE) && length(legend.pos)>1 && emptycells < ngridplots){
    warning("There were not enough empty cells for the legends. Legends were positioned automatically.")
    legend.pos <- "auto"
  }
  
  
  
  gridnrow <- nrow
  gridncol <- ncol
  
  emptycells <- gridnrow*gridncol-ngridplots
  
  # Legend titles
  if(!is.na(withlegend) && withlegend!=FALSE){
    if(!is.na(title.legend) && title.legend!=FALSE && !is.null(title.legend)){
      if(length(title.legend)>1 || (length(title.legend)==1 && title.legend!="auto")){
        if(length(title.legend)!=x[[1]]$nplots){
          warning("The length of the vector provided for title.legend does not match the number of legends. Argument title.legend was set to \"auto\".")
          title.legend=="auto"
        }
      }
      if(length(title.legend==1 && title.legend=="auto")){
        title.legend <- x[[1]]$ylab
      }
    }
  }
  
  # Legend positions
  if(!is.na(withlegend) && withlegend!=FALSE){
    if(length(legend.pos)>1){
      if(byrow==TRUE){
        if(max(legend.pos)-min(legend.pos)+1>length(legend.pos)){
          legendp <- "right"
        }else{
          legendp <- "bottom"
        }
        # byrow=FALSE
      }else{
        if(max(legend.pos)-min(legend.pos)+1>length(legend.pos)){
          legendp <- "bottom"
        }else{
          legendp <- "right"
        }
      }
    }      
    
    # Legends at bottom
    if(length(legend.pos)==1 && (legend.pos=="auto" || legend.pos=="bottom")){
      legendp <- "bottom"
      if(byrow==FALSE){
        if(emptycells==0){
          if(rcfixed=="both"){
            warning(paste0("Legend does not fit to requested grid with ", nrow, " nrow and ", ncol, " columns. A row was added."))
            gridnrow <- nrow+1
            emptycells <- gridnrow*gridncol-ngridplots
          }else if(rcfixed=="ncol" || rcfixed=="none"){
            gridnrow <- nrow+1
            emptycells <- gridnrow*gridncol-ngridplots
          }else{
            gridncol <- ncol+1
            emptycells <- gridnrow*gridncol-ngridplots
          }
        }
        if(emptycells<=gridncol){
          legend.pos <- c(c(gridncol:1)*gridnrow)[emptycells:1]
        }else{
          legend.pos <- c(gridncol:1)*gridnrow
        }
      }else{
        if(emptycells==0){
          if(rcfixed=="both"){
            warning(paste0("Legend does not fit to requested grid with ", nrow, " nrow and ", ncol, " columns. A row was added."))
            gridnrow <- nrow+1
            emptycells <- gridnrow*gridncol-ngridplots
          }else if(rcfixed=="ncol" || rcfixed=="none"){
            gridnrow <- nrow+1
            emptycells <- gridnrow*gridncol-ngridplots
          }else{
            gridncol <- ncol+1
            emptycells <- gridnrow*gridncol-ngridplots
          }
        }
        if(emptycells<=gridncol){
          legend.pos <- c((gridnrow*gridncol-emptycells+1):(gridnrow*gridncol))
        }else{
          legend.pos <- c((gridnrow*gridncol-floor(emptycells/gridncol)*gridncol+1):(gridnrow*gridncol))
        }
      }        
      if(length(ncol.legend)==1 && ncol.legend=="auto" && withlegend!=TRUE && 
         withlegend=="combined"){
        ncol.legend <- x[[1]]$nplots
      }
      # Legend at right
    }else if(length(legend.pos)==1 && legend.pos=="right"){
      legendp <- "right"
      if(byrow==FALSE){
        if(emptycells==0){
          if(rcfixed=="both"){
            warning(paste0("Legend does not fit to requested grid with ", nrow, " nrow and ", ncol, " columns. A column was added."))
            gridncol <- ncol+1
            emptycells <- gridnrow*gridncol-ngridplots
          }else if(rcfixed=="nrow" || rcfixed=="none"){
            gridncol <- ncol+1
            emptycells <- gridnrow*gridncol-ngridplots
          }else{
            gridnrow <- nrow+1
            emptycells <- gridnrow*gridncol-ngridplots
          }
        }
        if(emptycells<=gridnrow){
          legend.pos <- c((gridnrow*gridncol-emptycells+1):(gridnrow*gridncol))
        }else{
          legend.pos <- c((gridnrow*gridncol-floor(emptycells/gridnrow)*gridnrow+1):(gridnrow*gridncol))
        }
        # byrow=TRUE
      }else{
        if(emptycells==0){
          if(rcfixed=="both"){
            warning(paste0("Legend does not fit to requested grid with ", nrow, " nrow and ", ncol, " columns. A column was added."))
            gridncol <- ncol+1
            emptycells <- gridnrow*gridncol-ngridplots
          }else if(rcfixed=="nrow" || rcfixed=="none"){
            gridncol <- ncol+1
            emptycells <- gridnrow*gridncol-ngridplots
          }else{
            gridnrow <- nrow+1
            emptycells <- gridnrow*gridncol-ngridplots
          }
        }
        if(emptycells<=gridnrow){
          legend.pos <- c(c(gridnrow:1)*gridncol)[emptycells:1]
        }else{
          legend.pos <- c(gridnrow:1)*gridncol
        }
      }  
      if(length(ncol.legend)==1 && ncol.legend=="auto" && withlegend!=TRUE && 
         withlegend=="combined"){
        ncol.legend <- 1
      }
    }
    
    if(withlegend==TRUE || withlegend=="auto" || withlegend=="many"){
      if(length(ncol.legend)==1 && ncol.legend=="auto"){
        ncol.legend <- rep(1, x[[1]]$nplots)
        legend.nrow <- sapply(lapply(x[[1]]$obs, "alphabet"), "length")
      }else if(length(ncol.legend)==1 && x[[1]]$nplots>1){
        legend.nrow <- ceiling(sapply(lapply(x[[1]]$obs, "alphabet"), "length")/ncol.legend)
        ncol.legend <- rep(ncol.legend, x[[1]]$nplots)
      }else if(length(ncol.legend)<x[[1]]$nplots){
        warning(paste0("The length of ncol.legend does not match the number of requested plots. The last were arranged in 1 column."))
        ncol.legend <- c(ncol.legend, rep(1,(x[[1]]$nplots-length(ncol.legend))))
        legend.nrow <- ceiling(sapply(lapply(x[[1]]$obs, "alphabet"), "length")/ncol.legend)
      }else if(length(ncol.legend)>x[[1]]$nplots){
        warning(paste0("The length of ncol.legend does not match the number of requested plots. Only the first ", x[[1]]$nplots, " arguments of \"ncol.legend\" were used."))
        legend.nrow <- ceiling(sapply(lapply(x[[1]]$obs, "alphabet"), "length")/ncol.legend[1:x[[1]]$nplots])
      }else{
        legend.nrow <- ceiling(x$n_states/ncol.legend)
      }
      
      if(!is.na(title.legend) && title.legend!=FALSE && !is.null(title.legend)){
        legend.nrow <- legend.nrow+1
      }
    }
  }
  
  
  
  
  # Cells' proportions
  if(!is.numeric(row.prop) && row.prop=="auto"){
    row.prop <- rep(1/gridnrow, gridnrow)
  }
  if(!is.numeric(col.prop) && col.prop=="auto"){
    col.prop <- rep(1/gridncol, gridncol)
  }
  if(length(row.prop)!=gridnrow){
    warning("The length of the vector provided for row.prop does not match the number of nrow in the plot. Argument row.prop was changed to \"auto\".")
    row.prop <- rep(1/gridnrow, gridnrow)
  }
  if(length(col.prop)!=gridncol){
    warning("The length of the vector provided for col.prop does not match the number of columns in the plot. Argument col.prop was changed to \"auto\".")
    col.prop <- rep(1/gridncol, gridncol)
  }
  
  
  if(byrow==FALSE){
    plotnrow <- rep(c(1:gridnrow), times=gridncol)
    plotncol <- rep(c(1:gridncol), each=gridnrow)
    plotgrid <- cbind(plotnrow,plotncol)
    if(!is.na(withlegend) && withlegend!=FALSE){
      lpos <- matrix(c(1:nrow(plotgrid)), byrow=FALSE, nrow=gridnrow)
      legendplace <- matrix(c(lpos[,] %in% legend.pos), byrow=FALSE, nrow=gridnrow)
      plotplace <- matrix(c(!(lpos[,] %in% legend.pos)), byrow=FALSE, nrow=gridnrow)
      plotcells <- plotgrid[c(plotplace),,drop=FALSE][1:ngridplots,]
      plotlegend <- plotgrid[c(legendplace),,drop=FALSE]
      if(legendp=="bottom"){
        plotlegend <- plotlegend[order(plotlegend[,2]),,drop=FALSE]
      }
      if(length(ncol.legend)==1 && ncol.legend=="auto"){
        if(withlegend=="combined"){
          if(length(unique(plotlegend[,1]))>=length(unique(plotlegend[,2]))){
            ncol.legend <- 1
          }else{
            ncol.legend <- x[[1]]$nplots
          }
        }else{
          rep(1,x[[1]]$nplots)
        }
      }
    }else{
      plotcells <- plotgrid[1:ngridplots,,drop=FALSE]
    }
  }else{
    plotnrow <- rep(c(1:gridnrow), each=gridncol)
    plotncol <- rep(c(1:gridncol), times=gridnrow)
    plotgrid <- cbind(plotnrow,plotncol)
    if(!is.na(withlegend) && withlegend!=FALSE){
      lpos <- matrix(c(1:nrow(plotgrid)), byrow=TRUE, nrow=gridnrow)
      legendplace <- matrix(c(lpos[,] %in% legend.pos), byrow=FALSE, nrow=gridnrow)
      plotplace <- matrix(c(!(lpos[,] %in% legend.pos)), byrow=FALSE, nrow=gridnrow)
      plotcells <- plotgrid[c(t(plotplace)),,drop=FALSE][1:ngridplots,]
      plotlegend <- plotgrid[c(t(legendplace)),,drop=FALSE]
      if(legendp=="right"){
        plotlegend <- plotlegend[order(plotlegend[,2]),,drop=FALSE]
      }
      if(length(ncol.legend)==1 && ncol.legend=="auto"){
        if(withlegend=="combined"){
          if(length(unique(plotlegend[,1]))>=length(unique(plotlegend[,2]))){
            ncol.legend <- 1
          }else{
            ncol.legend <- x[[1]]$nplots
          }
        }else{
          rep(1,x[[1]]$nplots)
        }
      }
    }else{
      plotcells <- plotgrid[1:ngridplots,,drop=FALSE]
    }
  }
  
  
  multitop.vp <- viewport(layout=
                            grid.layout(gridnrow,gridncol,
                                        widths = do.call(unit,
                                                         args=list(col.prop, 
                                                                   rep("npc", 
                                                                       length(col.prop)))), 
                                        heights = do.call(unit, 
                                                          args=list(row.prop, 
                                                                    rep("npc", 
                                                                        length(row.prop))))), 
                          width=unit(1, "npc"))
  for(i in 1:ngridplots){
    assign(paste0("vpplot",i), viewport(layout.pos.row=plotcells[i,1], 
                                        layout.pos.col=plotcells[i,2], 
                                        name=paste0("vpplot",i)))
  }
  if(!is.na(withlegend) && withlegend!=FALSE){
    #     if(withlegend==TRUE || withlegend!="combined"){
    #       for(i in 1:x[[1]]$nplots){
    #         assign(paste0("vplegend",i), viewport(layout.pos.row=plotlegend[i,1], 
    #                                               layout.pos.col=plotlegend[i,2], 
    #                                               name=paste0("vplegend",i)))
    #       }
    #       vpall <- vpTree(multitop.vp, do.call(vpList, 
    #                                            args=mget(c(paste0("vpplot",1:ngridplots), 
    #                                                        paste0("vplegend",1:x[[1]]$nplots)))))
    #     }else{
    assign("vplegend", viewport(layout.pos.row=unique(plotlegend[,1]), 
                                layout.pos.col=unique(plotlegend[,2]), 
                                name="vplegend"))
    vpall <- vpTree(multitop.vp, do.call(vpList, 
                                         args=mget(c(paste0("vpplot",1:ngridplots), 
                                                     "vplegend"))))
    #     }
  }else{
    vpall <- vpTree(multitop.vp, do.call(vpList, args=mget(paste0("vpplot",1:ngridplots))))
  }
  
  
  pushViewport(vpall)
  
  upViewport()
  
  # Plots
  for(p in 1:ngridplots){
    downViewport(paste0("vpplot",p))
    do.call(SSPlotter,args=x[[p]])
    popViewport()
    upViewport()
  }
  
  for(i in 1:ngridplots){
    if(x[[i]]$nchannels == 1){
      x[[i]]$obs <- list(x[[i]]$obs)
    }
  }
  
  # Legends
  if(!is.na(withlegend) && withlegend!=FALSE){
    if(x[[1]]$plots=="both" || x[[1]]$plots=="mpp"){
      maxhs <- 0
      maxhsplot <- 1
      for(i in 1:ngridplots){
        if(length(alphabet(x[[i]]$mpp.seq))>maxhs){
          maxhs <- length(alphabet(x[[i]]$mpp.seq))
          maxhsplot <- i
        }
      }
    }
    ltext <- NULL
    cpal <- NULL
    if(x[[1]]$plots=="both" || x[[1]]$plots=="obs"){
      ltexts <- rep(list(NULL), x[[1]]$nchannels)
      cpals <- rep(list(NULL), x[[1]]$nchannels)
      anymissing <- FALSE
      for(i in 1:ngridplots){
        for(j in 1:x[[1]]$nchannels){
          ltexts[[j]] <- unique(c(c(ltexts[[j]], attr(x[[i]]$obs[[j]], "labels"))))
          cpals[[j]] <- unique(c(cpals[[j]], attr(x[[i]]$obs[[j]], "cpal")))
          if(any(x[[i]]$obs[[j]]=="*")){
            anymissing <- TRUE
          }
        }
      }
      for(j in 1:x[[1]]$nchannels){
        ltext <- c(ltext, unlist(ltexts[[j]]))
        cpal <- c(cpal, unlist(cpals[[j]]))
      }
    }
    if(x[[1]]$plots=="both" || x[[1]]$plots=="mpp"){
      ltext <- c(ltext, attr(x[[maxhsplot]]$mpp.seq, "labels"))
      cpal <- c(cpal, attr(x[[maxhsplot]]$mpp.seq, "cpal")) 
    }
    
    # Separate legends
    if(withlegend==TRUE || withlegend=="auto" || withlegend=="many"){
      downViewport("vplegend")
      # Vertical legends
      if(legendp=="right"){
        pushViewport(viewport(layout=
                                grid.layout(nrow=x[[1]]$nplots, ncol=1,
                                            heights=unit((legend.nrow/sum(legend.nrow)), 
                                                         "npc")), 
                              width=unit(0.95, "npc")))
        # Legends for channels
        if(x[[1]]$plots=="both" || x[[1]]$plots=="obs"){
          for(i in 1:x[[1]]$nchannels){
            pushViewport(viewport(layout.pos.col=1, layout.pos.row=i))
            par(plt=gridPLT(), new=TRUE)
            seqlegend(x[[1]]$obs[[i]], fontsize=cex.legend, position=legend.pos2, 
                      cpal=cpals[[i]], ltext=ltexts[[i]],
                      ncol=ncol.legend[i], with.missing=x[[1]]$with.missing.legend,
                      title=title.legend[i])
            popViewport()
          }
        }      
        # Legends for most probable paths
        if(x[[1]]$plots=="both" || x[[1]]$plots=="mpp"){
          pushViewport(viewport(layout.pos.col=1, layout.pos.row=x[[1]]$nplots))
          par(plt=gridPLT(), new=TRUE)
          seqlegend(x[[maxhsplot]]$mpp.seq, fontsize=cex.legend, 
                    position=legend.pos2, ncol=ncol.legend[length(ncol.legend)], 
                    with.missing=x[[maxhsplot]]$with.missing.legend,
                    title=title.legend[length(title.legend)])
          popViewport()
        }
        popViewport()
        # Horizontal legends
      }else{
        pushViewport(viewport(layout=
                                grid.layout(ncol=x[[1]]$nplots, nrow=1,
                                            widths=unit((legend.nrow/sum(legend.nrow)), 
                                                        "npc")), 
                              width=unit(0.95, "npc")))
        # Legends for channels
        if(x[[1]]$plots=="both" || x[[1]]$plots=="obs"){
          for(i in 1:x[[1]]$nchannels){
            pushViewport(viewport(layout.pos.col=i, layout.pos.row=1))
            par(plt=gridPLT(), new=TRUE)
            seqlegend(x[[1]]$obs[[i]], fontsize=cex.legend, position=legend.pos2, 
                      cpal=cpals[[i]], ltext=ltexts[[i]],
                      ncol=ncol.legend[i], with.missing=x[[1]]$with.missing.legend,
                      title=title.legend[i])
            popViewport()
          }
        }      
        # Legends for most probable paths
        if(x[[1]]$plots=="both" || x[[1]]$plots=="mpp"){
          pushViewport(viewport(layout.pos.col=x[[1]]$nplots, layout.pos.row=1))
          par(plt=gridPLT(), new=TRUE)
          seqlegend(x[[maxhsplot]]$mpp.seq, fontsize=cex.legend, 
                    position=legend.pos2, ncol=ncol.legend[length(ncol.legend)], 
                    with.missing=x[[maxhsplot]]$with.missing.legend,
                    title=title.legend[length(title.legend)])
          popViewport()
        }
        popViewport()
      }
      
      #       if(x[[1]]$plots=="both" || x[[1]]$plots=="obs"){
      #         for(i in 1:x[[1]]$nchannels){
      # #           downViewport(paste0("vplegend",i))
      #           par(plt=gridPLT(), new=TRUE)
      #           seqlegend(x[[1]]$obs[[i]], fontsize=cex.legend, position=legend.pos2, 
      #                     cpal=cpals[[i]], ltext=ltexts[[i]],
      #                     ncol=ncol.legend[i], with.missing=x[[1]]$with.missing.legend,
      #                     title=title.legend[i])
      #           popViewport()
      #         }
      #       }
      #       # Legends for most probable paths
      #       if(x[[1]]$plots=="both" || x[[1]]$plots=="mpp"){
      #         downViewport(paste0("vplegend",x[[1]]$nplots))
      #         par(plt=gridPLT(), new=TRUE)
      #         seqlegend(x[[maxhsplot]]$mpp.seq, fontsize=cex.legend, 
      #                   position=legend.pos2, ncol=ncol.legend[length(ncol.legend)], 
      #                   with.missing=x[[maxhsplot]]$with.missing.legend,
      #                   title=title.legend[length(title.legend)])
      #         popViewport()
      #       }
      #       popViewport()
      
      # Combined legends
    }else if(withlegend=="combined"){
      downViewport("vplegend")
      par(plt=gridPLT(), new=TRUE)
      pushViewport(viewport(width=unit(0.9, "npc")))
      
      seqlegend(x[[1]]$obs[[1]], fontsize=cex.legend, position=legend.pos2, 
                ncol=ncol.legend, cpal=cpal, ltext=ltext,
                with.missing=anymissing, 
                missing.color=attr(x[[1]]$obs[[1]],"missing.color"))
      
      popViewport()
    }
  }
  
  par(opar)
}

