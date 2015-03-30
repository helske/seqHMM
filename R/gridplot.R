#' Plot Multidimensional Sequence Plots in a Grid
#' 
#' Function \code{gridplot} plots multiple \code{ssp} objects to a
#' grid.
#' 
#' 
#' 
#' @export
#' @import gridBase
#' @import grid
#'   
#' @param x A list of \code{ssp} objects.
#'   
#' @param rows,cols Optional arguments to arrange plots.
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
#'   cell(s). One of \code{"bottomright"}, \code{"bottom"}, \code{"bottomleft"}, \code{"left"}, 
#'   \code{"topleft"}, \code{"top"}, \code{"topright"}, \code{"right"} and \code{"center"}.
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
#' data(biofam)
#' biofam <- biofam[1:500,]
#' 
#' ## Building one channel per type of event left, children or married
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
#' # Preparing plot for state distribution plots of observations for women
#' ssp_f <- ssp(list(child.seq[biofam$sex=="woman",], 
#' marr.seq[biofam$sex=="woman",], left.seq[biofam$sex=="woman",]),
#'                    type="d", plots="obs", title="Women", 
#'                    ylab=c("Children", "Married", "Left home"))
#' 
#' # Preparing plot for state distribution plots of observations for men
#' ssp_m <- ssp(list(child.seq[biofam$sex=="man",], 
#' marr.seq[biofam$sex=="man",], left.seq[biofam$sex=="man",]), 
#'                    type="d", plots="obs", title="Men", 
#'                    ylab=c("Children", "Married", "Left home"))
#' 
#' # Plotting state distribution plots of observations for women and men in two columns 
#' gridplot(list(ssp_f, ssp_m), cols=2, withlegend=FALSE)
#' 
#' 
#' # Preparing plots for state distributios and index plots of observations for women
#' ssp_f2 <- ssp(list(child.seq[biofam$sex=="woman",], 
#' marr.seq[biofam$sex=="woman",], left.seq[biofam$sex=="woman",]), 
#'                     type="d", plots="obs", title="Women", 
#'                     ylab=c("Children", "Married", "Left home"), withlegend=FALSE)
#' ssp_f3 <- ssp(list(child.seq[biofam$sex=="woman",], 
#' marr.seq[biofam$sex=="woman",], left.seq[biofam$sex=="woman",]), 
#'                     type="I", plots="obs", title="Women", 
#'                     ylab=c("Children", "Married", "Left home"), withlegend=FALSE)
#' 
#' # Preparing plots for state distributios and index plots of observations for men
#' ssp_m2 <- ssp(list(child.seq[biofam$sex=="man",], 
#' marr.seq[biofam$sex=="man",], left.seq[biofam$sex=="man",]), 
#'                     type="d", plots="obs", title="Men", 
#'                     ylab=c("Children", "Married", "Left home"), withlegend=FALSE)
#' ssp_m3 <- ssp(list(child.seq[biofam$sex=="man",], 
#' marr.seq[biofam$sex=="man",], left.seq[biofam$sex=="man",]),
#'                     type="I", plots="obs", title="Men", 
#'                     ylab=c("Children", "Married", "Left home"), withlegend=FALSE)
#' 
#' # Plotting state distributions and index plots of observations for women and men in two columns 
#' gridplot(list(ssp_f2, ssp_f3, ssp_m2, ssp_m3), cols=2, byrow=TRUE, 
#'                   withlegend="combined", legend.pos="bottom", row.prop=c(0.4,0.4,0.2))
#'                   
#'                   
#' @seealso \code{\link{ssp}} for defining the plot before using
#'   \code{gridplot}, and \code{\link{plot.ssp}} for plotting only one ssp object.

gridplot <- function(x, rows=NA, cols=NA, byrow=FALSE,
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
        warning("The number of requested plots is not the same in all for all requested plots. Legends could not be printed.")
        withlegend <- FALSE
        break()
      }
      nlegend <- x[[i]]$nplots
    }  
  }
  
  
  if(is.na(rows) && is.na(cols)){
    rows <- ceiling(sqrt(ngridplots))
    cols <- ceiling(ngridplots/rows)
    rcfixed <- "none"
  }else if(is.na(rows)){
    rows <- ceiling(ngridplots/cols)
    rcfixed <- "cols"
  }else if(is.na(cols)){
    cols <- ceiling(ngridplots/rows)
    rcfixed <- "rows"
  }else{
    rcfixed <- "both"
  }
  
  
  
  emptycells <- rows*cols-ngridplots
  
  if((is.na(withlegend) && withlegend!=FALSE) && length(legend.pos)>1 && emptycells < ngridplots){
    warning("There were not enough empty cells for the legends. Legends were positioned automatically.")
    legend.pos <- "auto"
  }
  
  
  
  gridrows <- rows
  gridcols <- cols
  
  emptycells <- gridrows*gridcols-ngridplots
  
  if(!is.numeric(row.prop) && length(row.prop)!=gridrows){
  }
  
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
    # # Combined legend box
    # if(withlegend!=TRUE && withlegend=="combined"){
    # Non-adjacent cells for combined legend box
#     if(withlegend!=TRUE && withlegend=="combined"){
#       if(length(legend.pos)>1 && (max(legend.pos)-min(legend.pos)+1)>length(legend.pos)){
#         
#         warning("The legend positions (cells) must be in one row/column. Argument legend.pos was set to \"auto\".")
#         legendp <- "bottom"
#         if(byrow==FALSE){
#           if(emptycells<=gridcols){
#             legend.pos <- c(c(gridcols:1)*gridrows)[emptycells:1]
#           }else{
#             legend.pos <- c(c(gridcols:1)*gridrows)
#           }
#         }else{
#           legend.pos <- c((gridrows*gridcols-x[[1]]$nplots+1):(gridrows*gridcols))
#         }        
#         if(length(ncol.legend)==1 && ncol.legend=="auto"){
#           ncol.legend <- x[[1]]$nplots
#         }
#       }
# #     }else 
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
            warning(paste0("Legend does not fit to requested grid with ", rows, " rows and ", cols, " columns. A row was added."))
            gridrows <- rows+1
            emptycells <- gridrows*gridcols-ngridplots
          }else if(rcfixed=="cols" || rcfixed=="none"){
            gridrows <- rows+1
            emptycells <- gridrows*gridcols-ngridplots
          }else{
            gridcols <- cols+1
            emptycells <- gridrows*gridcols-ngridplots
          }
        }
        if(emptycells<=gridcols){
          legend.pos <- c(c(gridcols:1)*gridrows)[emptycells:1]
        }else{
          #           v1 <- c(gridcols:1)*gridrows
          #           v2 <- 0:(floor(emptycells/gridcols)-1)
          #           legend.pos <- c(rep(v1, times=length(v2))-rep(v2, each=length(v1)))
          legend.pos <- c(gridcols:1)*gridrows
        }
        # byrow=TRUE
      }else{
        if(emptycells==0){
          if(rcfixed=="both"){
            warning(paste0("Legend does not fit to requested grid with ", rows, " rows and ", cols, " columns. A row was added."))
            gridrows <- rows+1
            emptycells <- gridrows*gridcols-ngridplots
          }else if(rcfixed=="cols" || rcfixed=="none"){
            gridrows <- rows+1
            emptycells <- gridrows*gridcols-ngridplots
          }else{
            gridcols <- cols+1
            emptycells <- gridrows*gridcols-ngridplots
          }
        }
        if(emptycells<=gridcols){
          legend.pos <- c((gridrows*gridcols-emptycells+1):(gridrows*gridcols))
        }else{
          legend.pos <- c((gridrows*gridcols-floor(emptycells/gridcols)*gridcols+1):(gridrows*gridcols))
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
            warning(paste0("Legend does not fit to requested grid with ", rows, " rows and ", cols, " columns. A column was added."))
            gridcols <- cols+1
            emptycells <- gridrows*gridcols-ngridplots
          }else if(rcfixed=="rows" || rcfixed=="none"){
            gridcols <- cols+1
            emptycells <- gridrows*gridcols-ngridplots
          }else{
            gridrows <- rows+1
            emptycells <- gridrows*gridcols-ngridplots
          }
        }
        if(emptycells<=gridrows){
          legend.pos <- c((gridrows*gridcols-emptycells+1):(gridrows*gridcols))
        }else{
          legend.pos <- c((gridrows*gridcols-floor(emptycells/gridrows)*gridrows+1):(gridrows*gridcols))
        }
        # byrow=TRUE
      }else{
        if(emptycells==0){
          if(rcfixed=="both"){
            warning(paste0("Legend does not fit to requested grid with ", rows, " rows and ", cols, " columns. A column was added."))
            gridcols <- cols+1
            emptycells <- gridrows*gridcols-ngridplots
          }else if(rcfixed=="rows" || rcfixed=="none"){
            gridcols <- cols+1
            emptycells <- gridrows*gridcols-ngridplots
          }else{
            gridrows <- rows+1
            emptycells <- gridrows*gridcols-ngridplots
          }
        }
        if(emptycells<=gridrows){
          legend.pos <- c(c(gridrows:1)*gridcols)[emptycells:1]
        }else{
          #           v1 <- c(gridrows:1)*gridcols
          #           v2 <- 0:(floor(emptycells/gridcols)-1)
          #           legend.pos <- c(rep(v1, times=length(v2))-rep(v2, each=length(v1)))
          legend.pos <- c(gridrows:1)*gridcols
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
        legend.rows <- sapply(lapply(x[[1]]$obs, "alphabet"), "length")
      }else if(length(ncol.legend)==1 && x[[1]]$nplots>1){
        legend.rows <- ceiling(sapply(lapply(x[[1]]$obs, "alphabet"), "length")/ncol.legend)
        ncol.legend <- rep(ncol.legend, x[[1]]$nplots)
      }else if(length(ncol.legend)<x[[1]]$nplots){
        warning(paste0("The length of ncol.legend does not match the number of requested plots. The last were arranged in 1 column."))
        ncol.legend <- c(ncol.legend, rep(1,(x[[1]]$nplots-length(ncol.legend))))
        legend.rows <- ceiling(sapply(lapply(x[[1]]$obs, "alphabet"), "length")/ncol.legend)
      }else if(length(ncol.legend)>x[[1]]$nplots){
        warning(paste0("The length of ncol.legend does not match the number of requested plots. Only the first ", x[[1]]$nplots, " arguments of \"ncol.legend\" were used."))
        legend.rows <- ceiling(sapply(lapply(x[[1]]$obs, "alphabet"), "length")/ncol.legends[1:x[[1]]$nplots])
      }
      
      if(!is.na(title.legend) && title.legend!=FALSE && !is.null(title.legend)){
        legend.rows <- legend.rows+1
      }
    }
  }
  
  
  
  
  # Cells' proportions
  if(!is.numeric(row.prop) && row.prop=="auto"){
    row.prop <- rep(1/gridrows, gridrows)
  }
  if(!is.numeric(col.prop) && col.prop=="auto"){
    col.prop <- rep(1/gridcols, gridcols)
  }
  if(length(row.prop)!=gridrows){
    warning("The length of the vector provided for row.prop does not match the number of rows in the plot. Argument row.prop was changed to \"auto\".")
    row.prop <- rep(1/gridrows, gridrows)
  }
  if(length(col.prop)!=gridcols){
    warning("The length of the vector provided for col.prop does not match the number of columns in the plot. Argument col.prop was changed to \"auto\".")
    col.prop <- rep(1/gridcols, gridcols)
  }
  

  if(byrow==FALSE){
    plotrows <- rep(c(1:gridrows), times=gridcols)
    plotcols <- rep(c(1:gridcols), each=gridrows)
    plotgrid <- cbind(plotrows,plotcols)
    if(!is.na(withlegend) && withlegend!=FALSE){
      lpos <- matrix(c(1:nrow(plotgrid)), byrow=FALSE, nrow=gridrows)
      legendplace <- matrix(c(lpos[,] %in% legend.pos), byrow=FALSE, nrow=gridrows)
      plotplace <- matrix(c(!(lpos[,] %in% legend.pos)), byrow=FALSE, nrow=gridrows)
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
    plotrows <- rep(c(1:gridrows), each=gridcols)
    plotcols <- rep(c(1:gridcols), times=gridrows)
    plotgrid <- cbind(plotrows,plotcols)
    if(!is.na(withlegend) && withlegend!=FALSE){
      lpos <- matrix(c(1:nrow(plotgrid)), byrow=TRUE, nrow=gridrows)
      legendplace <- matrix(c(lpos[,] %in% legend.pos), byrow=FALSE, nrow=gridrows)
      plotplace <- matrix(c(!(lpos[,] %in% legend.pos)), byrow=FALSE, nrow=gridrows)
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
                            grid.layout(gridrows,gridcols,
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
                                            heights=unit((legend.rows/sum(legend.rows)), 
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
                                            widths=unit((legend.rows/sum(legend.rows)), 
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

