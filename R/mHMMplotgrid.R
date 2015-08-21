mHMMplotgrid <- function(x, which.plots = NULL, rows=NA, cols=NA, byrow=FALSE,
                         row.prop="auto", col.prop="auto", layout="horizontal", pie=TRUE, 
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
                         withlegend="bottom", legend.pos="center", ltext=NULL, legend.prop=0.5, 
                         cex.legend=1, ncol.legend="auto", cpal="auto", ...){
  
  plot.new()
  opar <- par(no.readonly=TRUE)
  on.exit(opar)
  
  divmodels <- sepMixHMM(x)
  
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
  if(is.null(which.plots)){
    which.plots <- 1:x$numberOfClusters
  }
  ngridplots <- length(which.plots)
  
  
  if(is.na(rows) && is.na(cols)){
    rows <- ceiling(sqrt(ngridplots))
    cols <- ceiling(ngridplots/rows)
  }else if(is.na(rows)){
    rows <- ceiling(ngridplots/cols)
  }else if(is.na(cols)){
    cols <- ceiling(ngridplots/rows)
  }
  
  # Number of columns in legends
  if(!is.na(withlegend) && withlegend==TRUE){ 
    if(length(ncol.legend)==1 && ncol.legend=="auto"){
      ncol.legend <- rep(1, ngridplots)
    }else if(length(ncol.legend)==1 && x$numberOfClusters>1){
      ncol.legend <- rep(ncol.legend, ngridplots)
    }else if(length(ncol.legend)<ngridplots){
      warning(paste0("The length of ncol.legend does not match the number of requested plots. The last were arranged in 1 column."))
      ncol.legend <- c(ncol.legend, rep(1,(ngridplots-length(ncol.legend))))
    }else if(length(ncol.legend)>ngridplots){
      warning(paste0("The length of ncol.legend does not match the number of requested plots. Only the first ", x$numberOfClusters, " arguments of \"ncol.legend\" were used."))
      ncol.legend <- ncol.legend[1:ngridplots]
    }
  }
  
  
  
  
  
  # Cells' proportions
  if(!is.numeric(row.prop) && row.prop=="auto"){
    row.prop <- rep(1/rows, rows)
  }
  if(!is.numeric(col.prop) && col.prop=="auto"){
    col.prop <- rep(1/cols, cols)
  }
  if(length(row.prop)!=rows){
    warning("The length of the vector provided for row.prop does not match the number of rows in the plot. Argument row.prop was changed to \"auto\".")
    row.prop <- rep(1/rows, rows)
  }
  if(length(col.prop)!=cols){
    warning("The length of the vector provided for col.prop does not match the number of columns in the plot. Argument col.prop was changed to \"auto\".")
    col.prop <- rep(1/cols, cols)
  }
  
  # Plotting order for layout
  if(!is.na(withlegend) && withlegend!=FALSE){
    if(byrow==FALSE){
      plotlayout <- matrix(c(1:ngridplots, rep(0,rows*cols-ngridplots)), nrow=rows)
      legendlayout <- matrix(c((ngridplots+1):(2*ngridplots), rep(0,rows*cols-ngridplots)), nrow=rows)
      if(withlegend=="right"){
        # Matrix for layout
        lmatrix <- cbind(plotlayout[,1], legendlayout[,1])
        if(cols>1){
          for(i in 2:cols){
            lmatrix <- cbind(lmatrix,plotlayout[,i], legendlayout[,i])
          }
        }
        cprops <- c(col.prop[1]*(1-legend.prop),col.prop[1]*legend.prop)
        if(cols>1){
          for(i in 2:cols){
            cprops <- c(cprops,col.prop[i]*(1-legend.prop),col.prop[i]*legend.prop)
          }
        }
        rprops <- row.prop
      }else if(withlegend=="left"){
        lmatrix <- cbind(legendlayout[,1], plotlayout[,1])
        if(cols>1){
          for(i in 2:cols){
            lmatrix <- cbind(lmatrix,legendlayout[,i], plotlayout[,i])
          }
        }
        cprops <- c(col.prop[1]*legend.prop,col.prop[1]*(1-legend.prop))
        if(cols>1){
          for(i in 2:cols){
            cprops <- c(cprops,col.prop[i]*legend.prop,col.prop[i]*(1-legend.prop))
          }
        }
        rprops <- row.prop
      }else if(withlegend=="bottom"){
        lmatrix <- rbind(plotlayout[1,], legendlayout[1,])
        if(rows>1){
          for(i in 2:rows){
            lmatrix <- rbind(lmatrix, plotlayout[i,], legendlayout[i,])
          }
        }
        rprops <- c(row.prop[1]*(1-legend.prop),row.prop[1]*legend.prop)
        if(rows>1){
          for(i in 2:rows){
            rprops <- c(rprops,row.prop[i]*(1-legend.prop),row.prop[i]*legend.prop)
          }
        }
        cprops <- col.prop
        # withlegend=="top"
      }else{
        lmatrix <- rbind(legendlayout[1,], plotlayout[1,])
        if(rows>1){
          for(i in 2:rows){
            lmatrix <- rbind(lmatrix,legendlayout[i,], plotlayout[i,])
          }
        }
        rprops <- c(row.prop[1]*legend.prop,row.prop[1]*(1-legend.prop))
        if(rows>1){
          for(i in 2:rows){
            rprops <- c(rprops,row.prop[i]*legend.prop,row.prop[i]*(1-legend.prop))
          }
        }
        cprops <- col.prop
      }
      # byrow=TRUE
    }else{
      plotlayout <- matrix(c(1:ngridplots, rep(0,rows*cols-ngridplots)), nrow=rows, byrow=TRUE)
      legendlayout <- matrix(c((ngridplots+1):(2*ngridplots), rep(0,rows*cols-ngridplots)), nrow=rows, byrow=TRUE)
      if(rows*cols>ngridplots){
        plotlayout[plotlayout>ngridplots] <- 0
        legendlayout[legendlayout>(2*ngridplots)] <- 0
      }
      if(withlegend=="right"){
        # Matrix for layout
        lmatrix <- cbind(plotlayout[,1], legendlayout[,1])
        if(cols>1){
          for(i in 2:cols){
            lmatrix <- cbind(lmatrix,plotlayout[,i], legendlayout[,i])
          }
        }
        cprops <- c(col.prop[1]*(1-legend.prop),col.prop[1]*legend.prop)
        if(cols>1){
          for(i in 2:cols){
            cprops <- c(cprops,col.prop[i]*(1-legend.prop),col.prop[i]*legend.prop)
          }
        }
        rprops <- row.prop
      }else if(withlegend=="left"){
        lmatrix <- cbind(legendlayout[,1], plotlayout[,1])
        if(cols>1){
          for(i in 2:cols){
            lmatrix <- cbind(lmatrix,legendlayout[,i], plotlayout[,i])
          }
        }
        cprops <- c(col.prop[1]*legend.prop,col.prop[1]*(1-legend.prop))
        if(cols>1){
          for(i in 2:cols){
            cprops <- c(cprops,col.prop[i]*legend.prop,col.prop[i]*(1-legend.prop))
          }
        }
        rprops <- row.prop 
      }else if(withlegend=="bottom"){
        lmatrix <- rbind(plotlayout[1,], legendlayout[1,])
        if(rows>1){
          for(i in 2:rows){
            lmatrix <- rbind(lmatrix, plotlayout[i,], legendlayout[i,])
          }
        }
        rprops <- c(row.prop[1]*(1-legend.prop),row.prop[1]*legend.prop)
        if(rows>1){
          for(i in 2:rows){
            rprops <- c(rprops,row.prop[i]*(1-legend.prop),row.prop[i]*legend.prop)
          }
        }
        cprops <- col.prop 
        # "top"
      }else{
        lmatrix <- rbind(legendlayout[1,], plotlayout[1,])
        if(rows>1){
          for(i in 2:rows){
            lmatrix <- rbind(lmatrix,legendlayout[i,], plotlayout[i,])
          }
        }
        rprops <- c(row.prop[1]*legend.prop,row.prop[1]*(1-legend.prop))
        if(rows>1){
          for(i in 2:rows){
            rprops <- c(rprops,row.prop[i]*legend.prop,row.prop[i]*(1-legend.prop))
          }
        }
        cprops <- col.prop 
      }
    }
    # No legends
  }else{
    if(byrow==FALSE){
      plotlayout <- matrix(c(1:ngridplots, rep(0,rows*cols-ngridplots)), nrow=rows)
      legendlayout <- matrix(c((ngridplots+1):(2*ngridplots), rep(0,rows*cols-ngridplots)), nrow=rows)
      cprops <- col.prop
      rprops <- row.prop
      # byrow=TRUE
    }else{
      plotlayout <- matrix(c(1:ngridplots, rep(0,rows*cols-ngridplots)), nrow=rows, byrow=TRUE)
      legendlayout <- matrix(c((ngridplots+1):(2*ngridplots), rep(0,rows*cols-ngridplots)), nrow=rows, byrow=TRUE)
      cprops <- col.prop
      rprops <- row.prop
    }
  }
  
  graphics::layout(lmatrix, widths=cprops, heights=rprops)
  
  
  # Plotting arguments for graphs and legends
  HMMcalls <- list()
  length(HMMcalls) <- ngridplots
  for(p in which.plots){
    if(length(ncol.legend)>1){
      ncolleg <- ncol.legend[p]
    }else{
      ncolleg <- ncol.legend
    }
    HMMcalls[[p]] <- do.call(HMMplot,args=list(divmodels[[p]], ncol.legend=ncolleg, 
                                               legend.pos=legend.pos, 
                                               cex.legend=cex.legend, 
                                               withlegend=withlegend, ...))
  }
  
  # Plotting graphs
  for(p in which.plots){
    eval(HMMcalls[[p]]$plotcall)
  }
  
  # Plotting legends
  if(withlegend!=FALSE){
    for(p in which.plots){
      eval(HMMcalls[[p]]$legendcall)
    }
  }
  
  
  # par(opar)
  graphics::layout(1)
}
