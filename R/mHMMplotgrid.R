mHMMplotgrid <- function(x, which.plots = NULL, nrow=NA, ncol=NA, byrow=FALSE,
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
  on.exit(opar, add = TRUE)
  on.exit(graphics::layout(1), add = TRUE)
  
  divmodels <- separate_mhmm(x)
  
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
    which.plots <- 1:x$n_clusters
  }
  ngridplots <- length(which.plots)
  
  
  if(is.na(nrow) && is.na(ncol)){
    nrow <- ceiling(sqrt(ngridplots))
    ncol <- ceiling(ngridplots/nrow)
  }else if(is.na(nrow)){
    nrow <- ceiling(ngridplots/ncol)
  }else if(is.na(ncol)){
    ncol <- ceiling(ngridplots/nrow)
  }
  
  # Number of columns in legends
  if(!is.na(withlegend) && withlegend==TRUE){ 
    if(length(ncol.legend)==1 && ncol.legend=="auto"){
      ncol.legend <- rep(1, ngridplots)
    }else if(length(ncol.legend)==1 && x$n_clusters>1){
      ncol.legend <- rep(ncol.legend, ngridplots)
    }else if(length(ncol.legend)<ngridplots){
      warning(paste0("The length of ncol.legend does not match the number of requested plots. The last were arranged in 1 column."))
      ncol.legend <- c(ncol.legend, rep(1,(ngridplots-length(ncol.legend))))
    }else if(length(ncol.legend)>ngridplots){
      warning(paste0("The length of ncol.legend does not match the number of requested plots. Only the first ", x$n_clusters, " arguments of \"ncol.legend\" were used."))
      ncol.legend <- ncol.legend[1:ngridplots]
    }
  }
  
  
  
  
  
  # Cells' proportions
  if(!is.numeric(row.prop) && row.prop=="auto"){
    row.prop <- rep(1/nrow, nrow)
  }
  if(!is.numeric(col.prop) && col.prop=="auto"){
    col.prop <- rep(1/ncol, ncol)
  }
  if(length(row.prop)!=nrow){
    warning("The length of the vector provided for row.prop does not match the number of nrow in the plot. Argument row.prop was changed to \"auto\".")
    row.prop <- rep(1/nrow, nrow)
  }
  if(length(col.prop)!=ncol){
    warning("The length of the vector provided for col.prop does not match the number of columns in the plot. Argument col.prop was changed to \"auto\".")
    col.prop <- rep(1/ncol, ncol)
  }
  
  # Plotting order for layout
  if(!is.na(withlegend) && withlegend!=FALSE){
    if(byrow==FALSE){
      plotlayout <- matrix(c(1:ngridplots, rep(0,nrow*ncol-ngridplots)), nrow=nrow)
      legendlayout <- matrix(c((ngridplots+1):(2*ngridplots), rep(0,nrow*ncol-ngridplots)), nrow=nrow)
      if(withlegend=="right"){
        # Matrix for layout
        lmatrix <- cbind(plotlayout[,1], legendlayout[,1])
        if(ncol>1){
          for(i in 2:ncol){
            lmatrix <- cbind(lmatrix,plotlayout[,i], legendlayout[,i])
          }
        }
        cprops <- c(col.prop[1]*(1-legend.prop),col.prop[1]*legend.prop)
        if(ncol>1){
          for(i in 2:ncol){
            cprops <- c(cprops,col.prop[i]*(1-legend.prop),col.prop[i]*legend.prop)
          }
        }
        rprops <- row.prop
      }else if(withlegend=="left"){
        lmatrix <- cbind(legendlayout[,1], plotlayout[,1])
        if(ncol>1){
          for(i in 2:ncol){
            lmatrix <- cbind(lmatrix,legendlayout[,i], plotlayout[,i])
          }
        }
        cprops <- c(col.prop[1]*legend.prop,col.prop[1]*(1-legend.prop))
        if(ncol>1){
          for(i in 2:ncol){
            cprops <- c(cprops,col.prop[i]*legend.prop,col.prop[i]*(1-legend.prop))
          }
        }
        rprops <- row.prop
      }else if(withlegend=="bottom"){
        lmatrix <- rbind(plotlayout[1,], legendlayout[1,])
        if(nrow>1){
          for(i in 2:nrow){
            lmatrix <- rbind(lmatrix, plotlayout[i,], legendlayout[i,])
          }
        }
        rprops <- c(row.prop[1]*(1-legend.prop),row.prop[1]*legend.prop)
        if(nrow>1){
          for(i in 2:nrow){
            rprops <- c(rprops,row.prop[i]*(1-legend.prop),row.prop[i]*legend.prop)
          }
        }
        cprops <- col.prop
        # withlegend=="top"
      }else{
        lmatrix <- rbind(legendlayout[1,], plotlayout[1,])
        if(nrow>1){
          for(i in 2:nrow){
            lmatrix <- rbind(lmatrix,legendlayout[i,], plotlayout[i,])
          }
        }
        rprops <- c(row.prop[1]*legend.prop,row.prop[1]*(1-legend.prop))
        if(nrow>1){
          for(i in 2:nrow){
            rprops <- c(rprops,row.prop[i]*legend.prop,row.prop[i]*(1-legend.prop))
          }
        }
        cprops <- col.prop
      }
      # byrow=TRUE
    }else{
      plotlayout <- matrix(c(1:ngridplots, rep(0,nrow*ncol-ngridplots)), nrow=nrow, byrow=TRUE)
      legendlayout <- matrix(c((ngridplots+1):(2*ngridplots), rep(0,nrow*ncol-ngridplots)), nrow=nrow, byrow=TRUE)
      if(nrow*ncol>ngridplots){
        plotlayout[plotlayout>ngridplots] <- 0
        legendlayout[legendlayout>(2*ngridplots)] <- 0
      }
      if(withlegend=="right"){
        # Matrix for layout
        lmatrix <- cbind(plotlayout[,1], legendlayout[,1])
        if(ncol>1){
          for(i in 2:ncol){
            lmatrix <- cbind(lmatrix,plotlayout[,i], legendlayout[,i])
          }
        }
        cprops <- c(col.prop[1]*(1-legend.prop),col.prop[1]*legend.prop)
        if(ncol>1){
          for(i in 2:ncol){
            cprops <- c(cprops,col.prop[i]*(1-legend.prop),col.prop[i]*legend.prop)
          }
        }
        rprops <- row.prop
      }else if(withlegend=="left"){
        lmatrix <- cbind(legendlayout[,1], plotlayout[,1])
        if(ncol>1){
          for(i in 2:ncol){
            lmatrix <- cbind(lmatrix,legendlayout[,i], plotlayout[,i])
          }
        }
        cprops <- c(col.prop[1]*legend.prop,col.prop[1]*(1-legend.prop))
        if(ncol>1){
          for(i in 2:ncol){
            cprops <- c(cprops,col.prop[i]*legend.prop,col.prop[i]*(1-legend.prop))
          }
        }
        rprops <- row.prop 
      }else if(withlegend=="bottom"){
        lmatrix <- rbind(plotlayout[1,], legendlayout[1,])
        if(nrow>1){
          for(i in 2:nrow){
            lmatrix <- rbind(lmatrix, plotlayout[i,], legendlayout[i,])
          }
        }
        rprops <- c(row.prop[1]*(1-legend.prop),row.prop[1]*legend.prop)
        if(nrow>1){
          for(i in 2:nrow){
            rprops <- c(rprops,row.prop[i]*(1-legend.prop),row.prop[i]*legend.prop)
          }
        }
        cprops <- col.prop 
        # "top"
      }else{
        lmatrix <- rbind(legendlayout[1,], plotlayout[1,])
        if(nrow>1){
          for(i in 2:nrow){
            lmatrix <- rbind(lmatrix,legendlayout[i,], plotlayout[i,])
          }
        }
        rprops <- c(row.prop[1]*legend.prop,row.prop[1]*(1-legend.prop))
        if(nrow>1){
          for(i in 2:nrow){
            rprops <- c(rprops,row.prop[i]*legend.prop,row.prop[i]*(1-legend.prop))
          }
        }
        cprops <- col.prop 
      }
    }
    # No legends
  }else{
    if(byrow==FALSE){
      plotlayout <- matrix(c(1:ngridplots, rep(0,nrow*ncol-ngridplots)), nrow=nrow)
      legendlayout <- matrix(c((ngridplots+1):(2*ngridplots), rep(0,nrow*ncol-ngridplots)), nrow=nrow)
      cprops <- col.prop
      rprops <- row.prop
      # byrow=TRUE
    }else{
      plotlayout <- matrix(c(1:ngridplots, rep(0,nrow*ncol-ngridplots)), nrow=nrow, byrow=TRUE)
      legendlayout <- matrix(c((ngridplots+1):(2*ngridplots), rep(0,nrow*ncol-ngridplots)), nrow=nrow, byrow=TRUE)
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
  
  graphics::layout(1)
}
