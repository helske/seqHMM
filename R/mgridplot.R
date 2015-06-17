#' Plot Mixture Hidden Markov Models in a Grid
#' 
#' Function \code{mgridplot} plots multiple \code{mixHMModel} graphs to a
#' grid.
#' 
#' 
#' 
#' @export
#'   
#' @param x A \code{mixHMModel} object.
#'   
#' @param rows,cols Optional arguments to arrange cluster plots.
#'   
#' @param byrow Controls the order of plotting. Defaults to \code{FALSE}, i.e. plots
#'   are arranged columnwise.
#'   
#' @param withlegend defines if and where the legend of the state colors is 
#'   plotted. Possible values include \code{"bottom"} (the default), 
#'   \code{"top"}, \code{"left"}, and \code{"right"}. 
#'   \code{FALSE} omits the legend.
#'   
#' @param legend.pos Defines the positions of the legend boxes relative to the
#'   cluster graphs. One of \code{"bottomright"}, \code{"bottom"}, \code{"bottomleft"}, 
#'   \code{"left"}, \code{"topleft"}, \code{"top"}, \code{"topright"}, 
#'   \code{"right"} and \code{"center"} (the default).
#'   
#' @param legend.prop The proportion of legends
#'   
#' @param title.legend The titles for the legend boxes. The default \code{"auto"} takes
#'   the titles from the cluster names in x. \code{NA} prints no title.
#'   
#' @param ncol.legend (A vector of) the number of columns for the legend(s). The
#'   default \code{"auto"} creates one column for each legend.
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
#' @param ... Other parameters passed on to \code{\link{plot.HMModel}}.
#'
#' @seealso \code{\link{ssp}} for defining the plot before using
#'   \code{gridplot}, and \code{\link{plot.ssp}} for plotting only one ssp object.
#'   
#' @examples 
#' require(TraMineR)
#' 
#' data(biofam)
#' biofam <- biofam[complete.cases(biofam[c(2:4)]),]
#' biofam <- biofam[1:500,]
#' 
#' ## Building one channel per type of event left, children or married
#' bf <- as.matrix(biofam[, 10:25])
#' children <-  bf==4 | bf==5 | bf==6
#' married <- bf == 2 | bf== 3 | bf==6
#' left <- bf==1 | bf==3 | bf==5 | bf==6 | bf==7
#' 
#' children[children==TRUE] <- "Children"
#' children[children==FALSE] <- "Childless"
#' # Divorced parents
#' div <- bf[(rowSums(bf==7)>0 & rowSums(bf==5)>0) | 
#'             (rowSums(bf==7)>0 & rowSums(bf==6)>0),]
#' children[rownames(bf) %in% rownames(div) & bf==7] <- "Children"
#' 
#' married[married==TRUE] <- "Married"
#' married[married==FALSE] <- "Single"
#' married[bf==7] <- "Divorced"
#' 
#' left[left==TRUE] <- "Left home"
#' left[left==FALSE] <- "With parents"
#' # Divorced living with parents (before divorce)
#' wp <- bf[(rowSums(bf==7)>0 & rowSums(bf==2)>0 & rowSums(bf==3)==0 &  rowSums(bf==5)==0 &  rowSums(bf==6)==0) | 
#'            (rowSums(bf==7)>0 & rowSums(bf==4)>0 & rowSums(bf==3)==0 &  rowSums(bf==5)==0 &  rowSums(bf==6)==0),]
#' left[rownames(bf) %in% rownames(wp) & bf==7] <- "With parents"
#' 
#' ## Building sequence objects
#' child.seq <- seqdef(children, start=15)
#' marr.seq <- seqdef(married, start=15)
#' left.seq <- seqdef(left, start=15)
#' 
#' ## Starting values for emission probabilities
#' 
#' # Cluster 1
#' alphabet(child.seq) # Checking for the order of observed states
#' B1_child <- matrix(c(0.99, 0.01, # High probability for childless
#'                      0.99, 0.01,
#'                      0.99, 0.01,
#'                      0.99, 0.01), nrow=4, ncol=2, byrow=TRUE)
#' 
#' alphabet(marr.seq)                      
#' B1_marr <- matrix(c(0.01, 0.01, 0.98, # High probability for single
#'                     0.01, 0.01, 0.98,
#'                     0.01, 0.98, 0.01, # High probability for married
#'                     0.98, 0.01, 0.01), # High probability for divorced
#'                     nrow=4, ncol=3, byrow=TRUE)                   
#' 
#' alphabet(left.seq)
#' B1_left <- matrix(c(0.01, 0.99, # High probability for living with parents
#'                     0.99, 0.01, # High probability for having left home
#'                     0.99, 0.01
#'                     0.99, 0.01), nrow=4, ncol=2, byrow=TRUE)
#' 
#' B2_marr <- matrix(c(0.01, 0.01, 0.98, # High probability for single
#'                      0.01, 0.01, 0.98,
#'                      0.01, 0.98, 0.01, # High probability for married
#'                      0.29, 0.7, 0.01),
#'                    nrow=4, ncol=3, byrow=TRUE)                   
#' 
#' B2_left <- matrix(c(0.01, 0.99, # High probability for living with parents
#'                      0.99, 0.01,
#'                      0.99, 0.01,
#'                      0.99, 0.01), nrow=4, ncol=2, byrow=TRUE) 
#' 
#' # Sinkkuvanhemmat ja kotona asuvat yhdessä
#' B3_child <- matrix(c(0.99, 0.01, # High probability for childless
#'                       0.99, 0.01,
#'                       0.01, 0.99,
#'                       0.99, 0.01,
#'                       0.01, 0.99,
#'                       0.01, 0.99), nrow=6, ncol=2, byrow=TRUE)
#' 
#' B3_marr <- matrix(c(0.01, 0.01, 0.98, # High probability for single
#'                      0.01, 0.01, 0.98,
#'                      0.01, 0.01, 0.98,
#'                      0.01, 0.98, 0.01,
#'                      0.01, 0.98, 0.01, # High probability for married
#'                      0.98, 0.01, 0.01), # High probability for divorced
#'                    nrow=6, ncol=3, byrow=TRUE)                   
#' 
#' B3_left <- matrix(c(0.01, 0.99, # High probability for living with parents
#'                      0.99, 0.01,
#'                      0.50, 0.50,
#'                      0.01, 0.99,
#'                      0.99, 0.01,
#'                      0.99, 0.01), nrow=6, ncol=2, byrow=TRUE) 
#' 
#' # Initial values for transition matrices
#' A1 <- matrix(c(0.8,   0.16, 0.03, 0.01,
#'                0,    0.9, 0.07, 0.03, 
#'                0,      0,  0.9,  0.1, 
#'                0,      0,    0,    1), 
#'              nrow=4, ncol=4, byrow=TRUE)
#' 
#' A2 <- matrix(c(0.8, 0.10, 0.05,  0.03, 0.01, 0.01,
#'                0,    0.7,  0.1,   0.1, 0.05, 0.05,
#'                0,      0,  0.85, 0.01,  0.1, 0.04,
#'                0,      0,    0,   0.9, 0.05, 0.05,
#'                0,      0,    0,     0,  0.9,  0.1,
#'                0,      0,    0,     0,    0,    1), 
#'              nrow=6, ncol=6, byrow=TRUE)
#' 
#' # Initial values for initial state probabilities 
#' initialProbs1 <- c(0.9, 0.07, 0.02, 0.01)
#' initialProbs2 <- c(0.9, 0.04, 0.03, 0.01, 0.01, 0.01)
#' 
#' # Creating covariate swiss
#' bio$swiss <- bio$nat_1_02=="Switzerland"
#' bio$swiss[bio$swiss==TRUE] <- "Swiss"
#' bio$swiss[bio$swiss==FALSE] <- "Other"
#' 
#' # Build mixture HMM
#' bmHMM <- buildMixHMM(observations=list(child.seq, marr.seq, left.seq), 
#'                        transitionMatrix=list(A1,A2,A1), 
#'                        emissionMatrix=list(list(B1_child, B1_marr, B1_left),
#'                                            list(B2_child, B2_marr, B2_left),
#'                                            list(B3_child, B3_marr, B3_left)),
#'                        initialProbs=list(initialProbs1, initialProbs2,
#'                                          initialProbs1), 
#'                        formula=~sex*birthyr+sex*swiss, data=bio,
#'                        clusterNames=c("Cluster 1", "Cluster 2", "Cluster 3"),
#'                        channelNames=c("Parenthood", "Marriage", "Left home"))
#' 
#' mHMM <- fitMixHMM(bmHMM)
#' 
#' mgridplot(mHMM$model)
#' 
#' mgridplot(mHMM$model, cols=2, withlegend="right", xlim=c(0,4), legend.prop=0.3)
#'                   

mgridplot <- function(x, rows=NA, cols=NA, byrow=FALSE,
                      withlegend="bottom", legend.pos="center", 
                      legend.prop=0.5, title.legend="auto",
                      ncol.legend="auto", 
                      with.missing.legend="auto",                       
                      row.prop="auto", col.prop="auto", cex.legend=1, ...){
  
  # Checking for the class of x
  if(!inherits(x, "mixHMModel")){
    stop("Your object provided for the x argument is not a mixHMModel object. Use buildMixHMM and fitMixHMM to create one.")
  }
  divmodels <- sepMixHMM(x)
  
  plot.new()
  opar <- par(no.readonly=TRUE)
  on.exit(opar)
  
  
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
  ngridplots <- x$numberOfClusters
  
  
  if(is.na(rows) && is.na(cols)){
    rows <- ceiling(sqrt(ngridplots))
    cols <- ceiling(ngridplots/rows)
  }else if(is.na(rows)){
    rows <- ceiling(ngridplots/cols)
  }else if(is.na(cols)){
    cols <- ceiling(ngridplots/rows)
  }
  
  # Legend titles
  if(!is.na(withlegend) && withlegend!=FALSE){
    if(!is.na(title.legend) && title.legend!=FALSE && !is.null(title.legend)){
      if(length(title.legend)>1 || (length(title.legend)==1 && title.legend!="auto")){
        if(length(title.legend)!=x$numberOfClusters){
          warning("The length of the vector provided for title.legend does not match the number of clusters (graphs). Argument title.legend was set to \"auto\".")
          title.legend=="auto"
        }
      }
      if(length(title.legend==1 && title.legend=="auto")){
        title.legend <- x$clusterNames
      }
    }
  }
  
  # Number of columns in legends
  if(!is.na(withlegend) && withlegend==TRUE){ 
    if(length(ncol.legend)==1 && ncol.legend=="auto"){
      ncol.legend <- rep(1, x$numberOfClusters)
    }else if(length(ncol.legend)==1 && x$numberOfClusters>1){
      ncol.legend <- rep(ncol.legend, x$numberOfClusters)
    }else if(length(ncol.legend)<x$numberOfClusters){
      warning(paste0("The length of ncol.legend does not match the number of clusters The last were arranged in 1 column."))
      ncol.legend <- c(ncol.legend, rep(1,(x$numberOfClusters-length(ncol.legend))))
    }else if(length(ncol.legend)>x$numberOfClusters){
      warning(paste0("The length of ncol.legend does not match the number of clusters. Only the first ", x$numberOfClusters, " arguments of \"ncol.legend\" were used."))
      ncol.legend <- ncol.legend[1:x$numberOfClusters]
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
  for(p in 1:ngridplots){
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
  for(p in 1:ngridplots){
    eval(HMMcalls[[p]]$plotcall)
  }
  
  # Plotting legends
  if(withlegend!=FALSE){
    for(p in 1:ngridplots){
      eval(HMMcalls[[p]]$legendcall)
    }
  }
  
  
  # par(opar)
  graphics::layout(1)
}
