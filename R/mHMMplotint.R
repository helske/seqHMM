mHMMplotint <- function(x, ask = FALSE, which.plots = NULL, layout="horizontal", pie=TRUE, 
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
  
  oldPar <- par(no.readonly=TRUE)
  on.exit(par(oldPar))
  on.exit(par(mfrow=c(1,1)))
  
  divmodels <- sepMixHMM(x)
  
  if (is.null(which.plots) && !ask){
    which.plots <- 1:x$numberOfModels
  }
  
  if(!is.null(which.plots)){
    if(any(!is.numeric(which.plots)) || any(!(which.plots %in% 1:x$numberOfModels))){
      stop(paste0("The which.plot argument only accepts numerical values between 1 and ", x$numberOfModels, "."))
    }
  }else if(!ask && is.null(which.plots)){
    which.plots <- 1:x$numberOfModels
  }
  
  if (ask && is.null(which.plots)) {
    tmenu <- x$modelNames
    repeat {
      pick <- menu(tmenu, title = "\nMake a model selection (or 0 to exit):\n")
      if(pick==0){
        return(invisible())
      }else{
        plot.HMModel(divmodels[[pick]], layout, pie, 
                     vertex.size, vertex.label, 
                     vertex.label.dist, vertex.label.pos,
                     vertex.label.family,
                     loops, edge.curved, edge.label, 
                     edge.width, cex.edge.width, 
                     edge.arrow.size, edge.label.family,
                     label.signif, label.scientific, label.max.length,
                     trim, 
                     combine.slices, combined.slice.color, 
                     combined.slice.label,
                     withlegend, ltext, legend.prop, 
                     cex.legend, ncol.legend, cpal, ...)
      }
    }
  }else if (ask && !is.null(which.plots)) {
    tmenu <- which.plots
    repeat {
      pick <- menu(tmenu, title = "\nMake a model selection (or 0 to exit):\n")
      if(pick==0){
        return(invisible())
      }else{
        plot.HMModel(divmodels[[pick]], layout, pie, 
                     vertex.size, vertex.label, 
                     vertex.label.dist, vertex.label.pos,
                     vertex.label.family,
                     loops, edge.curved, edge.label, 
                     edge.width, cex.edge.width, 
                     edge.arrow.size, edge.label.family,
                     label.signif, label.scientific, label.max.length,
                     trim, 
                     combine.slices, combined.slice.color, 
                     combined.slice.label,
                     withlegend, ltext, legend.prop, 
                     cex.legend, ncol.legend, cpal, ...)
      }
    }
  }else {
    ask <- length(which.plots) > 1
    for (i in which.plots) {
      plot.HMModel(divmodels[[i]], layout, pie, 
                   vertex.size, vertex.label, 
                   vertex.label.dist, vertex.label.pos,
                   vertex.label.family,
                   loops, edge.curved, edge.label, 
                   edge.width, cex.edge.width, 
                   edge.arrow.size, edge.label.family,
                   label.signif, label.scientific, label.max.length,
                   trim, 
                   combine.slices, combined.slice.color, 
                   combined.slice.label,
                   withlegend, ltext, legend.prop, 
                   cex.legend, ncol.legend, cpal, ...)
      if (ask) {
        op <- par(ask = TRUE)
      }
    }
  }
  invisible()
  par(oldPar)
}