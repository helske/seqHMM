SSPlotter <- function(obs, nchannels, channel_names, nplots, 
                       legend.c.prop, legend.r.prop,
                       ylab.space, xaxis.space, xt.space,
                       mpp.seq=NULL, orderv=NULL,
                       plotxaxis,
                       mpp=NULL, plots="both", type="I", 
                       sortv=NULL, sort.channel=1, 
                       with.missing=FALSE,
                       title=NA, title.n=TRUE, cex.title=1, title.pos=1,
                       withlegend="auto", ncol.legend="auto", 
                       with.missing.legend="auto",                         
                       legend.prop=0.3, cex.legend=1,
                       mpp.color="auto", mpp.labels="auto",
                       xaxis=TRUE, xlab=NA, xtlab=NULL, xlab.pos=1,
                       ylab="auto", hidden.states.title="Hidden states", 
                       ylab.pos="auto", 
                       cex.lab=1, cex.axis=1, new=FALSE, call = match.call(), ...
){
  
  # Grid for plotting regions
  top.vp <- viewport(
    layout=grid.layout(nrow=4, ncol=3,
                       widths=unit(c(ylab.space, 1, legend.c.prop), 
                                   c("lines", "null", "npc")),
                       heights=unit(c((title.pos*cex.title+1.5), 1, 
                                      (xlab.pos*cex.lab+xaxis.space+xt.space+0.5), 
                                      legend.r.prop), 
                                    c("lines", "null", "lines", "npc"))))
  
  vptitle <- viewport(layout.pos.row=1, layout.pos.col=2, name="vptitle")
  vpylab <- viewport(layout.pos.row=2, layout.pos.col=1, name="vpylab")
  vpplot <- viewport(layout.pos.row=2, layout.pos.col=2, name="vpplot")
  if(withlegend==TRUE || withlegend=="auto" || withlegend=="right"
     || withlegend=="right.many"){
    vplegend <- viewport(layout.pos.row=2, layout.pos.col=3, name="vplegend")
  }else if(withlegend=="bottom" || withlegend=="bottom.many"){
    vplegend <- viewport(layout.pos.row=4, layout.pos.col=2, name="vplegend")
  }
  
  vpxaxis <- viewport(layout.pos.row=3, layout.pos.col=2, name="vpxaxis")
  
  if(withlegend==FALSE){
    splot <- vpTree(top.vp, vpList(vptitle, vpylab, vpplot, vpxaxis))
  }else{
    splot <- vpTree(top.vp, vpList(vptitle, vpylab, vpplot, vplegend, vpxaxis))
  }
  
  pushViewport(splot)
  
  # Grid for plots
  upViewport()
  downViewport("vpplot")
  pushViewport(viewport(layout=grid.layout(nrow=nplots, ncol=1), 
                        width=unit(0.95, "npc")))
  
  # Plotting observations
  if(plots=="both" || plots=="obs"){ 
    if(type=="I"){ 
      if(is.null(sortv)){
        if(nchannels>1){
          for(i in 1:(nchannels-1)){
            pushViewport(viewport(layout.pos.col=1, layout.pos.row=i))
            par(plt=gridPLT(), new=TRUE)
            seqplot(obs[[i]], type="I", withlegend=FALSE,
                    use.layout=FALSE, yaxis=FALSE, axes=FALSE, ylab=NA, 
                    xtlab=xtlab, ...)
            popViewport()
          }
          pushViewport(viewport(layout.pos.col=1, layout.pos.row=nchannels))
          par(plt=gridPLT(), new=TRUE)
          seqplot(obs[[nchannels]], type="I", withlegend=FALSE,
                  use.layout=FALSE, yaxis=FALSE, xaxis=plotxaxis, ylab=NA, 
                  xtlab=xtlab, cex.plot=cex.axis, ...)
          popViewport()
        }else{
          pushViewport(viewport(layout.pos.col=1, layout.pos.row=nchannels))
          par(plt=gridPLT(), new=TRUE)
          seqplot(obs, type="I", withlegend=FALSE,
                  use.layout=FALSE, yaxis=FALSE, xaxis=plotxaxis, ylab=NA, 
                  xtlab=xtlab, cex.plot=cex.axis, ...)
          popViewport()
        }
        
      }else if(length(sortv)==1 && (sortv=="from.start" || sortv=="from.end")){
        if(nchannels>1){
          for(i in 1:(nchannels-1)){
            pushViewport(viewport(layout.pos.col=1, layout.pos.row=i))
            par(plt=gridPLT(), new=TRUE)
            seqplot(obs[[i]][orderv,], type="I", withlegend=FALSE,
                    use.layout=FALSE, yaxis=FALSE, axes=FALSE, ylab=NA, 
                    xtlab=xtlab, ...)
            popViewport()
          }
          pushViewport(viewport(layout.pos.col=1, layout.pos.row=nchannels))
          par(plt=gridPLT(), new=TRUE)
          seqplot(obs[[nchannels]][orderv,], type="I", withlegend=FALSE,
                  use.layout=FALSE, yaxis=FALSE, xaxis=plotxaxis, ylab=NA, 
                  xtlab=xtlab, cex.plot=cex.axis, ...)
          popViewport()
        }else{
          pushViewport(viewport(layout.pos.col=1, layout.pos.row=nchannels))
          par(plt=gridPLT(), new=TRUE)
          seqplot(obs[orderv,], type="I", withlegend=FALSE,
                  use.layout=FALSE, yaxis=FALSE, xaxis=plotxaxis, ylab=NA, 
                  xtlab=xtlab, cex.plot=cex.axis, ...)
          popViewport()
        }

        
      }else if(length(sortv)>1){
        if(nchannels>1){
          for(i in 1:(nchannels-1)){
            pushViewport(viewport(layout.pos.col=1, layout.pos.row=i))
            par(plt=gridPLT(), new=TRUE)
            seqplot(obs[[i]], type="I", sortv=sortv, withlegend=FALSE,
                    use.layout=FALSE, yaxis=FALSE, axes=FALSE, ylab=NA, 
                    xtlab=xtlab, ...)
            popViewport()
          }
          pushViewport(viewport(layout.pos.col=1, layout.pos.row=nchannels))
          par(plt=gridPLT(), new=TRUE)
          seqplot(obs[[nchannels]], type="I", sortv=sortv, withlegend=FALSE,
                  use.layout=FALSE, yaxis=FALSE, xaxis=plotxaxis, ylab=NA, 
                  xtlab=xtlab, cex.plot=cex.axis, ...)
          popViewport()
        }else{
          pushViewport(viewport(layout.pos.col=1, layout.pos.row=nchannels))
          par(plt=gridPLT(), new=TRUE)
          seqplot(obs, type="I", sortv=sortv, withlegend=FALSE,
                  use.layout=FALSE, yaxis=FALSE, xaxis=plotxaxis, ylab=NA, 
                  xtlab=xtlab, cex.plot=cex.axis, ...)
          popViewport()
        }
        
      }
      if(plots=="obs"){
        # Close viewport "vpplot"
        popViewport(2)
      }
    }else{  # if type=="d"
      if(nchannels>1){
        for(i in 1:(nchannels-1)){
          pushViewport(viewport(layout.pos.col=1, layout.pos.row=i))
          par(plt=gridPLT(), new=TRUE)
          seqplot(obs[[i]], type="d", withlegend=FALSE,
                  use.layout=FALSE, yaxis=FALSE, axes=FALSE, ylab=NA, 
                  with.missing=with.missing, xtlab=xtlab, ...)
          popViewport()
        }
        pushViewport(viewport(layout.pos.col=1, layout.pos.row=nchannels))
        par(plt=gridPLT(), new=TRUE)
        seqplot(obs[[nchannels]], type="d", withlegend=FALSE, xaxis=plotxaxis,
                use.layout=FALSE, yaxis=FALSE, xaxis=plotxaxis, ylab=NA, 
                with.missing=with.missing, xtlab=xtlab, cex.plot=cex.axis, ...)
        popViewport()
      }else{
        pushViewport(viewport(layout.pos.col=1, layout.pos.row=nchannels))
        par(plt=gridPLT(), new=TRUE)
        seqplot(obs, type="d", withlegend=FALSE, xaxis=plotxaxis,
                use.layout=FALSE, yaxis=FALSE, xaxis=plotxaxis, ylab=NA, 
                with.missing=with.missing, xtlab=xtlab, cex.plot=cex.axis, ...)
        popViewport()
      }
      
      if(plots=="obs"){ 
        # Close viewport "vpplot"
        popViewport(2)
      }
    }
  }
  
  # Plotting the most probable paths
  if(plots=="both" || plots=="mpp"){  
    if(type=="I"){    
      if(is.null(sortv)){
        pushViewport(viewport(layout.pos.col=1, layout.pos.row=nplots))
        par(plt=gridPLT(), new=TRUE)
        seqplot(mpp.seq, type=type, sortv=sortv, withlegend=FALSE,
                use.layout=FALSE, yaxis=FALSE, axes=xaxis, ylab=NA,
                xtlab=xtlab, cex.plot=cex.axis, ...)
        popViewport()
      }else if(length(sortv)==1 && (sortv=="from.start" || sortv=="from.end")){
        pushViewport(viewport(layout.pos.col=1, layout.pos.row=nplots))
        par(plt=gridPLT(), new=TRUE)
        seqplot(mpp.seq[orderv,], type=type, withlegend=FALSE,
                use.layout=FALSE, yaxis=FALSE, axes=xaxis, ylab=NA,
                xtlab=xtlab, cex.plot=cex.axis, ...)
        popViewport()
        #       }else if(length(sortv)==1 && sortv=="mds.mpp"){
        #         pushViewport(viewport(layout.pos.col=1, layout.pos.row=nplots))
        #         par(plt=gridPLT(), new=TRUE)
        #         seqplot(mpp.seq, type=type, sortv=mds.mppscore, withlegend=FALSE,
        #                 use.layout=FALSE, yaxis=FALSE, axes=xaxis, ylab=NA,
        #                 xtlab=xtlab, cex.plot=cex.axis, ...)
        #         popViewport()
      }else if(length(sortv)>1){
        pushViewport(viewport(layout.pos.col=1, layout.pos.row=nplots))
        par(plt=gridPLT(), new=TRUE)
        seqplot(mpp.seq, type=type, sortv=sortv, withlegend=FALSE,
                use.layout=FALSE, yaxis=FALSE, axes=xaxis, ylab=NA,
                xtlab=xtlab, cex.plot=cex.axis, ...)
        popViewport()
      }   
    }else{    
      pushViewport(viewport(layout.pos.col=1, layout.pos.row=nplots))
      par(plt=gridPLT(), new=TRUE)
      seqplot(mpp.seq, type=type, withlegend=FALSE,
              use.layout=FALSE, yaxis=FALSE, axes=xaxis, ylab=NA,
              xtlab=xtlab, cex.plot=cex.axis, ...)
      popViewport()
    }
    # Close viewport "vpplot"
    popViewport(2)
  }
  
  # Plotting x label
  if(!is.na(xlab)){
    downViewport("vpxaxis")
    grid.text(xlab, y = unit(1, "lines"), 
              gp = gpar(cex = cex.lab))
    popViewport()
  }
  
  # Plotting y labels (channels and/or hidden states)
  if(length(ylab)>1 || (length(ylab)==1 && !is.na(ylab) && ylab!=FALSE)){
    downViewport("vpylab")
    pushViewport(viewport(layout=grid.layout(nrow=nplots, ncol=1)))
    if(plots=="both" || plots=="obs"){
      for(i in 1:nchannels){
        pushViewport(viewport(layout.pos.row=i, 
                              layout.pos.col=1))
        grid.text(ylab[i], x = unit(ylab.space/cex.lab-
                                      ylab.pos[i]+1, "lines"), 
                  gp = gpar(cex = cex.lab), rot=90, vjust=0.5, hjust=0.5)
        popViewport()
      }
    }
    if(plots=="both" || plots=="mpp"){
      pushViewport(viewport(layout.pos.row=nplots, 
                            layout.pos.col=1))
      grid.text(hidden.states.title, x = unit(ylab.space/cex.lab-
                                               ylab.pos[nplots]+1, 
                                             "lines"), 
                gp = gpar(cex = cex.lab), rot=90, vjust=0.5, hjust=0.5)
      popViewport()
    }
    # Close viewport "vpylab"
    popViewport(2)
  }  
  
  # Legends
  if(withlegend==TRUE || withlegend=="auto" || withlegend=="right.many" 
     || withlegend=="bottom.many"){
    # Grid for legends
    if(withlegend==TRUE || withlegend=="auto" || withlegend=="right.many"){    
      upViewport()
      downViewport("vplegend")
      pushViewport(viewport(layout=grid.layout(nrow=nplots, ncol=1)))
      lposrow <- 1:nplots
      lposcol <- rep(1, nplots)
    }else{
      upViewport()
      downViewport("vplegend")
      pushViewport(viewport(layout=grid.layout(nrow=1, ncol=nplots)))
      lposrow <- rep(1, nplots)
      lposcol <- 1:nplots
    }
    
    # Legends for observations
    if(plots=="both" || plots=="obs"){
      for(i in 1:nchannels){
        pushViewport(viewport(layout.pos.row=lposrow[i], 
                              layout.pos.col=lposcol[i]))
        pushViewport(viewport(width=unit(0.9, "npc")))
        par(plt=gridPLT(), new=TRUE)
        if(nchannels>1){
          if(withlegend=="bottom.many"){
            seqlegend(obs[[i]], fontsize=cex.legend, position="top", 
                      ncol=ncol.legend[i], with.missing=with.missing.legend)
          }else{
            seqlegend(obs[[i]], fontsize=cex.legend, position="left", 
                      ncol=ncol.legend[i], with.missing=with.missing.legend)
          }
        }else{
          if(withlegend=="bottom.many"){
            seqlegend(obs, fontsize=cex.legend, position="top", 
                      ncol=ncol.legend[i], with.missing=with.missing.legend)
          }else{
            seqlegend(obs, fontsize=cex.legend, position="left", 
                      ncol=ncol.legend[i], with.missing=with.missing.legend)
          }
        }
        popViewport(2)
      }
    }
    # Legends for most probable paths
    if(plots=="both" || plots=="mpp"){
      pushViewport(viewport(layout.pos.row=lposrow[nplots], 
                            layout.pos.col=lposcol[nplots]))
      pushViewport(viewport(width=unit(0.9, "npc")))
      par(plt=gridPLT(), new=TRUE)
      if(withlegend=="bottom.many"){
        seqlegend(mpp.seq, fontsize=cex.legend, position="top", 
                  ncol=ncol.legend[length(ncol.legend)], with.missing=with.missing.legend)
      }else{
        seqlegend(mpp.seq, fontsize=cex.legend, position="left", 
                  ncol=ncol.legend[length(ncol.legend)], 
                  with.missing=with.missing.legend)
      }
      popViewport(2)
    }
    popViewport()
    # Combined legends  
  }else if(withlegend=="right" || withlegend=="bottom"){
    ltext <- NULL
    cpal <- NULL
    if(plots=="both" || plots=="obs"){
      if(nchannels>1){
        for(i in 1:nchannels){
          ltext <- c(ltext, attr(obs[[i]], "labels"))
          cpal <- c(cpal, attr(obs[[i]], "cpal")) 
        }
      }else{
        ltext <- c(ltext, attr(obs, "labels"))
        cpal <- c(cpal, attr(obs, "cpal")) 
      }
    }
    if(plots=="both" || plots=="mpp"){
      ltext <- c(ltext, attr(mpp.seq, "labels"))
      cpal <- c(cpal, attr(mpp.seq, "cpal")) 
    }
    anymissing <- FALSE
    if(nchannels>1){
      for(i in 1:nchannels){
        if(any(obs[[i]]=="*")){
          anymissing <- TRUE
          break()
        }
      }
    }else{
      if(any(obs=="*")){
        anymissing <- TRUE
      }
    }
    upViewport()
    downViewport("vplegend")
    par(plt=gridPLT(), new=TRUE)
    pushViewport(viewport(width=unit(0.9, "npc")))
    if(plots=="both" || plots=="obs"){
      if(nchannels>1){
        if(withlegend=="right"){
          seqlegend(obs[[1]], fontsize=cex.legend, position="left", 
                    ncol=ncol.legend, cpal=cpal, ltext=ltext,
                    with.missing=anymissing, 
                    missing.color=attr(obs[[1]],"missing.color"))
        }else{ # withlegend=="bottom"
          seqlegend(obs[[1]], fontsize=cex.legend, position="top", 
                    ncol=ncol.legend, cpal=cpal, ltext=ltext,
                    with.missing=anymissing, 
                    missing.color=attr(obs[[1]],"missing.color"))
        }
      }else{
        if(withlegend=="right"){
          seqlegend(obs, fontsize=cex.legend, position="left", 
                    ncol=ncol.legend, cpal=cpal, ltext=ltext,
                    with.missing=anymissing, 
                    missing.color=attr(obs,"missing.color"))
        }else{ # withlegend=="bottom"
          seqlegend(obs, fontsize=cex.legend, position="top", 
                    ncol=ncol.legend, cpal=cpal, ltext=ltext,
                    with.missing=anymissing, 
                    missing.color=attr(obs,"missing.color"))
        }
      }
    }else{
      if(nchannels>1){
        if(withlegend=="right"){
          seqlegend(mpp, fontsize=cex.legend, position="left", 
                    ncol=ncol.legend, cpal=cpal, ltext=ltext,
                    with.missing=with.missing.legend, 
                    missing.color=attr(obs[[1]],"missing.color"))
        }else{ # withlegend=="bottom"
          seqlegend(mpp, fontsize=cex.legend, position="top", 
                    ncol=ncol.legend, cpal=cpal, ltext=ltext,
                    with.missing=with.missing.legend, 
                    missing.color=attr(obs[[1]],"missing.color"))
        }
      }else{
        if(withlegend=="right"){
          seqlegend(mpp, fontsize=cex.legend, position="left", 
                    ncol=ncol.legend, cpal=cpal, ltext=ltext,
                    with.missing=with.missing.legend, 
                    missing.color=attr(obs,"missing.color"))
        }else{ # withlegend=="bottom"
          seqlegend(mpp, fontsize=cex.legend, position="top", 
                    ncol=ncol.legend, cpal=cpal, ltext=ltext,
                    with.missing=with.missing.legend, 
                    missing.color=attr(obs,"missing.color"))
        }
      }
    }
    popViewport()
  }
  
  # Title
  if(!is.na(title)){
    if(title!=FALSE){
      upViewport()
      downViewport("vptitle")
      if(title.n==TRUE){
        if(nchannels>1){
          grid.text(paste0(title,", n=",dim(obs[[1]])[1]), y = unit(title.pos, "lines"), 
                    gp = gpar(cex = cex.title))
        }else{
          grid.text(paste0(title,", n=",dim(obs)[1]), y = unit(title.pos, "lines"), 
                    gp = gpar(cex = cex.title))  
        }
        popViewport()
      }else{
        grid.text(title, y = unit(title.pos, "lines"), 
                  gp = gpar(cex = cex.title))
        popViewport()
      }
    }
  }else if(title.n==TRUE){
    upViewport()
    downViewport("vptitle")
    if(nchannels>1){
      grid.text(paste0("n=",dim(obs[[1]])[1]), y = unit(title.pos, "lines"), 
                gp = gpar(cex = cex.title))
    }else{
      grid.text(paste0("n=",dim(obs)[1]), y = unit(title.pos, "lines"), 
                gp = gpar(cex = cex.title))
    }
    popViewport()
  }
  popViewport()
}