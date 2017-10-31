SSPlotter <- function(obs, nchannels, nplots,
  legend.c.prop, legend.r.prop, ylab.space, xaxis.space, xt.space,
  orderv, plotxaxis,
  hidden.paths, plots, type, n.seq, sortv, sort.channel,
  with.missing, missing.color, title, title.n, cex.title, title.pos,
  with.legend, ncol.legend, with.missing.legend,
  legend.prop, cex.legend, hidden.states.colors, hidden.states.labels,
  xaxis, xlab, xtlab, xlab.pos, yaxis, ylab,
  hidden.states.title, ylab.pos,
  cex.lab, cex.axis, new, call = match.call(), ...){

  if (xlab.pos * cex.lab + xaxis.space + xt.space > 0) {
    xspace <- xlab.pos * cex.lab + xaxis.space + xt.space + 0.5
  } else {
    xspace <- 0
  }

  # Grid for plotting regions
  top.vp <- viewport(layout = grid.layout(nrow = 4, ncol = 3,
    widths = unit(c(ylab.space, 1, legend.c.prop),
      c("lines", "null", "npc")),
    heights = unit(c((title.pos*cex.title+1.5), 1,
      xspace,
      legend.r.prop),
      c("lines", "null", "lines", "npc"))))

  vptitle <- viewport(layout.pos.row = 1, layout.pos.col = 2, name = "vptitle")
  vpylab <- viewport(layout.pos.row = 2, layout.pos.col = 1, name = "vpylab")
  vpplot <- viewport(layout.pos.row = 2, layout.pos.col = 2, name = "vpplot")
  if (with.legend == "right" || with.legend == "right.combined") {
    vplegend <- viewport(layout.pos.row = 2, layout.pos.col = 3, name = "vplegend")
  } else if (with.legend == "bottom" || with.legend == "bottom.combined") {
    vplegend <- viewport(layout.pos.row = 4, layout.pos.col = 2, name = "vplegend")
  }

  vpxaxis <- viewport(layout.pos.row = 3, layout.pos.col = 2, name = "vpxaxis")

  if (with.legend == FALSE) {
    splot <- vpTree(top.vp, vpList(vptitle, vpylab, vpplot, vpxaxis))
  } else {
    splot <- vpTree(top.vp, vpList(vptitle, vpylab, vpplot, vplegend, vpxaxis))
  }

  pushViewport(splot)

  # Grid for plots
  upViewport()
  downViewport("vpplot")
  pushViewport(viewport(layout = grid.layout(nrow = nplots, ncol = 1),
    width = unit(0.95, "npc")))

  # Plotting observations
  if (plots == "both" || plots == "obs") {
    if (type == "I") {
      if (is.null(sortv)) {
        if (nchannels > 1) {
          for(i in 1:(nchannels-1)) {
            pushViewport(viewport(layout.pos.col = 1, layout.pos.row = i))
            par(plt = gridPLT(), new = TRUE)
            seqplot(obs[[i]], type = "I", with.legend = FALSE,
              use.layout = FALSE, yaxis = yaxis, axes = FALSE, ylab = NA,
              xtlab = xtlab, missing.color = missing.color, ...)
            popViewport()
          }
          pushViewport(viewport(layout.pos.col = 1, layout.pos.row = nchannels))
          par(plt = gridPLT(), new = TRUE)
          seqplot(obs[[nchannels]], type = "I", with.legend = FALSE,
            use.layout = FALSE, yaxis = yaxis, xaxis = plotxaxis, ylab = NA,
            xtlab = xtlab, cex.axis = cex.axis, missing.color = missing.color, ...)
          popViewport()
        } else {
          pushViewport(viewport(layout.pos.col = 1, layout.pos.row = nchannels))
          par(plt = gridPLT(), new = TRUE)
          seqplot(obs, type = "I", with.legend = FALSE,
            use.layout = FALSE, yaxis = yaxis, xaxis = plotxaxis, ylab = NA,
            xtlab = xtlab, cex.axis = cex.axis, missing.color = missing.color, ...)
          popViewport()
        }

      } else if (length(sortv) == 1 && (sortv == "from.start" || sortv == "from.end")) {
        if (nchannels > 1) {
          for (i in 1:(nchannels-1)) {
            pushViewport(viewport(layout.pos.col = 1, layout.pos.row = i))
            par(plt = gridPLT(), new = TRUE)
            seqplot(obs[[i]][orderv,], type = "I", with.legend = FALSE,
              use.layout = FALSE, yaxis = yaxis, axes = FALSE, ylab = NA,
              xtlab = xtlab, missing.color = missing.color, ...)
            popViewport()
          }
          pushViewport(viewport(layout.pos.col = 1, layout.pos.row = nchannels))
          par(plt = gridPLT(), new = TRUE)
          seqplot(obs[[nchannels]][orderv,], type = "I", with.legend = FALSE,
            use.layout = FALSE, yaxis = yaxis, xaxis = plotxaxis, ylab = NA,
            xtlab = xtlab, cex.axis = cex.axis, missing.color = missing.color, ...)
          popViewport()
        } else {
          pushViewport(viewport(layout.pos.col = 1, layout.pos.row = nchannels))
          par(plt = gridPLT(), new = TRUE)
          seqplot(obs[orderv,], type = "I", with.legend = FALSE,
            use.layout = FALSE, yaxis = yaxis, xaxis = plotxaxis, ylab = NA,
            xtlab = xtlab, cex.axis = cex.axis, missing.color = missing.color, ...)
          popViewport()
        }


      } else if(length(sortv) > 1) {
        if (nchannels > 1) {
          for (i in 1:(nchannels-1)) {
            pushViewport(viewport(layout.pos.col = 1, layout.pos.row = i))
            par(plt = gridPLT(), new = TRUE)
            seqplot(obs[[i]], type = "I", sortv = sortv, with.legend = FALSE,
              use.layout = FALSE, yaxis = yaxis, axes = FALSE, ylab = NA,
              xtlab = xtlab, missing.color = missing.color, ...)
            popViewport()
          }
          pushViewport(viewport(layout.pos.col = 1, layout.pos.row = nchannels))
          par(plt = gridPLT(), new = TRUE)
          seqplot(obs[[nchannels]], type = "I", sortv = sortv, with.legend = FALSE,
            use.layout = FALSE, yaxis = yaxis, xaxis = plotxaxis, ylab = NA,
            xtlab = xtlab, cex.axis = cex.axis, missing.color = missing.color, ...)
          popViewport()
        } else {
          pushViewport(viewport(layout.pos.col = 1, layout.pos.row = nchannels))
          par(plt = gridPLT(), new = TRUE)
          seqplot(obs, type = "I", sortv = sortv, with.legend = FALSE,
            use.layout = FALSE, yaxis = yaxis, xaxis = plotxaxis, ylab = NA,
            xtlab = xtlab, cex.axis = cex.axis, missing.color = missing.color, ...)
          popViewport()
        }

      }
      if(plots == "obs"){
        # Close viewport "vpplot"
        popViewport(2)
      }
    } else{  # if type == "d"
      if (nchannels > 1) {
        for (i in 1:(nchannels-1)) {
          pushViewport(viewport(layout.pos.col = 1, layout.pos.row = i))
          par(plt = gridPLT(), new = TRUE)
          seqplot(obs[[i]], type = "d", with.legend = FALSE,
            use.layout = FALSE, yaxis = yaxis, axes = FALSE, ylab = NA,
            with.missing = with.missing, xtlab = xtlab, missing.color = missing.color, ...)
          popViewport()
        }
        pushViewport(viewport(layout.pos.col = 1, layout.pos.row = nchannels))
        par(plt = gridPLT(), new = TRUE)
        seqplot(obs[[nchannels]], type = "d", with.legend = FALSE, xaxis = plotxaxis,
          use.layout = FALSE, yaxis = yaxis, xaxis = plotxaxis, ylab = NA,
          with.missing = with.missing, xtlab = xtlab, cex.axis = cex.axis, 
          missing.color = missing.color, ...)
        popViewport()
      } else {
        pushViewport(viewport(layout.pos.col = 1, layout.pos.row = nchannels))
        par(plt = gridPLT(), new = TRUE)
        seqplot(obs, type = "d", with.legend = FALSE, xaxis = plotxaxis,
          use.layout = FALSE, yaxis = yaxis, xaxis = plotxaxis, ylab = NA,
          with.missing = with.missing, xtlab = xtlab, cex.axis = cex.axis, 
          missing.color = missing.color, ...)
        popViewport()
      }

      if(plots == "obs"){
        # Close viewport "vpplot"
        popViewport(2)
      }
    }
  }

  # Plotting the most probable paths
  if (plots == "both" || plots == "hidden.paths") {
    if (type == "I") {
      if (is.null(sortv)) {
        pushViewport(viewport(layout.pos.col = 1, layout.pos.row = nplots))
        par(plt = gridPLT(), new = TRUE)
        seqplot(hidden.paths, type = type, sortv = sortv, with.legend = FALSE,
          use.layout = FALSE, yaxis = yaxis, axes = xaxis, ylab = NA,
          xtlab = xtlab, cex.axis = cex.axis, missing.color = missing.color, ...)
        popViewport()
      } else if (length(sortv) == 1 && (sortv == "from.start" || sortv == "from.end")) {
        pushViewport(viewport(layout.pos.col = 1, layout.pos.row = nplots))
        par(plt = gridPLT(), new = TRUE)
        seqplot(hidden.paths[orderv,], type = type, with.legend = FALSE,
          use.layout = FALSE, yaxis = yaxis, axes = xaxis, ylab = NA,
          xtlab = xtlab, cex.axis = cex.axis, missing.color = missing.color, ...)
        popViewport()
      } else if(length(sortv) > 1) {
        pushViewport(viewport(layout.pos.col = 1, layout.pos.row = nplots))
        par(plt = gridPLT(), new = TRUE)
        seqplot(hidden.paths, type = type, sortv = sortv, with.legend = FALSE,
          use.layout = FALSE, yaxis = yaxis, axes = xaxis, ylab = NA,
          xtlab = xtlab, cex.axis = cex.axis, missing.color = missing.color, ...)
        popViewport()
      }
    } else {
      pushViewport(viewport(layout.pos.col = 1, layout.pos.row = nplots))
      par(plt = gridPLT(), new = TRUE)
      seqplot(hidden.paths, type = type, with.legend = FALSE,
        use.layout = FALSE, yaxis = yaxis, axes = xaxis, ylab = NA,
        xtlab = xtlab, cex.axis = cex.axis, missing.color = missing.color, ...)
      popViewport()
    }
    # Close viewport "vpplot"
    popViewport(2)
  }

  # Plotting x label
  if (!is.na(xlab)) {
    downViewport("vpxaxis")
    grid.text(xlab, y = unit(1, "lines"),
      gp = gpar(cex = cex.lab))
    popViewport()
  }

  # Plotting y labels (channels and/or hidden states)
  if (length(ylab) > 1 || (length(ylab) == 1 && !is.na(ylab) && ylab != FALSE)) {
    downViewport("vpylab")
    pushViewport(viewport(layout = grid.layout(nrow = nplots, ncol = 1)))
    if(plots == "both" || plots == "obs"){
      for(i in 1:nchannels){
        pushViewport(viewport(layout.pos.row = i,
          layout.pos.col = 1))
        grid.text(ylab[i], x = unit(ylab.space/cex.lab-
            ylab.pos[i]+1, "lines"),
          gp = gpar(cex = cex.lab), rot = 90, vjust = 0.5, hjust = 0.5)
        popViewport()
      }
    }
    if (plots == "both" || plots == "hidden.paths") {
      pushViewport(viewport(layout.pos.row = nplots,
        layout.pos.col = 1))
      grid.text(hidden.states.title, x = unit(ylab.space/cex.lab-
          ylab.pos[nplots]+1,
        "lines"),
        gp = gpar(cex = cex.lab), rot = 90, vjust = 0.5, hjust = 0.5)
      popViewport()
    }
    # Close viewport "vpylab"
    popViewport(2)
  }

  # Legends
  if (with.legend == "right" || with.legend == "bottom") {
    # Grid for legends
    if (with.legend == "right") {
      upViewport()
      downViewport("vplegend")
      pushViewport(viewport(layout = grid.layout(nrow = nplots, ncol = 1)))
      lposrow <- 1:nplots
      lposcol <- rep(1, nplots)
    } else {
      upViewport()
      downViewport("vplegend")
      pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = nplots)))
      lposrow <- rep(1, nplots)
      lposcol <- 1:nplots
    }

    # Legends for observations
    if (plots == "both" || plots == "obs") {
      for (i in 1:nchannels) {
        pushViewport(viewport(layout.pos.row = lposrow[i],
          layout.pos.col = lposcol[i]))
        pushViewport(viewport(width = unit(0.9, "npc")))
        par(plt = gridPLT(), new = TRUE)
        if (nchannels > 1) {
          if (with.legend == "bottom") {
            seqlegend(obs[[i]], cex = cex.legend, position = "top",
              ncol = ncol.legend[i], with.missing = with.missing.legend, 
              missing.color = missing.color)
          } else {
            seqlegend(obs[[i]], cex = cex.legend, position = "left",
              ncol = ncol.legend[i], with.missing = with.missing.legend, 
              missing.color = missing.color)
          }
        } else {
          if(with.legend == "bottom"){
            seqlegend(obs, cex = cex.legend, position = "top",
              ncol = ncol.legend[i], with.missing = with.missing.legend, 
              missing.color = missing.color)
          } else {
            seqlegend(obs, cex = cex.legend, position = "left",
              ncol = ncol.legend[i], with.missing = with.missing.legend, 
              missing.color = missing.color)
          }
        }
        popViewport(2)
      }
    }
    # Legends for most probable paths
    if (plots == "both" || plots == "hidden.paths") {
      pushViewport(viewport(layout.pos.row = lposrow[nplots],
        layout.pos.col = lposcol[nplots]))
      pushViewport(viewport(width = unit(0.9, "npc")))
      par(plt = gridPLT(), new = TRUE)
      if (with.legend == "bottom") {
        seqlegend(hidden.paths, cex = cex.legend, position = "top",
          ncol = ncol.legend[length(ncol.legend)], with.missing = with.missing.legend, 
          missing.color = missing.color, 
          ltext = hidden.states.labels)
      } else {
        seqlegend(hidden.paths, cex = cex.legend, position = "left",
          ncol = ncol.legend[length(ncol.legend)],
          with.missing = with.missing.legend, 
          missing.color = missing.color, 
          ltext = hidden.states.labels)
      }
      popViewport(2)
    }
    popViewport(2)
    # Combined legends
  } else if (with.legend == "right.combined" || with.legend == "bottom.combined") {
    ltext <- NULL
    cpal <- NULL
    if (plots == "both" || plots == "obs") {
      if (nchannels > 1) {
        for(i in 1:nchannels){
          ltext <- c(ltext, attr(obs[[i]], "labels"))
          cpal <- c(cpal, attr(obs[[i]], "cpal"))
        }
      } else {
        ltext <- c(ltext, attr(obs, "labels"))
        cpal <- c(cpal, attr(obs, "cpal"))
      }
    }
    if (plots == "both" || plots == "hidden.paths") {
      if(length(hidden.states.labels) == 1 && hidden.states.labels == "auto"){
        ltext <- c(ltext, attr(hidden.paths, "labels"))
      } else {
        ltext <- c(ltext, hidden.states.labels)
      }
      if(length(hidden.states.colors) == 1 && hidden.states.colors == "auto"){
        cpal <- c(cpal, attr(hidden.paths, "cpal"))
      } else {
        cpal <- c(cpal, hidden.states.colors)
      }
    }
    anymissing <- FALSE
    if (nchannels > 1) {
      for (i in 1:nchannels) {
        if (any(obs[[i]] == "*")) {
          anymissing <- TRUE
          break()
        }
      }
    } else {
      if (any(obs == "*")) {
        anymissing <- TRUE
      }
    }
    upViewport()
    downViewport("vplegend")
    par(plt = gridPLT(), new = TRUE)
    pushViewport(viewport(width = unit(0.9, "npc")))
    if (plots == "both" || plots == "obs") {
      if (nchannels > 1) {
        if (with.legend == "right.combined") {
          seqlegend(obs[[1]], cex = cex.legend, position = "left",
            ncol = ncol.legend, cpal = cpal, ltext = ltext,
            with.missing = anymissing,
            missing.color = ifelse(is.null(missing.color), 
                                   attr(obs[[1]],"missing.color"), missing.color))
        } else { # with.legend == "bottom.combined"
          seqlegend(obs[[1]], cex = cex.legend, position = "top",
            ncol = ncol.legend, cpal = cpal, ltext = ltext,
            with.missing = anymissing,
            missing.color = ifelse(is.null(missing.color), 
                                   attr(obs[[1]],"missing.color"), missing.color))
        }
      } else {
        if (with.legend == "right.combined") {
          seqlegend(obs, cex = cex.legend, position = "left",
            ncol = ncol.legend, cpal = cpal, ltext = ltext,
            with.missing = anymissing,
            missing.color = ifelse(is.null(missing.color), 
                                   attr(obs,"missing.color"), missing.color))
        } else { # with.legend == "bottom.combined"
          seqlegend(obs, cex = cex.legend, position = "top",
            ncol = ncol.legend, cpal = cpal, ltext = ltext,
            with.missing = anymissing,
            missing.color = ifelse(is.null(missing.color), 
                                   attr(obs,"missing.color"), missing.color))
        }
      }
    } else {
      if (nchannels > 1) {
        if (with.legend == "right.combined") {
          seqlegend(hidden.paths, cex = cex.legend, position = "left",
            ncol = ncol.legend, cpal = cpal, ltext = ltext,
            with.missing = with.missing.legend,
            missing.color = ifelse(is.null(missing.color), 
                                   attr(obs[[1]],"missing.color"), missing.color))
        } else { # with.legend == "bottom.combined"
          seqlegend(hidden.paths, cex = cex.legend, position = "top",
            ncol = ncol.legend, cpal = cpal, ltext = ltext,
            with.missing = with.missing.legend,
            missing.color = ifelse(is.null(missing.color), 
                                   attr(obs[[1]],"missing.color"), missing.color))
        }
      } else {
        if(with.legend == "right.combined"){
          seqlegend(hidden.paths, cex = cex.legend, position = "left",
            ncol = ncol.legend, cpal = cpal, ltext = ltext,
            with.missing = with.missing.legend,
            missing.color = ifelse(is.null(missing.color), 
                                   attr(obs,"missing.color"), missing.color))
        } else { # with.legend == "bottom.combined"
          seqlegend(hidden.paths, cex = cex.legend, position = "top",
            ncol = ncol.legend, cpal = cpal, ltext = ltext,
            with.missing = with.missing.legend,
            missing.color = ifelse(is.null(missing.color), 
                                   attr(obs,"missing.color"), missing.color))
        }
      }
    }
    popViewport(2)
  }

  # Title
  if (!is.na(title)) {
    if (title != FALSE) {
      upViewport()
      downViewport("vptitle")
      if (title.n) {
        if (nchannels > 1) {
          grid.text(paste0(title,", n = ", n.seq), y = unit(title.pos, "lines"),
            gp = gpar(cex = cex.title))
        } else {
          grid.text(paste0(title,", n = ", n.seq), y = unit(title.pos, "lines"),
            gp = gpar(cex = cex.title))
        }
        popViewport()
      } else {
        grid.text(title, y = unit(title.pos, "lines"),
          gp = gpar(cex = cex.title))
        popViewport()
      }
    }
  } else if (title.n) {
    upViewport()
    downViewport("vptitle")
    if (nchannels > 1) {
      grid.text(paste0("n = ", n.seq), y = unit(title.pos, "lines"),
        gp = gpar(cex = cex.title))
    } else {
      grid.text(paste0("n = ", n.seq), y = unit(title.pos, "lines"),
        gp = gpar(cex = cex.title))
    }
    popViewport()
  }
  popViewport()
}
