#' Sort sequences in a sequence object
#' 
#' @param x A sequence object or a list of sequence objects
#' @param sort_by A character string specifying the sorting criterion.
#' @param sort_channel An integer or character string specifying the channel to 
#' sort by.
#' @param dist_method A character string specifying the distance method to use.
#' @export
sort_sequences <- function(
    x, sort_by = "start", sort_channel = 1, dist_method = "OM") {
  
  n_channels <- if(inherits(x, "stslist")) 1 else length(x)
  if (sort_by %in% c("start", "end")) {
    if (n_channels == 1) {
      idx <- seq_len(max(seqlength(x)))
      if (sort_by == "start") {
        ordering <- do.call(order, x[, idx])
      } else {
        ordering <- do.call(order, x[, rev(idx)])
      }
    } else {
      idx <- seq_len(max(seqlength(x)))
      if (sort_by == "start") {
        ordering <- do.call(order, x[[sort_channel]][, idx])
      } else {
        ordering <- do.call(order, x[[sort_channel]][, rev(idx)])
      }
    }
  } else {
    # Sort sequences according to multidimensional scaling score of hidden.paths
    if (n_channels == 1) {
      distances <- suppressWarnings(suppressMessages(
        TraMineR::seqdist(
          x, method = dist_method, sm = "TRATE", with.missing = TRUE
        )
      ))
    } else {
      distances <- suppressWarnings(suppressMessages(
        TraMineR::seqMD(
          x, method = dist_method, sm = "TRATE", what = "diss", 
          with.missing = NULL
        )
      ))
    }
    ordering <- order(drop(cmdscale(distances, k = 1)))
  }
  if (n_channels == 1) {
    x <- x[ordering,]
  } else {
    for (i in seq_len(n_channels)) {
      x[[i]] <- x[[i]][ordering, ]
    }
  }
  x
}