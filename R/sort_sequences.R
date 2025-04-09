#' Sort sequences in a sequence object
#' 
#' @param x An `stslist` object or a list of of such objects of same size, 
#' typically created with [TraMineR::seqdef()] or [data_to_stslist()].
#' @param sort_by A character string specifying the sorting criterion. Options
#' are `"none"` (no sorting), `"start"` (sort by the first state), `"end"` (sort by 
#' last state), and `"mds"` (sort by the multidimensional scaling).
#' @param sort_channel An integer or character string specifying the channel to 
#' sort by (unless `sort_by = "mds` in which case all channels are used for 
#' defining the sorting).
#' @param dist_method A character string specifying the distance method to use 
#' when sorting by the multidimensional scaling. Passed to 
#' [TraMineR::seqdist()], or [TraMineR::seqMD()] in the multichannel case.
#' @export
sort_sequences <- function(
    x, sort_by = "start", sort_channel = 1, dist_method = "OM") {
  
  if (identical(sort_by, "none")) return(x)
  n_channels <- if(inherits(x, "stslist")) 1 else length(x)
  if (n_channels == 1) {
    n <- nrow(x)
  } else {
    n <- nrow(x[[1]])
  }
  sl <- length(sort_by)
  sort_by <- try(
    match.arg(sort_by, c("none", "start", "end", "mds")), 
    silent = TRUE
  )
  stopifnot_(
    !inherits(sort_by, "try-error") || sl == n,
    "Argument {.arg sort_by} must be {.val none}, {.val start}, 
    {.val end}, {.val mds}, or an integer vector of length {n}."
  )

  if (sort_by %in% c("start", "end")) {
    if (n_channels == 1) {
      if (sort_by == "start") {
        ordering <- do.call(order, x)
      } else {
        ordering <- do.call(order, x[, ncol(x):1])
      }
    } else {
      if (sort_by == "start") {
        ordering <- do.call(order, x[[sort_channel]])
      } else {
        ordering <- do.call(order, x[[sort_channel]][, ncol(x[[1]]):1])
      }
    }
  } else {
    # Sort sequences according to multidimensional scaling score
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
    ordering <- order(drop(stats::cmdscale(distances, k = 1)))
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
