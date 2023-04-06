
is_list_of_lists <- function(x) {
  if (!is.list(x)) {
    return(FALSE)
  } else {
    if (all(unlist(lapply(x, is.list)))) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
}

is_multichannel <- function(x, what = "observations") {
  if (TraMineR::is.stslist(x)) {
    multichannel <- FALSE
  } else {
    if (is_list_of_lists(x) && all(unlist(lapply(x, TraMineR::is.stslist)))) {
      multichannel <- TRUE
      if (length(unique(sapply(x, nrow))) > 1) {
        stop("The number of subjects (rows) is not the same in all channels.")
      }
      if (length(unique(sapply(x, ncol))) > 1) {
        stop("The length of the sequences (number of columns) is not the same in all channels.")
      }
    } else {
      stop(
        paste0(
          "Argument '", what, "' should a 'stslist' object ",
          "created with 'seqdef' function, or a list of such objects in case ",
          "of multichannel data."
        )
      )
    }
  }
  multichannel
}
