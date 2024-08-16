
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
    stopifnot_(
      is_list_of_lists(x) && all(unlist(lapply(x, TraMineR::is.stslist))),
      "{.arg {what}} should a {.cls stslist} object created with 
      {.fn seqdef} or, in a case of multichannel data, a {.cls list} of such 
      objects."
    )
    multichannel <- TRUE
    stopifnot_(
      length(unique(sapply(x, nrow))) == 1,
      "The number of subjects (rows) is not the same in all channels."
    )
    stopifnot_(
      length(unique(sapply(x, ncol))) == 1,
      "The length of the sequences (columns) is not the same in all channels."
    )
  }
  multichannel
}
