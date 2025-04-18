#' Transform TraMineR's state sequence object to data.table and vice versa
#'
#' @param x For `data_to_stslist`, a `data.frame` type of object in long format, 
#' or a model object of class `nhmm` or `mnhmm`. 
#' For `stslist_to_data`, an object of class `stslist` or list of such objects.
#' @param time A character string specifying the time variable. Ignored if `x` 
#' is NHMM.
#' @param id A character string specifying the id variable. Ignored if `x` 
#' is NHMM.
#' @param responses A character vector specifying the name(s) of the response 
#' variable(s). Ignored if `x` is NHMM.
#' @param seqdef_args A list of additional arguments to [TraMineR::seqdef()] in 
#' case of  `data_to_stslist`. In case of `length(responses) > 1`, a list of 
#' lists. Ignored in `stslist_to_data`.
#' @param ... Ignored
#' @rdname data_to_stslist
#' @export
data_to_stslist <- function(x, id, time, responses, seqdef_args = NULL, ...) {
  
  stopifnot_(
    !missing(responses) && checkmate::test_character(x = responses), 
    "Argument {.arg responses} must be a character vector defining the response 
    variable(s) in the {.arg x}."
  )
  stopifnot_(
    length(responses) == length(unique(responses)), 
    "Response names in {.arg responses} should be unique."
  )
  
  if (inherits(x, "nhmm") || inherits(x, "mnhmm")) {
    responses <- x$responses
    x <- x$data[, c(x$id_variable, x$time_variable)]
  } else {
    cols <- c(id, time, responses)
    x <- as.data.table(x)
    x[, (responses) := lapply(.SD, as.factor), .SDcols = responses]
    x <- .check_data(x, id, time, responses)[, cols, env = list(cols = I(cols))]
    x <- fill_time(x, id, time)
  }
  sequences <- vector("list", length(responses))
  names(sequences) <- responses
  colnames(x)[1:2] <- c("id", "time")
  if (!is.null(seqdef_args)) {
    if (length(responses) > 1) {
      stopifnot_(
        is_list_of_lists(seqdef_args) && length(seqdef_args) == length(responses),
        "Argument {.arg seqdef_args} should be a list of lists in case of 
        multiple responses."
      )
    } else {
      seqdef_args <- list(seqdef_args)
    }
    names(seqdef_args) <- responses
  }
  for (y in responses) {
    wide_data <- dcast(x, id ~ time, value.var = y, drop = FALSE)
    sequences[[y]] <- do.call(
      TraMineR::seqdef, 
      c(
        list(
          data = wide_data[, -1],
          informat = "STS",
          alphabet = levels(x[[y]]),
          id = wide_data[["id"]]
        ),
        seqdef_args[[y]]
      )
    )
  }
  if (length(responses) == 1) sequences[[1]] else sequences
}
#' @rdname data_to_stslist
#' @export
stslist_to_data <- function(x, id, time, responses, ...) {
  stopifnot_(
    !missing(responses) && checkmate::test_character(x = responses), 
    "Argument {.arg responses} must be a character vector defining the names of 
    the response variable(s)."
  )
  stopifnot_(
    length(responses) == length(unique(responses)), 
    "Response names in {.arg responses} should be unique."
  )
  if (TraMineR::is.stslist(x)) {
    x <- list(x)
  } else {
    stopifnot_(
      is_list_of_lists(x) && all(unlist(lapply(x, TraMineR::is.stslist))),
      "{.arg observations} should be a {.cls stslist} object created with 
        {.fn seqdef}, or a {.cls list} of {.cls stslist} objects in a 
        multichannel case."
    )
    stopifnot_(
      length(unique(sapply(x, nrow))) == 1,
      "The number of subjects (rows) is not the same in all channels."
    )
    stopifnot_(
      length(unique(sapply(x, ncol))) == 1,
      "The length of the sequences (columns) is not the same in all channels."
    )
  }
  times <- colnames(x[[1]])
  na_times <- suppressWarnings(any(is.na(timenames <- as.numeric(times))))
  if (na_times) {
    na_times <- suppressWarnings(any(is.na(timenames <- as.numeric(sub('.', '', times)))))
    if (na_times) {
      warning_(
        paste0(
          "Time indices (column names) of sequences are not coarceable to ",
          "numeric. Replacing them with integers."
        )
      )
      timenames <- seq_len(ncol(x[[1]]))
    }
  }
  stopifnot_(
    identical(sort(timenames), timenames),
    paste0(
      "The numeric time indices based on column names of sequence object ", 
      "are not numerically sorted. Please recode the column names.")
  )
  ids <- factor(rownames(x[[1]]), levels = rownames(x[[1]]))
  data <- data.table(
    id = rep(ids, times = length(timenames)), 
    time = rep(timenames, each = length(ids))
  )
  for (i in seq_along(x)) {
    data[, (responses[i]) := factor(unlist(x[[i]]), levels = alphabet(x[[i]]))]
  }
  setkeyv(data, c("id", "time"))
  setnames(data, c("id", "time"), c(id, time))
  data[]
}

