#' Merge Multiple Sequence Objects into One (from Multichannel to Single Channel Data)
#'
#' Function `mc_to_sc_data` combines observed states of multiple
#'   sequence objects into one, time point by time point.
#'
#' @export
#' @param data A list of state sequence objects (`stslist`s)
#'   created with the [seqdef()] function.
#' @param combine_missing Controls whether combined states of observations
#'   at time t are coded missing (coded with * in `stslist`s)
#'   if one or more of the channels include missing information at time t.
#'   Defaults to `TRUE`. `FALSE` keeps missing states
#'   as they are, producing more states in data; e.g. single/childless/*
#'   where the observation in channel 3 is missing.
#' @param all_combinations Controls whether all possible combinations of
#'   observed states are included in the single channel representation or
#'   only combinations that are found in the data. Defaults to `FALSE`,
#'   i.e. only actual observations are included.
#' @param cpal The color palette used for the new combined symbols. Optional in
#'   a case where the number of symbols is less or equal to 200 (in which case
#'   the `seqHMM::colorpalette` is used).
#' @seealso [mc_to_sc()] for transforming multichannel `hmm`
#' or `mhmm` objects into single-channel representations;
#' [stacked_sequence_plot] for plotting multiple sequence data sets in the
#' same plot; and [seqdef()] for creating state sequence objects.
#'
#' @examples
#' # Load three-channel sequence data
#' data("biofam3c")
#'
#' # Building sequence objects
#' marr_seq <- seqdef(biofam3c$married,
#'   start = 15,
#'   alphabet = c("single", "married", "divorced"),
#'   cpal = c("violetred2", "darkgoldenrod2", "darkmagenta")
#' )
#' child_seq <- seqdef(biofam3c$children,
#'   start = 15,
#'   alphabet = c("childless", "children"),
#'   cpal = c("darkseagreen1", "coral3")
#' )
#' left_seq <- seqdef(biofam3c$left,
#'   start = 15,
#'   alphabet = c("with parents", "left home"),
#'   cpal = c("lightblue", "red3")
#' )
#'
#' # Converting multichannel data to single-channel data
#' sc_data <- mc_to_sc_data(list(marr_seq, child_seq, left_seq))
#'
#' # 10 combined states
#' alphabet(sc_data)
#'
#' # Colors for combined states
#' attr(sc_data, "cpal") <- colorpalette[[14]][1:10]
#'
#' # Plotting sequences for the first 10 subjects
#' stacked_sequence_plot(
#'   list(
#'     "Marriage" = marr_seq, "Parenthood" = child_seq,
#'     "Residence" = left_seq, "Combined" = sc_data
#'   ),
#'   type = "i",
#'   ids = 1:10
#' )
#'
#'
#' # Including all combinations (whether or not available in data)
#' sc_data_all <- mc_to_sc_data(list(marr_seq, child_seq, left_seq),
#'   all_combinations = TRUE
#' )
#'
#' # 12 combined states, 2 with no observations in data
#' seqstatf(sc_data_all)
#'
mc_to_sc_data <- function(data, combine_missing = TRUE, all_combinations = FALSE, cpal) {
  stopifnot_(
    n_unique(vapply(data, nrow, 1L)) == 1L,
    "The number of subjects (rows) is not the same in all channels."
  )
  stopifnot_(
    n_unique(vapply(data, ncol, 1L)) == 1L,
    "The length of the sequences (columns) is not the same in all channels."
  )
  alph <- apply(expand.grid(lapply(data, alphabet)), 1, paste0, collapse = "/")
  datax <- data[[1]]
  for (i in 2:length(data)) {
    datax <- as.data.frame(mapply(paste, datax,
                                  data[[i]],
                                  USE.NAMES = FALSE, SIMPLIFY = FALSE,
                                  MoreArgs = list(sep = "/")
    ))
  }
  names(datax) <- names(data[[1]])
  if (combine_missing == TRUE) {
    datax[Reduce("|", lapply(
      data,
      function(x) {
        x == attr(data[[1]], "nr") |
          x == attr(data[[1]], "void") |
          is.na(x)
      }
    ))] <- NA
  }
  if (missing(cpal) || identical(cpal, "auto")) {
    stopifnot_(
      length(alph) <= 200,
      "The number of observed states is {length(alph)}, which is more than 
      supported by the default color palette (200). Specify your own color 
      palette with the argument {.arg cpal}."
    )
    cpal <- seqHMM::colorpalette[[length(alph)]]
  }
  stopifnot_(
    length(alph) == length(cpal),
    "The number of observed states is {length(alph)} but the supplied 
    color palette contains only {length(cpal)} colours."
  )
  if (all_combinations == TRUE) {
    datax <- suppressWarnings(suppressMessages(seqdef(datax, alphabet = alph)))
  } else {
    datax <- suppressWarnings(suppressMessages((seqdef(datax))))
  }
  attr(datax, "xtstep") <- attr(data[[1]], "xtstep")
  attr(datax, "missing.color") <- attr(data[[1]], "missing.color")
  attr(datax, "nr") <- attr(data[[1]], "nr")
  attr(datax, "void") <- attr(data[[1]], "void")
  attr(datax, "missing") <- attr(data[[1]], "missing")
  attr(datax, "start") <- attr(data[[1]], "start")
  TraMineR::cpal(datax) <- cpal[alph %in% alphabet(datax)]
  
  datax
}
