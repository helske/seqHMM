test_that("'stacked_sequence_plot' works", {
  data("leaves")
  
  expect_error(
    stacked_sequence_plot(leaves, plots = "abc"),
    "Argument `plots` must be \"obs\", \"hidden_paths\", or \"both\"."
  )
  expect_error(
    stacked_sequence_plot(leaves),
    "`x` is not an <stslist> object, a list of such objects, or an object of class <hmm>, <mhmm>, <nhmm>, or <mnhmm>."
  )
  expect_error(
    stacked_sequence_plot(hmm_biofam$observations, type = "abc"),
    "Argument `type` must be \"distribution\" or \"index\"."
  )
  expect_error(
    stacked_sequence_plot(hmm_biofam$observations, ids = -1:10),
    "Argument `ids` should be a vector of integers between 1 and 2000."
  )
  expect_error(
    stacked_sequence_plot(hmm_biofam$observations, type = "i", sort_by = "abc"),
    "Argument `sort_by` must be \"none\", \"start\", \"end\", \"mds\", or an integer vector of length 2000."
  )
  expect_error(
    stacked_sequence_plot(
      hmm_biofam$observations, ids = 1:10, sort_by = 10:1, type = "index"
    ),
    NA
  )
  expect_error(
    stacked_sequence_plot(
      hmm_biofam$observations, ids = 1:10, sort_by = "end", type = "d"
    ),
    NA
  )
  expect_error(
    stacked_sequence_plot(
      hmm_biofam$observations, ids = 1:10, sort_by = "mds", type = "d"
    ),
    NA
  )
})
