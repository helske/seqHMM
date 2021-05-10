.onAttach <- function(libname, pkgname) {
  note <- "Please cite seqHMM in publications by using: \n
  Helske S, Helske J (2019). Mixture Hidden Markov Models for Sequence Data: The seqHMM Package in R. Journal of
Statistical Software, 88(3), 1-32. doi: 10.18637/jss.v088.i03."
  packageStartupMessage(note)
}
