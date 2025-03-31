.onAttach <- function(libname, pkgname) {
  note <- "Please cite seqHMM appropriately, see `citation('seqHMM')` for details."
  packageStartupMessage(note)
}
