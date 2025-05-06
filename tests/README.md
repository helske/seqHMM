### Running extended tests

The testthat directory contains long-running tests in `test-extended`which are 
not run by default. These tests can be switched on by defining an 
environmental variables `SEQHMM_EXTENDED_TESTS = "true"` (e.g., 
`Sys.setenv("SEQHMM_EXTENDED_TESTS" = "true")` in R), and on GitHub Actions by 
adding `run-extended` to the commit message.
