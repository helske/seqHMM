[![R-CMD-check](https://github.com/helske/seqHMM/workflows/R-CMD-check/badge.svg)](https://github.com/helske/seqHMM/actions)
[![Codecov test coverage](https://codecov.io/gh/helske/seqHMM/branch/main/graph/badge.svg)](https://app.codecov.io/gh/helske/seqHMM?branch=main)
[![cran version](https://www.r-pkg.org/badges/version/seqHMM)](https://CRAN.R-project.org/package=seqHMM)
[![downloads](https://cranlogs.r-pkg.org/badges/seqHMM)](https://cranlogs.r-pkg.org/badges/seqHMM)

seqHMM: Hidden Markov Models for Life Sequences and Other Multivariate, Multichannel Categorical Time Series
====================================================================================================

The `seqHMM` package is designed for fitting hidden (or latent) Markov models (HMMs) and mixture hidden Markov models (MHMMs) for social sequence data and other categorical time series. Also some more restricted versions of these type of models are now available: Markov models, mixture Markov models, and latent class models. 

The package supports models for one or multiple subjects with one or multiple parallel sequences (channels). External covariates can be added to explain cluster membership in mixture models. The package provides functions for evaluating and comparing models, as well as functions for easy plotting of multichannel sequence data and hidden Markov models.

Maximum likelihood estimation via EM algorithm and direct numerical maximization with analytical gradients is supported. All main algorithms are written in C++ with parallel computation support via OpenMP.

When using the package, please cite:

Helske, Satu and Helske, Jouni (2019). Mixture hidden Markov models for sequence data: the seqHMM package in R. *Journal of Statistical Software, 88*(3). [doi:10.18637/jss.v088.i03](https://dx.doi.org/10.18637/jss.v088.i03)

If you find bugs, please add a new issue here in GitHub. You can also contact Satu Helske (firstname.lastname@utu.fi) or Jouni Helske (firstname.lastname@iki.fi). We would be happy to hear your feedback.

The package is available on CRAN. Install it via

```R
install.packages("seqHMM")
```

If you want to try the development version of the `seqHMM` package, install it from Github using the `remotes` package:

```R
install.packages("remotes")
library("remotes")
install_github("helske/seqHMM")
```

