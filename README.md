LifeSequenceHMM: Hidden Markov Models for Life Sequences and Other Multivariate, Multichannel Categorical Time Series
====================================================================================================

LifeSequenceHMM is an R package with focus on fast algorithms for estimation and inference of hidden 
Markov models where both the hidden state space and the symbol space of observations is discrete, 
and observations consists of multiple sequences as channels such as life calendar data with multiple 
life domains. LifeSequenceHMM supports maximum likelihood estimation of HMMs via numerical 
optimization routines nlm and bobyqa. Main functions of LifeSequenceHMM are written in Fortran for 
improved performance.

Package is still in testing phase, but will be available in CRAN later in 2014.