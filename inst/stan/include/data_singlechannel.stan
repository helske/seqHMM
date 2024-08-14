int<lower=1> N; // number of individuals
int<lower=1> T; // number of time points
int<lower=2> M; // number of observed symbols
int<lower=1> S; // number of hidden states
array[T, N] int<lower=0, upper=M> obs; // observations
