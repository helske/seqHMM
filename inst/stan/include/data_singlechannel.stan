  int<lower=1> N; // number of individuals
  int<lower=2> max_T; // (maximum) number of time points
  array[N] int<lower=2,upper=max_T> T; // number of time points per individual
  int<lower=2> M; // number of observed symbols
  int<lower=1> S; // number of hidden states
  array[max_T, N] int<lower=0, upper=M + 1> obs; // observations
