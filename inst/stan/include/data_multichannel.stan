  int<lower=1> N; // number of individuals
  int<lower=2> max_T; // (maximum) number of time points
  array[N] int<lower=2,upper=max_T> T; // number of time points per individual
  int<lower=2> max_M; // (maximum) number of observed symbols
  int<lower=1> S; // number of hidden states
  int<lower=2> C; // number of response channels
  array[C] int<lower=2, upper=max_M> M; // number of observed symbols per channel
  array[C, max_T, N] int<lower=0,upper=max_M + 1> obs; // observations
