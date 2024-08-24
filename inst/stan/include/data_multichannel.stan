  int<lower=1> N; // number of individuals
  int<lower=1> T; // number of time points
  int<lower=2> max_M; // (maximum) number of observed symbols
  int<lower=1> S; // number of hidden states
  int<lower=2> C; // number of response channels
  array[C] int<lower=2, upper=max_M> M; // number of observed symbols per channel
  array[C, T, N] int<lower=0,upper=max_M + 1> obs; // observations
