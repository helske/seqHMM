  int<lower=1> N; // number of individuals
  int<lower=1> T; // number of time points
  int<lower=2> M; // (maximum) number of observed symbols
  int<lower=1> S; // number of hidden states
  int<lower=2> C; // number of response channels
  array[C] int<lower=2, upper=M> n_M; // number of observed symbols per channel
  array[C, T, N] int<lower=0,upper=M + 1> obs; // observations
