// marginal likelihood via forward algorithm, with time-varying transition matrix
// log_A is T x S x S, logarithm of A
// log_py is S x T matrix of log-pmf values
real loglik(vector log_Pi, array[] matrix log_A, matrix log_py) {
  int S = rows(log_py);
  int T = cols(log_py);
  vector[S] log_alpha;
  vector[S] log_alpha_new;
  log_alpha = log_Pi + log_py[, 1];
  for (t in 2:T) {
    for (k in 1:S) {
      log_alpha_new[k] = log_sum_exp(log_alpha + log_A[t, , k] + log_py[k, t]);
    }
    log_alpha = log_alpha_new;
  }
  return log_sum_exp(log_alpha);
}

real loglik(vector log_Pi, matrix log_A, matrix log_py) {
  int S = rows(log_py);
  int T = cols(log_py);
  vector[S] log_alpha;
  vector[S] log_alpha_new;
  log_alpha = log_Pi + log_py[, 1];
  for (t in 2:T) {
    for (k in 1:S) {
      log_alpha_new[k] = log_sum_exp(log_alpha + log_A[, k] + log_py[k, t]);
    }
    log_alpha = log_alpha_new;
  }
  return log_sum_exp(log_alpha);
}
