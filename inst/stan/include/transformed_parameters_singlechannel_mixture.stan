real log_lik;
real prior = 0;
{
  matrix[N, D] ll;
  vector[S] log_Pi;
  matrix[S, S] log_A;
  matrix[S, M + 1] log_B;
  vector[S] log_py;
  vector[S] log_alpha;
  vector[S] log_alpha_new;
  matrix[S, K_i] beta_i;
  array[S] matrix[S, K_s] beta_s;
  array[S] matrix[M, K_o] beta_o;
  matrix[D, K_d] theta = append_row(zeros_row_K_d, theta_raw);
  for(d in 1:D) {
    // priors for (very weak) regularisation
    prior += normal_lpdf(to_vector(beta_i_raw[d]) | 0, 5);
    for(s in 1:S) {
      prior += normal_lpdf(to_vector(beta_s_raw[d, s]) | 0, 5);
      prior += normal_lpdf(to_vector(beta_o_raw[d, s]) | 0, 5);
    }
    beta_i = append_row(zeros_row_K_i, beta_i_raw[d]);
    for(s in 1:S) {
      beta_s[s] = append_row(zeros_row_K_s, beta_s_raw[d, s]);
      beta_o[s] = append_row(zeros_row_K_o, beta_o_raw[d, s]);
    }
    for(ii in 1:N_sample) {
      int i = ids[ii];
      log_Pi = log_softmax(beta_i * X_i[i]);
      for(s in 1:S) {
        log_B[s, 1:M] = log_softmax(beta_o[s] * X_o[1, i])';
        log_B[s, M + 1] = 0;
      }
      log_py = log_B[, obs[1, i]];
      log_alpha = log_Pi + log_py;
      for (t in 2:T[i]) {
        for(s in 1:S) {
          log_B[s, 1:M] = log_softmax(beta_o[s] * X_o[t, i])';
          log_B[s, M + 1] = 0;
        }
        log_py = log_B[, obs[t, i]];
        for(s in 1:S) {
          log_A[s, ] = log_softmax(beta_s[s] * X_s[t - 1, i])';
        }
        for (k in 1:S) {
          log_alpha_new[k] = log_sum_exp(log_alpha + log_A[, k] + log_py[k]);
        }
        log_alpha = log_alpha_new;
      }
      ll[i, d] = log_sum_exp(log_alpha);
    }
  }
  prior += normal_lpdf(to_vector(theta_raw) | 0, 5);
  vector[N_sample] ll_i;
  for(ii in 1:N_sample) {
    int i = ids[ii];
    vector[D] log_omega = log_softmax(theta * X_d[i]);
    ll_i[i] = log_sum_exp(log_omega + ll[i, ]');
  }
  log_lik = sum(ll_i);
}
