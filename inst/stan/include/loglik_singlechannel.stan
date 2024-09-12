
real loglik_sc(matrix beta_i_raw, array[] matrix beta_s_raw, 
array[] matrix beta_o_raw, array[,] int obs, int M,
array[] int T, int N_sample, array[] int ids, array[,] vector X_s,
array[,] vector X_o, array[] vector X_i) {
  
  int S = rows(beta_i_raw) + 1;
  int N = size(X_i);
  int K_i = num_elements(X_i[1]);
  int K_s = num_elements(X_s[1, 1]);
  int K_o = num_elements(X_o[1, 1]);
  
  row_vector[K_i] zeros_row_K_i = rep_row_vector(0, K_i);
  row_vector[K_s] zeros_row_K_s = rep_row_vector(0, K_s);
  row_vector[K_o] zeros_row_K_o = rep_row_vector(0, K_o);
  vector[S] zeros_S = rep_vector(0, S);
  
  vector[N] ll;
  vector[S] log_Pi;
  matrix[S, S] log_A;
  matrix[S, M + 1] log_B;
  vector[S] log_py;
  vector[S] log_alpha;
  vector[S] log_alpha_new;
  matrix[S, K_i] beta_i;
  array[S] matrix[S, K_s] beta_s;
  array[S] matrix[M, K_o] beta_o;
  beta_i = append_row(zeros_row_K_i, beta_i_raw);
  for(s in 1:S) {
    beta_s[s] = append_row(zeros_row_K_s, beta_s_raw[s]);
    beta_o[s] = append_row(zeros_row_K_o, beta_o_raw[s]);
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
    ll[i] = log_sum_exp(log_alpha);
  }
  return sum(ll);
}
