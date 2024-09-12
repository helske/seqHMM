
real loglik_sc_mix(array[] matrix beta_i_raw, array[,] matrix beta_s_raw, 
array[,] matrix beta_o_raw, matrix theta_raw, array[,] int obs, int M,
array[] int T, int N_sample, array[] int ids, array[,] vector X_s,
array[,] vector X_o, array[] vector X_i, array[] vector X_d) {
  
  int D = size(beta_i_raw);
  int S = rows(beta_i_raw[1]) + 1;
  int N = size(X_i);
  int K_i = num_elements(X_i[1]);
  int K_s = num_elements(X_s[1, 1]);
  int K_o = num_elements(X_o[1, 1]);
  int K_d = num_elements(X_d[1]);
  
  row_vector[K_i] zeros_row_K_i = rep_row_vector(0, K_i);
  row_vector[K_s] zeros_row_K_s = rep_row_vector(0, K_s);
  row_vector[K_o] zeros_row_K_o = rep_row_vector(0, K_o); 
  row_vector[K_d] zeros_row_K_d = rep_row_vector(0, K_d);
  vector[S] zeros_S = rep_vector(0, S);
  
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
  vector[N_sample] ll_i;
  for(ii in 1:N_sample) {
    int i = ids[ii];
    vector[D] log_omega = log_softmax(theta * X_d[i]);
    ll_i[i] = log_sum_exp(log_omega + ll[i, ]');
  }
  return sum(ll_i);
}
