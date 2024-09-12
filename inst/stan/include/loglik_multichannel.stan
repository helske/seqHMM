
real loglik_mc(matrix beta_i_raw, array[] matrix beta_s_raw, 
row_vector beta_o_raw, array[,,] int obs, array[] int M,
array[] int T, int N_sample, array[] int ids, array[,] vector X_s,
array[,] vector X_o, array[] vector X_i) {
  
  int S = rows(beta_i_raw) + 1;
  int N = size(X_i);
  int K_i = num_elements(X_i[1]);
  int K_s = num_elements(X_s[1, 1]);
  int K_o = num_elements(X_o[1, 1]);
  int C = size(M);
  int max_M = max(M);
  
  row_vector[K_i] zeros_row_K_i = rep_row_vector(0, K_i);
  row_vector[K_s] zeros_row_K_s = rep_row_vector(0, K_s);
  matrix[max_M, K_o] zeros_mat_M_K_o = rep_matrix(0, max_M, K_o);  
  vector[S] zeros_S = rep_vector(0, S);
  array[C + 1] int cumsum_M = append_array({0}, cumulative_sum(M));
  
  vector[N] ll;
  vector[S] log_Pi;
  matrix[S, S] log_A;
  matrix[S, max_M + 1] log_B;
  vector[S] log_py;
  vector[S] log_alpha;
  vector[S] log_alpha_new;
  matrix[S, K_i] beta_i;
  array[S] matrix[S, K_s] beta_s;
  array[C, S] matrix[max_M, K_o] beta_o;
  beta_i = append_row(zeros_row_K_i, beta_i_raw);
  for(s in 1:S) {
    beta_s[s] = append_row(zeros_row_K_s, beta_s_raw[s]);
  }
  int idx = 1;
  for (c in 1:C) {
    for (s in 1:S) {
      beta_o[c, s, , ] = zeros_mat_M_K_o;
      for (m in 2:M[c]) {
        beta_o[c, s, m, ] = beta_o_raw[idx:(idx + K_o - 1)];
        idx += K_o;
      }
    }
  }
  for(ii in 1:N_sample) {
    int i = ids[ii];
    log_Pi = log_softmax(beta_i * X_i[i]);
    log_py = zeros_S;
    for(c in 1:C) {
      for(s in 1:S) {
        log_B[s, 1:M[c]] = log_softmax(beta_o[c, s, 1:M[c], ] * X_o[1, i])';
        log_B[s, M[c] + 1] = 0;
      }
      log_py += log_B[, obs[c, 1, i]];
    }
    log_alpha = log_Pi + log_py;
    for (t in 2:T[i]) {
      log_py = zeros_S;
      for(c in 1:C) {
        for(s in 1:S) {
          log_B[s, 1:M[c]] = log_softmax(beta_o[c, s, 1:M[c], ] * X_o[t, i])';
          log_B[s, M[c] + 1] = 0;
        }
        log_py += log_B[, obs[c, t, i]];
      }
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
