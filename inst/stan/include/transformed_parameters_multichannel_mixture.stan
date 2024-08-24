  array[D] matrix[S, K_i] beta_i;
  array[D, S] matrix[S, K_s] beta_s;
  array[D, C, S] matrix[max_M, K_o] beta_o;
  matrix[D, K_d] theta = append_row(zeros_row_K_d, theta_raw);
  for(d in 1:D) {
    beta_i[d] = append_row(zeros_row_K_i, beta_i_raw[d]);
    for (s in 1:S) {
      beta_s[d, s] = append_row(zeros_row_K_s, beta_s_raw[d, s]);
    }
    int idx = 1;
    for (c in 1:C) {
      for (s in 1:S) {
        beta_o[d, c, s, , ] = zeros_mat_M_K_o;
        for (m in 2:M[c]) {
          beta_o[d, c, s, m, ] = beta_o_raw[d, idx:(idx + K_o - 1)];
          idx += K_o;
        }
      }
    }
  }
