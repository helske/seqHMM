  matrix[S, K_i] beta_i = append_row(zeros_row_K_i, beta_i_raw);
  array[S] matrix[S, K_s] beta_s;
  array[C, S] matrix[max_M, K_o] beta_o;
  for (s in 1:S) {
    beta_s[s] = append_row(zeros_row_K_s, beta_s_raw[s]);
  }
  {
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
  }
