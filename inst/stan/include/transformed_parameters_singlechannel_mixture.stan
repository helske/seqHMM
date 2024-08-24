  array[D] matrix[S, K_i] beta_i;
  array[D, S] matrix[S, K_s] beta_s;
  array[D, S] matrix[M, K_o] beta_o;
  matrix[D, K_d] theta = append_row(zeros_row_K_d, theta_raw);
  for (d in 1:D) {
    beta_i[d] = append_row(zeros_row_K_i, beta_i_raw[d]);
    for (s in 1:S) {
      beta_s[d, s] = append_row(zeros_row_K_s, beta_s_raw[d, s]);
      beta_o[d, s] = append_row(zeros_row_K_o, beta_o_raw[d, s]);
    }
  }
