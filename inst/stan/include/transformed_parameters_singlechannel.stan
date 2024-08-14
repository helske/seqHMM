matrix[S, K_i] beta_i = append_row(zeros_row_K_i, beta_i_raw);
array[S] matrix[S, K_s] beta_s;
array[S] matrix[M, K_o] beta_o;

for(s in 1:S) {
  beta_s[s] = append_row(zeros_row_K_s, beta_s_raw[s]);
  beta_o[s] = append_row(zeros_row_K_o, beta_o_raw[s]);
}
