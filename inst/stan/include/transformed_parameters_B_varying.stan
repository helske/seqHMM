array[S] matrix[M, K_o] beta_o;
for(s in 1:S) {
  beta_o[s] = append_row(rep_row_vector(0, K_o), beta_o_raw[s]);
}
