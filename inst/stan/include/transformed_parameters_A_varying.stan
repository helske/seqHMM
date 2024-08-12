array[S] matrix[S, K_s] beta_s;
for(s in 1:S) {
  beta_s[s] = append_row(rep_row_vector(0, K_s), beta_s_raw[s]);
}
