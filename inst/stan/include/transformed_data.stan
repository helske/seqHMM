  int compute_ploglik_N = 0;
  if (N != N_sample) {
    compute_ploglik_N = 1;
  } else {
    for (i in 1:N) {
      if (ids[i] != i) {
        compute_ploglik_N = 1;
        break;
      }
    }
  }

