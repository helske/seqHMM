        for(c in 1:C) {
          for(s in 1:S) {
            log_B[c, s, 1:M] = log_softmax(beta_o[c, s] * X_o[t, i])';
            log_B[c, s, M + 1] = 0;
          }
        }
