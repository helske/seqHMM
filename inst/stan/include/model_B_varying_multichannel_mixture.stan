        for(c in 1:C) {
          for(s in 1:S) {
            log_B[c, s, ] = zeros_M;
            log_B[c, s, 1:M[c]] = log_softmax(beta_o[d, c, s, 1:M[c], ] * X_o[t, i])';
          }
        }
