    // calc zi & yreps
    vector[N] pr_zi;
    if (use_zi) {
      pr_g0 = beta_rng(hp_g0[1], hp_g0[2]);
      pr_g1 = hp_g1[1] * std_halfnormal_rng();
      pr_zi = inv_logit(fma(-pr_g1 / pr_H, pr_mu, logit(pr_g0)));
      pr_zi_rep = inv_logit(fma(-pr_g1 / pr_H, pr_mu_rep, logit(pr_g0)));
    }
