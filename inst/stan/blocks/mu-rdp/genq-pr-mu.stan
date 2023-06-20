    // dont store these
    vector[N] pr_mu;

    // generate pars first for log_prior then transform
    real pr_zH = beta_rng(hp_H[1], hp_H[2]);
    log_prior += beta_lpdf(pr_zH | hp_H[1], hp_H[2]);
    pr_H = pr_zH * y_max;
    
    real pr_zm = beta_rng(hp_m[1], hp_m[2]);
    log_prior += beta_lpdf(pr_zm | hp_m[1], hp_m[2]);
    pr_m = fma(x_pos_range, pr_zm, x_pos_min);

    real pr_zs = std_normal_rng();
    log_prior += std_normal_lpdf(pr_zs);
    pr_s = exp(fma(hp_s[2], pr_zs, hp_s[1])) * x_pos_range;
    
    if (use_r > 0) {
      pr_r = beta_rng(hp_r[1], hp_r[2]);
      log_prior += beta_lpdf(pr_r | hp_r[1], hp_r[2]);
    } else {
      pr_r = fix_r;
    }
    
    if (use_d > 0) {
      real pr_zd = std_normal_rng();
      log_prior += std_normal_lpdf(pr_zd);
      pr_d = exp(fma(hp_d[2], pr_zd, hp_d[1]));
    } else {
      pr_d = fix_d;
    }
    
    if (use_p > 0) {
      real pr_zp = beta_rng(hp_p[1], hp_p[2]);
      log_prior += beta_lpdf(pr_zp | hp_p[1], hp_p[2]);
      pr_p = fma(2.95, pr_zp, 0.05);
    } else {
      pr_p = fix_p;
    }
    
    // maybe we can just calc this and derive the evenly gridded versions after?
    pr_mu = exp(leff + lmskt((x - pr_m) / pr_s, pr_H, pr_r, pr_d, pr_p));
    pr_mu_rep = exp(lmskt((xrep - pr_m) / pr_s, pr_H, pr_r, pr_d, pr_p));
