  // extract z rv's and get on desired scale
  // fma => hp_sd * z + hp_mu
  // https://mc-stan.org/docs/stan-users-guide/reparameterizations.html
  real H = zH * y_max;
  real m = fma(x_pos_range, zm, x_pos_min);
  real s = exp(fma(hp_s[2], zs, hp_s[1])) * x_pos_range;

  real r;
  if (use_r > 0) {
    r = zr[1];
  } else {
    r = fix_r;
  }

  real d;
  if (use_d > 0) {
    d = exp(fma(hp_d[2], zd[1], hp_d[1]));
  } else {
    d = fix_d;
  }

  real p;
  if (use_p > 0) {
    p = fma(2.95, zp[1], + 0.05);
  } else {
    p = fix_p;
  }
