  // extract z rv's and get on desired scale
  // fma => hp_sd * z + hp_mu
  // https://mc-stan.org/docs/stan-users-guide/reparameterizations.html
  real H = zH * y_max;
  real m = fma(x_pos_range, zm, x_pos_min);
  // this can sometimes reach zero...
  // a better solution might be to sample zs from half-normal
  // then use expm1? 
  // real s = exp(fma(hp_s[2], zs, hp_s[1])) * x_pos_range;
  real s = exp_hp_s1 * expm1(hp_s[2] * zs) * x_pos_range;
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
    p = fma(1.95, zp[1], + 0.05);
  } else {
    p = fix_p;
  }
  vector[N] mu = eff .* mskt((x - m) / s, H, r, d, p);
  vector[N] lmu = log(mu);
  // for (n in 1:N) {
  //   if (is_nan(lmu[n])) {
  //     print("eff:", eff[n],
  //           ",x:", x[n],
  //           ",H:", H, 
  //           ",m:", m,
  //           ",s:", s, 
  //           ",r:", r, 
  //           ",d:", d, 
  //           ",p:", p);
  //   }
  // }
