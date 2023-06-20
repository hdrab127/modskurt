  target += beta_lupdf(zH | hp_H[1], hp_H[2]);
  target += beta_lupdf(zm | hp_m[1], hp_m[2]);
  target += std_normal_lupdf(zs);
  if (use_r > 0) {
    target += beta_lupdf(zr[1] | hp_r[1], hp_r[2]);
  }
  if (use_d > 0) {
    target += std_normal_lupdf(zd[1]);
  }
  if (use_p > 0) {
    target += beta_lupdf(zp[1] | hp_p[1], hp_p[2]);
  }
