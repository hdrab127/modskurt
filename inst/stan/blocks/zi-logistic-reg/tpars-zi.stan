  array[use_zi] real g0;
  array[use_zi] real g1;
  vector[N] zi;
  if (use_zi) {
    // logistic function of the relative mean
    g0[1] = fma(hp_g0[2], zg0[1], hp_g0[1]);
    g1[1] = hp_g1[1] * zg1[1];
    zi = inv_logit(fma(-g1[1] / H, exp(lmu), g0[1]));
  }
