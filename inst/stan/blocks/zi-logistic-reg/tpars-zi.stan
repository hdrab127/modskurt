  array[use_zi] real g1;
  vector[N] zi;
  if (use_zi) {
    // logistic function of the relative mean
    g1[1] = hp_g1[1] * zg1[1];
    zi = inv_logit(fma(-g1[1] / H, exp(lmu), logit(g0[1])));
  }
