  /*
  ln-modskurt mean function with args ...
  */
  vector lmskt(vector zms, real H, real r, real d, real p) {
    return log(H) - lmultiply(p, fma(
      r,
      exp( (zms ./ r - d) / p ),
      fma(
        1 - r,
        exp( -(zms ./ (1 - r) + d) / p ),
        - exp( -d / p ) + 1  
      )
    ));
  }
  real std_halfnormal_rng() {
    real u = uniform_rng(0.5, 1);
    real y = inv_Phi(u);
    return y;
  }
