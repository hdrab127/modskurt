  /*
  Log of the "mskt" unimodal mean function
  @param x predictor variable (real number vector)
  @param H height of curve at mode (real number > 0)
  @param m mode (real number)
  @param s scale (real number > 0)
  @param r asymmetry (real number in (0, 1))
  @param d flatness (real number >= 0)
  @param p tail exaggeration (real number in (0, 2))
  @return vector of mean values (real number vector > 0)
  */
  vector lmskt(vector xms, real H, real r, real d, real p) {
    real invp = 1 / p;
    real mdop = - d * invp;
    return log(fma(H,
                 pow(fma(r,
                         exp(fma(invp, xms / r, mdop)),
                         fma(1 - r,
                             exp(fma(-invp, xms / (1 - r), mdop)),
                             exp(- mdop) + 1)),
                     -p),
                 machine_precision()));
  }
  // return log(H) - lmultiply(p, fma(
  //   r,
  //   exp( (xms ./ r - d) / p ),
  //   fma(
  //     1 - r,
  //     exp( -(xms ./ (1 - r) + d) / p ),
  //     - exp( -d / p ) + 1  
  //   )
  // ));
  real std_halfnormal_rng() {
    real u = uniform_rng(0.5, 1);
    real y = inv_Phi(u);
    return y;
  }
