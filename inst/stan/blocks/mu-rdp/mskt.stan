  /*
  The "mskt" unimodal mean function
  @param xms ( predictor variable (real number vector) - mode (m, real number) ) 
              / scale (s, real number > 0)
  @param H height of curve at mode (real number > 0)
  @param r asymmetry (real number in (0, 1))
  @param d flatness (real number >= 0)
  @param p tail exaggeration (real number in (0.05, 1.95))
  @return vector of mean values (real number vector > 0)
  */
  vector mskt(vector xms, real H, real r, real d, real p) {
    real invp = 1 / p;
    real mdop = - d * invp;
    return fma(H,
               pow(fma(r, exp(fma(invp, xms / r, mdop)),
                       fma(1 - r, exp(fma(-invp, xms / (1 - r), mdop)),
                           - exp(mdop) + 1)),
                   -p),
               // avoid underflow for really small mu
               machine_precision());
    // unintuitively log scale creates more arithmetic precision issues
    // return log(H) - lmultiply(p, fma(
    //   r,
    //   exp( (xms ./ r - d) / p ),
    //   fma(
    //     1 - r,
    //     exp( -(xms ./ (1 - r) + d) / p ),
    //     - exp( -d / p ) + 1  
    //   )
    // ));
  }
  /*
  The "zinbl" zero-inflation link function
  @param mu mean negative binomial abundance (real number vector)
  @param H height of curve at mode (real number > 0)
  @param g0 log-odds probability of excess zero when mu = 0
  @param g1 rate at which log-odds probability of excess zero decreases as mu increases
  @return vector of excess zero probabilities (real number vector in [0, 1])
  */
  vector zilink(vector mu, real H, real g0, real g1) {
    return inv_logit(fma(-g1 / H, mu, g0));
  }
  real std_halfnormal_rng() {
    real u = uniform_rng(0.5, 1);
    real y = inv_Phi(u);
    return y;
  }
