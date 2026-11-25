# Stan code for the discrete model

The following code in the Stan probabilistic programming language
defines a model where the discrete species abundance, \\y\\, follows a
zero-inflated negative binomial distribution with probability of
excess-zero and mean abundance both linked to the data through
`modskurt` mean functions of the environmental gradient, \\x\\.

``` stan
functions {
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
}
data {
  int<lower=1> N;
  int<lower=1> Nz;
  vector[N] x;

  // log effort (can be 0 if no offset required)
  vector[N] leff;

  // make predictions over grid
  int Nrep;
  vector[Nrep] xrep;

  // re-scaling
  real<lower=0> y_max;
  real<lower=0> x_pos_range;
  real x_pos_min;

  // hyperpriors
  vector[2] hp_H;
  vector[2] hp_m;
  vector[2] hp_s;
  vector[2] hp_r;
  vector[2] hp_d;
  vector[2] hp_p;

  // toggle optional vars H,m,s[,r,d,p,a]
  int<lower=0,upper=1> use_r;
  int<lower=0,upper=1> use_d;
  int<lower=0,upper=1> use_p;

  // specify fixed values for parameters if not being used^
  real<lower=0,upper=1>   fix_r;
  real<lower=0>           fix_d;
  real<lower=0.1>         fix_p;

  // 0: no prior, 1: with prior, 2: only prior
  int<lower=0,upper=2> sample_prior;
  // base distribution
  array[N] int<lower=0> y;
  real<lower=0> hp_kap;

  int<lower=0,upper=1> use_zi;
  vector[2] hp_g0;
  vector[1] hp_g1;
}
transformed data {
  // data type specific perhaps
  int Nzp = Nz + 1;

  // standardise
  // real z_min = x_min / x_sd;
  // real z_range = x_range / x_sd;
  // vector[N] z = x / x_sd;
  // vector[Nrep] zrep = xrep / x_sd;
}
parameters {
  // standard-normal rvs for non-centred par
  real<lower=0,upper=1> zH;
  real<lower=0,upper=1> zm;
  real zs;
  array[use_r] real<lower=0,upper=1> zr;
  array[use_d] real                  zd;
  array[use_p] real<lower=0,upper=1> zp;
  array[use_zi] real<lower=0,upper=1> g0;
  array[use_zi] real<lower=0> zg1;
  real<lower=0> kap;
}
transformed parameters {
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
  vector[N] lmu = leff + lmskt((x - m) / s, H, r, d, p);
  array[use_zi] real g1;
  vector[N] zi;
  if (use_zi) {
    // logistic function of the relative mean
    g1[1] = hp_g1[1] * zg1[1];
    zi = inv_logit(fma(-g1[1] / H, exp(lmu), logit(g0[1])));
  }
  real phi = pow(kap, -2);
} 
model {
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
  target += exponential_lupdf(kap | hp_kap);

  if (use_zi) {
    target += beta_lupdf(g0[1] | hp_g0[1], hp_g0[2]);
    target += std_normal_lupdf(g1[1]);
    if (Nz > 0) {
      target += log_sum_exp(bernoulli_lupmf(1 | zi[1:Nz]),
                            bernoulli_lupmf(0 | zi[1:Nz]) +
                            neg_binomial_2_log_lupmf(0 | lmu[1:Nz], phi));
    }
    target += bernoulli_lupmf(0 | zi[Nzp:N]) +
              neg_binomial_2_log_lupmf(y[Nzp:N] | lmu[Nzp:N], phi);
  } else {
    target += neg_binomial_2_log_lupmf(y | lmu, phi);
  }
}
generated quantities {
  vector[N] mu = exp(lmu);
  vector[N] log_lik;
  vector[Nrep] mu_rep;
  vector[Nrep] pr_mu_rep;
  real pr_H;
  real pr_m;
  real pr_s;
  real pr_r;
  real pr_d;
  real pr_p;
  // TODO: this should be log_prior_mu
  // and have one for log_prior_phi
  // and log_prior_zi for other dist checks
  real log_prior = 0;
  real pr_g0;
  real pr_g1;
  vector[Nrep] zi_rep;
  vector[Nrep] pr_zi_rep;
  array[N] int y_rep;
  array[N] int pr_y_rep;

  real pr_kap;
  real pr_phi;

  if (sample_prior < 2) {
    mu_rep = exp(lmskt((xrep - m) / s, H, r, d, p));
    // calc log lik and yreps
    if (use_zi) {
      zi_rep = inv_logit(fma(-g1[1] / H, mu_rep, logit(g0[1])));
      for (n in 1:N) {
        if (n < Nzp) {
          log_lik[n] = log_sum_exp(bernoulli_lpmf(1 | zi[n]),
                                   bernoulli_lpmf(0 | zi[n]) +
                                   neg_binomial_2_log_lpmf(0 | lmu[n], phi));
        } else {
          log_lik[n] = bernoulli_lpmf(0 | zi[n]) +
                       neg_binomial_2_log_lpmf(y[n] | lmu[n], phi);
        }
        if (zi[n] > machine_precision() && bernoulli_rng(zi[n]) > 0.5) {
          y_rep[n] = 0;
        } else {
          y_rep[n] = neg_binomial_2_log_rng(lmu[n], phi);
        }
      }
    } else {
      for (n in 1:N) {
        log_lik[n] = neg_binomial_2_log_lpmf(y[n] | lmu[n], phi);
        if (mu[n] == 0) {
          y_rep[n] = 0;
        } else {
          y_rep[n] = neg_binomial_2_rng(mu[n], phi);
        }
      }
    }
  }

  if (sample_prior > 0) {
    // calcs log_prior for mean pars now 
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
    // calc zi & yreps
    vector[N] pr_zi;
    if (use_zi) {
      pr_g0 = beta_rng(hp_g0[1], hp_g0[2]);
      pr_g1 = hp_g1[1] * std_halfnormal_rng();
      pr_zi = inv_logit(fma(-pr_g1 / pr_H, pr_mu, logit(pr_g0)));
      pr_zi_rep = inv_logit(fma(-pr_g1 / pr_H, pr_mu_rep, logit(pr_g0)));
    }
    pr_kap = exponential_rng(hp_kap);
    pr_phi = pow(pr_kap, -2);
    for (n in 1:N) {
      if (is_nan(pr_mu[n])) {
        pr_y_rep[n] = -999;
      } else if (pr_mu[n] <= machine_precision() ||
          (!is_nan(pr_zi[n]) &&
           pr_zi[n] > machine_precision() &&
           bernoulli_rng(pr_zi[n]) > 0.5)) {
        pr_y_rep[n] = 0;
      } else {
        pr_y_rep[n] = neg_binomial_2_rng(pr_mu[n], pr_phi);
      }
    }
  }
}
```
