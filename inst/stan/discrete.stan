functions {
#include mu-rdp/mskt.stan
}
data {
#include mu-rdp/data.stan
  // integer abundance y 
  array[N] int<lower=0> y;
  real<lower=0> hp_kap;

#include zi-logistic-reg/data-zi.stan
}
transformed data {
#include mu-rdp/tdata.stan
}
parameters {
#include mu-rdp/pars.stan
#include zi-logistic-reg/pars-zi.stan
  // hopefully that upper bound restricts phi >= eps and avoids
  // underflow complaints MH tries kap -> inf
  real<lower=machine_precision(),upper=1/sqrt(machine_precision())> kap;
}
transformed parameters {
#include mu-rdp/tpars.stan
  // more numerically stable than pow(kap, -2)?
  real phi = 1 / square(kap);
#include zi-logistic-reg/tpars-zi.stan
}
model {
#include mu-rdp/model.stan
  target += exponential_lupdf(kap | hp_kap);
  if (use_zi) {
#include zi-logistic-reg/model-zi.stan
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
#include mu-rdp/genq.stan
#include zi-logistic-reg/genq-zi.stan
  array[N] int y_rep;
  array[N] int pr_y_rep;

  real pr_kap;
  real pr_phi;

  if (sample_prior < 2) {
    mu_rep = mskt((xrep - m) / s, H, r, d, p);
    // calc log lik and yreps
    if (use_zi) {
#include zi-logistic-reg/genq-zi-rep.stan
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
#include mu-rdp/genq-pr-mu.stan
#include zi-logistic-reg/genq-pr-zi.stan
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
