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
