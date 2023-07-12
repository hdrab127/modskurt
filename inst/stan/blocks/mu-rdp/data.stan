  int<lower=1> N;
  int<lower=1> Nz;
  vector[N] x;

  // effort (can be 1 if no offset required)
  vector[N] eff;

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
