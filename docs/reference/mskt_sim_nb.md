# Simulate negative binomial data with a modskurt mean curve

Simulate negative binomial data with a modskurt mean curve

## Usage

``` r
mskt_sim_nb(
  x,
  n_per_x = 1,
  H,
  m,
  s,
  r = 0.5,
  d = 0,
  p = 1,
  kap,
  g0 = -Inf,
  g1 = 0
)
```

## Arguments

- x:

  vector of environmental gradient values

- n_per_x:

  number of y replicates to simulate for each x value

- H:

  height

- m:

  middle

- s:

  spread

- r:

  assymetry

- d:

  flatness

- p:

  persistence

- kap:

  negative binomial dispersion

- g0:

  constant rate of zero-inflation

- g1:

  decay of zero-inflation as mean increases

## Value

data.frame with simulated mu, zi, and y for each x
