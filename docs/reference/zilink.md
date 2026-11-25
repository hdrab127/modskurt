# ZINBL zero-inflation link function

ZINBL zero-inflation link function

## Usage

``` r
zilink(mu, H, g0, g1)
```

## Arguments

- mu:

  mean negative binomial abundance (real number vector)

- H:

  height of curve at mode (real number \> 0)

- g0:

  log-odds probability of excess zero when mu = 0

- g1:

  rate at which log-odds probability of excess zero decreases as mu
  increases

## Value

vector of excess zero probabilities (real number vector in \[0, 1\])
