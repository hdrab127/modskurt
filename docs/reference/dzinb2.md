# Zero-inflated Negative Binomial Probability Mass

Zero-inflated Negative Binomial Probability Mass

## Usage

``` r
dzinb2(x, mu, kappa, zi, log = FALSE)
```

## Arguments

- x:

  vector of (non-negative integer) values

- mu:

  mean of the distribution

- kappa:

  dispersion (inverse square root size)

- zi:

  vector of (real lying in \\0, 1\\) zero-inflation parameters

- log:

  logical; if TRUE, probabilities, p, are given as log(p)

## Value

The (log) probability mass at x, given mu, kappa and zi
