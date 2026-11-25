# The Negative Binomial Quantile Function (alternative parameterisation)

The Negative Binomial Quantile Function (alternative parameterisation)

## Usage

``` r
qnbinom2(p, mu, kappa, lower.tail = TRUE, log.p = FALSE)
```

## Arguments

- p:

  vector of probabilities

- mu:

  mean

- kappa:

  dispersion (inverse square root size)

- lower.tail:

  logical; if TRUE (default), probablities are \\P\[X \<= x\]\\,
  otherwise, \\\[X \> x\]\\

- log.p:

  logical; if TRUE, probabilities are given as log(p)

## Value

A vector of quantiles, each of which correspond to a probability in p
