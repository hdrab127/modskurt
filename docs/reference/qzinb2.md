# Zero-inflated Negative Binomial Quantile Function

Zero-inflated Negative Binomial Quantile Function

## Usage

``` r
qzinb2(p, mu, kappa, zi, lower.tail = TRUE, log.p = FALSE)
```

## Arguments

- p:

  vector of quantile probabilities

- mu:

  mean of the distribution

- kappa:

  dispersion (inverse square root size)

- zi:

  vector of (in \\\[0, 1\]\\) zero-inflation parameters

- lower.tail:

  logical; if TRUE (default), probabilities are \\P\[X \<= x\]\\,
  otherwise, \\P\[X \> x\]\\

- log.p:

  logical; if TRUE, probabilities p. This doesn't affect zi

## Value

A vector of quantiles, each of which coincide with the respective
probability in p
