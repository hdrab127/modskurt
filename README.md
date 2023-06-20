# modskurt

## this package has five core functions

[] specify Bayesian modskurt models for discrete abundance data (nb or zinbl)
[] set and check prior specification
[] compute and validate fit
[] check posterior specification
[] use posterior (visual relationship, thresholds, intervals)

## getting started

- Install `cmdstanr` (instructions here)
- Install `modskurt`:

```
# install development version of package from github
remotes::install_github('modskurt')

# one-off compilation of stan files for Bayesian modskurt models
# optionally edit path...
modskurt::compile_stanmodels()

# test the basics work
library(modskurt)
fit <-
  fit_modskurt(data = data.frame(x = rnorm(50),
                                 y = rnbinom(50, mu = 5, size = 2)),
               shape = 'rdp',
               dist = 'nb')
pp_check(fit)
```
