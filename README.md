# modskurt <a href="https://hdrab127.github.io/modskurt/"></a>

Bayesian modelling of species abundance distributions along environmental gradients

## Installing `modskurt`

The `modskurt` package utilises the Stan Bayesian modelling software implemented in R through the `cmdstanr` package.

First, install the `cmdstanr` interface (if not already) by running

```r
install.packages('cmdstanr', repos = c('https://mc-stan.org/r-packages/', getOption('repos')))
```

Then to install `cmdstan` itself, try

```r
library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE)
install_cmdstan(cores = 2)
```

For any hiccups in the above, see https://mc-stan.org/cmdstanr/articles/cmdstanr.html.

Now to install `modskurt` (currently beta version on github only), try running

```r
# install.packages('pak')
pak::pkg_install('hdrab/modskurt')
```

And finally, to compile the `modskurt` Stan models specifically for your computer, run

```r
library(modskurt)
# only needs running after new package install or package update
compile_stanmodels()
```

## Fitting a model

To check the installation worked, try fitting a simple negative binomial regression with assymetric modskurt mean

```r
# model specification with fake data
spec <-
  mskt_spec(data = mskt_sim_nb(x = 0:50, 
                               H = 100,
                               m = 35,
                               s = 40,
                               r = 0.9,
                               kap = 0.5),
            y = c('Abundance (count)' = 'y'),
            x = c('Env gradient' = 'x'),
            dist = 'nb')
# fit a minimal model (with less than optimal iterations)
fit <- mskt_fit(spec,
                iter_warmup = 200,
                iter_sampling = 100)

# check the posterior summary and computation diagnostics
check_computation(fit)

# plot the abundance distribution
abundance_dist(fit)
```

## How to use the `modskurt` package

Statistical modelling often requires more than just computing algorithms, see [Getting started](./articles/getting-started.html) for a recommended workflow that takes species-environment data through a Bayesian workflow to incorporate ecological knowledge and instil trust in results.

## Roadmap and contribution

The current development version of this `modskurt` package is very much in its infancy. Planned additions in a vague order of priority are:

- [ ] &nbsp; More vignettes and examples
- [ ] &nbsp; Model comparison (e.g. between nb and zinbl)
- [ ] &nbsp; Continuous abundance (like biomass)
- [ ] &nbsp; Binary package for Cran or drat
- [ ] &nbsp; Random effects
- [ ] &nbsp; Temporal effects
- [ ] &nbsp; Multiple gradients (2d modskurts)
- [ ] &nbsp; Multiple species
- [ ] &nbsp; ...
