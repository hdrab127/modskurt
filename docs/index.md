# modskurt [](https://hdrab127.github.io/modskurt/)

Bayesian modelling of species abundance distributions along
environmental gradients

## Installing `modskurt`

The `modskurt` package utilises Stan Bayesian modelling software
implemented in R through the `cmdstanr` package.

In the early stage of package development, `modskurt` is focused on the
latest advancements, rather than backwards compatibility. This means
you’ll need the latest version of R (\>=4.3.1), Rtools (\>=4.3) and
development version of `cmdstanr` (\>=0.5.3) to get going.

Check your `R.version` and Rtools version
(`devtools::has_devel(debug = TRUE)`) numbers and update if needed.

Then to install the development version of the `cmdstanr` interface, try
running:

``` r

# install.packages('pak') # <- if needed
pak::pkg_install('stan-dev/cmdstanr')
```

The `cmdstan` program that `cmdstanr` provides an interface to is
installed using:

``` r

library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE)
install_cmdstan(cores = 2)
```

For any hiccups and errors, follow the guide on the `cmdstanr` website:
<https://mc-stan.org/cmdstanr/articles/cmdstanr.html>, or log an
[issue](https://github.com/hdrab127/modskurt/issues) on this `modskurt`
github.

Now to install the alpha version of `modskurt`:

``` r

pak::pkg_install('hdrab127/modskurt')
```

And finally, to compile the `modskurt` Stan models on your computer:

``` r

library(modskurt)
# only compiles when needed, e.g. after package install or update
compile_stanmodels()
```

## Checking the installation

To check the installation worked, try fitting this simple negative
binomial regression with an assymetric `modskurt` mean shape. Don’t
worry about the parameters or functions just yet - everything will be
explained in the [Getting
started](https://hdrab127.github.io/modskurt/articles/getting-started.html)
guide.

``` r

# model specification with fake data
spec <-
  mskt_spec(data = mskt_sim_nb(x = 0:50, 
                               H = 100,
                               m = 35,
                               s = 40,
                               r = 0.9,
                               kap = 1),
            y = c('Abundance (count)' = 'y'),
            x = c('Env gradient' = 'x'),
            dist = 'nb')

# fit a minimal model (with less than optimal iterations)
fit <- 
  mskt_fit(spec,
           iter_warmup = 200,
           iter_sampling = 100)

# check the posterior summary and computation diagnostics
check_computation(fit)

# plot the abundance distribution
abundance_dist(fit)
```

## How to use the `modskurt` package

See [Getting
started](https://hdrab127.github.io/modskurt/articles/getting-started.html)
for a worked example and recommended workflow for analysing these
nonlinear `modskurt` models of the distribution of a species’ abundance
along environmental gradients.

More articles will follow.

## Roadmap and contribution

The current development version of this `modskurt` package is very much
in its infancy. Potential additions in a vague order of priority are:

  More vignettes and examples

  Model comparison (e.g. between nb and zinbl)

  Continuous abundance (like biomass)

  Binary package for Cran or drat

  Random effects

  Temporal effects

  Multiple gradients (2d modskurts)

  Multiple species

  …
