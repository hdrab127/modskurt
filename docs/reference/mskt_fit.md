# Fit the specified modskurt model by sampling its posterior distribution

Fit the specified modskurt model by sampling its posterior distribution

## Usage

``` r
mskt_fit(
  spec,
  use_subset = FALSE,
  seed = NULL,
  init = NULL,
  chains = 4,
  parallel_chains = getOption("mc.cores", 1),
  iter_warmup = 1000,
  iter_sampling = 1000,
  max_treedepth = NULL,
  adapt_delta = NULL,
  show_messages = TRUE,
  show_exceptions = FALSE,
  ...
)
```

## Arguments

- spec:

  the modskurt model specification returned by \`mskt_spec\`

- use_subset:

  whether to use specified subset or all data

- seed:

  a seed for the Stan random number generator

- init:

  initialisation method, see ?cmdstanr::\`model-method-sample\`

- chains:

  the number of Markov chains

- parallel_chains:

  how many cores to distribute chains over

- iter_warmup:

  number of iterations used to learn posterior curvature

- iter_sampling:

  number of posterior samples drawn

- max_treedepth, adapt_delta:

  tuning parameters, see details

- show_messages:

  additional verbosity

- show_exceptions:

  additional verbosity

- ...:

  additional arguments to ?cmdstanr::\`model-method-sample\`

## Value

A \`cmdstanr::CmdStanMCMC\` object
