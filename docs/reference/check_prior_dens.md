# Check marginal prior densities are encoded correctly

Check marginal prior densities are encoded correctly

## Usage

``` r
check_prior_dens(
  spec,
  pars = attr(spec, "pars"),
  use_subset = TRUE,
  prob_inner = 0.5,
  prob_outer = 0.99,
  prob_annotate = 0.95,
  plot = TRUE
)
```

## Arguments

- spec:

  the modskurt model specification returned by \`mskt_spec\`

- pars:

  parameters to plot densities for, default is all

- use_subset:

  whether to use specified subset or all data

- prob_inner:

  inner probability (high density interval) to highlight

- prob_outer:

  outer probability range to plot the density between

- prob_annotate:

  optional, add markers at certain probability interval

- plot:

  boolean, plot together, or return list of individual ggplots

## Value

a patchwork of ggplots
