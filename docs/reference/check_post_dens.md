# visual display of chain mixing and prior-data conflicts

visual display of chain mixing and prior-data conflicts

## Usage

``` r
check_post_dens(
  fit,
  pars = attr(spec, "pars"),
  prob_inner = 0.5,
  prob_outer = 0.99,
  by_chain = FALSE,
  chain_legend_ncol = 1,
  plot = TRUE
)
```

## Arguments

- fit:

  fitted model object from \`mskt_fit\`

- pars:

  parameters to plot densities for, default is all

- prob_inner:

  inner probability (high density interval) to highlight

- prob_outer:

  outer probability range to plot the density between

- by_chain:

  plot chains separately to identify multiple modes and mixing issues

- chain_legend_ncol:

  for the thesis figure

- plot:

  boolean, plot together, or return list of individual ggplots
