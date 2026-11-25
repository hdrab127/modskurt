# Specify modskurt model for the distribution of \`y\` over \`x\`

Specify modskurt model for the distribution of \`y\` over \`x\`

## Usage

``` r
mskt_spec(
  data,
  y,
  x,
  effort = NULL,
  dist = "nb",
  shape = "rdp",
  hyperparams = list(),
  pred_grid = NULL,
  subset_prop = 1,
  subset_seed = NULL
)
```

## Arguments

- data:

  data.frame or tibble

- y:

  species count, optionally named - see details

- x:

  environmental gradient, optionally named - see details

- effort:

  sampling effort used to observe y

- dist:

  distribution of y over x

- shape:

  shape of distribution mean - see ?mskt_shape()

- hyperparams:

  list of hyper-parameters to adjust in-built priors

- pred_grid:

  vector of x values to make y predictions over

- subset_prop:

  proportion of data to use for early model checks

- subset_seed:

  a seed for the subset sample proportion

## Value

list of standata for prior verification and posterior estimation

## Details

\## Forms of the discrete (count) data distribution - \`"nb"\`: Negative
Binomial - \`"zinbl"\`: Zero-inflated Negative Binomial (Linked)

An interactive graph of these count distributions and possible
parameters is available at
\[https://salt-ecology.shinyapps.io/nb-zinbl-prior-specs/\]().

\## Shapes for the mean function - default: height (H), mode (m), scale
(s) - optional extras + \`"r"\` assymmetry or skew + \`"d"\` flatness or
peakedness + \`"p"\` tail exaggeration or pinching

These can be combined as desired, e.g. \`"rp"\`, \`"rdp"\`, ...

An interactive graph of the modskurt mean function and possible
parameters is available at
\[https://salt-ecology.shinyapps.io/modskurt-prior-specs/\]().
