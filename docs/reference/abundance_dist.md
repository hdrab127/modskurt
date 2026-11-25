# Distribution of abundance along the environmental gradient

Distribution of abundance along the environmental gradient

## Usage

``` r
abundance_dist(
  fit,
  seed = NULL,
  ndraws = 50,
  summaries = c("mean", "q90"),
  include_zero_inflation = TRUE,
  point_alpha = 0.2,
  line_alpha = 0.1
)
```

## Arguments

- fit:

  modskurt_fit object

- seed:

  optional seed for ndraws sample

- ndraws:

  optional subset of posterior draws

- summaries:

  character vector of distribution summarise, see details

- include_zero_inflation:

  whether to include zero-inflation in summaries

- point_alpha:

  opacity of data points, default 0.1 is good for ~ 1000 obs

- line_alpha:

  opacity of summary lines, default 0.1 is good for ~ 50 draws

## Value

a ggplot object
