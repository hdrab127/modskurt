#' Specify modskurt model for the distribution of `y` over `x`
#'
#' @param data data.frame or tibble
#' @param y species count, optionally named - see details
#' @param x environmental gradient, optionally named - see details
#' @param effort sampling effort used to observe y
#' @param dist distribution of y over x
#' @param shape shape of distribution mean - see ?mskt_shape()
#' @param hyperparams list of hyper-parameters to adjust in-built priors
#' @param pred_grid vector of x values to make y predictions over
#' @param subset_prop proportion of data to use for early model checks
#' @param subset_seed a seed for the subset sample proportion
#'
#' @details
#'
#' ## Forms of the discrete (count) data distribution
#' - `"nb"`: Negative Binomial
#' - `"zinbl"`: Zero-inflated Negative Binomial (Linked)
#'
#' An interactive graph of these count distributions and possible parameters is
#' available at [https://salt-ecology.shinyapps.io/nb-zinbl-prior-specs/]().
#'
#' ## Shapes for the mean function
#' - default: height (H), mode (m), scale (s)
#' - optional extras
#'   + `"r"` assymmetry or skew
#'   + `"d"` flatness or peakedness
#'   + `"p"` tail exaggeration or pinching
#'
#' These can be combined as desired, e.g. `"rp"`, `"rdp"`, ...
#'
#' An interactive graph of the modskurt mean function and possible parameters is
#' available at [https://salt-ecology.shinyapps.io/modskurt-prior-specs/]().
#'
#' @return list of standata for prior verification and posterior estimation
#' @export
#'
mskt_spec <- function(data,
                      y,
                      x,
                      effort = NULL,
                      dist = 'nb',
                      shape = 'rdp',
                      hyperparams = list(),
                      pred_grid = NULL,
                      subset_prop = 1,
                      subset_seed = NULL) {
  # only discrete data currently
  d_all <- data.frame(y = as.integer(data[[y]]),
                      x = data[[x]],
                      eff = 1)
  if (!is.null(effort)) {
    # this just affects rate, not offset
    # https://discourse.mc-stan.org/t/
    # scaling-of-the-overdispersion-in-negative-binomial-models/15581
    d_all$eff <- data[[effort]]
  }

  # check missing
  okay <- stats::complete.cases(d_all)
  if (any(!okay)) {
    rlang::warn(paste('Removed', sum(!okay), 'rows containing missing values.'))
    d_all <- d_all[okay, ]
  }

  # set prior parameters
  hp <- hyperparams
  spec <-
    #    # dont sample prior in fit mode
    list(sample_prior = 0L,
         # mean function
         hp_H = hp$H %||% c( 2.0, 3.0), # ~ beta(a, b)
         hp_m = hp$m %||% c( 1.0, 1.0), # ~ beta(a, b)
         hp_s = hp$s %||% c(-2.0, 0.6), # ~ log_normal(mu, sd)
         hp_r = hp$r %||% c( 1.2, 1.2), # ~ beta(a, b)
         hp_d = hp$d %||% c( 0.5, 1.0), # ~ log_normal(mu, sd)
         hp_p = hp$p %||% c( 1.2, 1.2), # ~ beta(a, b)
         # use or set fixed "off" positions
         use_r = as.integer(grepl('r', shape)),
         fix_r = 0.5,
         use_d = as.integer(grepl('d', shape)),
         fix_d = 0.0,
         use_p = as.integer(grepl('p', shape)),
         fix_p = 1.0,
         # distribution
         use_zi = as.integer(dist == 'zinbl'),
         hp_kap = hp$kap %||% c(0.5), # ~ exponential(rate)
         hp_g0  =  hp$g0 %||% c(3.0, 1.5), # ~ normal(mu, sd)
         hp_g1  =  hp$g1 %||% c(5.0)) # ~ half_normal(mu, sd)

  # default prediction grid over all x
  if (is.null(pred_grid)) {
    spec$Nrep <- 50L
    spec$xrep <- seq(min(d_all$x), max(d_all$x), length.out = spec$Nrep)
  } else {
    spec$Nrep <- length(pred_grid)
    spec$xrep <- pred_grid
  }

  # standardise y, x over prop of data set
  # enclose seed for consistent subsets if not specified
  curr_seed <- .Random.seed
  get_subset <- function(prop = 1, seed = NULL) {
    # optionally sample a proportion of training data
    if (!missing(seed) & !is.null(seed)) {
      set.seed(seed)
    } else {
      set.seed(curr_seed)
    }
    idx <- sample(nrow(d_all), round(prop * nrow(d_all)))
    d <- d_all
    d$set <- 'test'
    d$set[idx] <- 'train'
    train <- d[idx, ]

    # arrange for faster stan zinbl
    train <- train[order(train$y), ]
    train <- as.list(train)
    train <- list(y = train$y,
                  x = train$x,
                  eff = train$eff,
                  N = length(train$y),
                  Nz = sum(train$y == 0),
                  x_range = diff(range(train$x)),
                  x_pos_range = diff(range(train$x[train$y > 0])),
                  x_min = min(train$x),
                  x_pos_min = min(train$x[train$y > 0]),
                  y_max = max(train$y))


    # check zeroes
    prop_zero <- (train$Nz / train$N)
    if (prop_zero > 0.90 & dist != 'zinbl') {
      rlang::warn(paste0(round(prop_zero * 100, 1), '% of y are zero. ',
                         'Maybe try a zero-inflated model instead?'))
    }
    list(train = train, all = d, prop = prop, seed = seed)
  }
  spec$subset <- function() { get_subset(subset_prop, subset_seed) }
  spec$full <- get_subset

  pars <- c('H', 'm', 's', strsplit(shape, '')[[1]], 'kap')
  if (spec$use_zi) {
    pars <- c(pars, 'g0', 'g1')
  }
  attr(spec, 'pars') <- pars
  attr(spec, 'nms') <- list(y = names(y) %||% y,
                            x = names(x) %||% x,
                            eff = names(effort) %||% effort)
  attr(spec, 'dist') <- dist
  attr(spec, 'shape') <- shape
  class(spec) <- c('mskt-spec', class(spec))
  spec
}
