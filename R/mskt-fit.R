#' Fit the specified modskurt model by sampling its posterior distribution
#'
#' @param spec the modskurt model specification returned by `mskt_spec`
#' @param use_subset whether to use specified subset or all data
#' @param seed a seed for the Stan random number generator
#' @param init initialisation method, see ?cmdstanr::`model-method-sample`
#' @param chains the number of Markov chains
#' @param parallel_chains how many cores to distribute chains over
#' @param iter_warmup number of iterations used to learn posterior curvature
#' @param iter_sampling number of posterior samples drawn
#' @param max_treedepth,adapt_delta tuning parameters, see details
#' @param show_messages additional verbosity
#' @param show_exceptions additional verbosity
#' @param ... additional arguments to ?cmdstanr::`model-method-sample`
#'
#' @return A `cmdstanr::CmdStanMCMC` object
#' @export
#'
mskt_fit <- function(spec,
                     use_subset = FALSE,
                     seed = NULL,
                     init = NULL,
                     chains = 4,
                     parallel_chains = getOption('mc.cores', 1),
                     iter_warmup = 1000,
                     iter_sampling = 1000,
                     max_treedepth = NULL,
                     adapt_delta = NULL,
                     show_messages = TRUE,
                     show_exceptions = FALSE,
                     ...) {
  # retrieve cmdstan_model
  pkgdir <- path.package('modskurt')
  standir <- file.path(pkgdir, 'stan')

  # TODO: remove this from testing
  if (!dir.exists(standir)) {
    standir <- file.path(gsub('OneDrive - Massey University',
                              'msc',
                              pkgdir,
                              fixed = TRUE), 'inst/stan')
  }

  stanexe <- file.path(standir, 'discrete.exe')
  if (!file.exists(stanexe)) {
    rlang::abort('modskurt models not compiled, see `?compile_stanmodels()')
  }

  suppressMessages({
    csm <-
      cmdstanr::cmdstan_model(stan_file = file.path(standir, 'discrete.stan'),
                              exe_file = stanexe,
                              include_paths = file.path(standir, '/blocks'))
  })

  if (!missing(seed) & !is.null(seed)) {
    set.seed(seed)
  }

  if (use_subset) {
    sdat <- c(spec, spec$subset()$train)
  } else {
    sdat <- c(spec, spec$full()$train)
  }
  sdat$subset <- NULL
  sdat$full <- NULL

  if (is.null(init)) {
    init <- mskt_inits(sdat, chains)
  }
  fit <-
    csm$sample(data = sdat,
               seed = seed,
               init = init,
               chains = chains,
               parallel_chains = parallel_chains,
               iter_warmup = iter_warmup,
               iter_sampling = iter_sampling,
               max_treedepth = max_treedepth,
               adapt_delta = adapt_delta,
               show_messages = show_messages,
               show_exceptions = show_exceptions,
               ...)
  attr(fit, 'spec') <- spec
  attr(fit, 'is_subset') <- use_subset
  fit
}

#' internal only for initialising from prior
#' @keywords internal
mskt_inits <- function(sdat, chains, from_data = TRUE) {
  # https://discourse.mc-stan.org/t/undesirable-behavior-in-initial-values-of-hierarchical-model/18961/8?u=hdrab127
  # "Another example is multimodal likelihood functions – priors can suppress
  # extraneous modes but they can’t remove them entirely. Initializations near
  # an irrelevant mode will tend to get the Markov chains stuck near that mode
  # regardless of how little probability that mode captures. Custom
  # initializations can avoid these extraneous modes but again it’s up to the
  # user to verify that they are in fact extraneous."
  ypos <- sdat$y[sdat$y > 0] / sdat$eff[sdat$y > 0]
  xpos <- sdat$x[sdat$y > 0]
  xwtd <- unlist(lapply(seq_along(xpos), \(j) rep(xpos[j], ypos[j])))

  # sample from prior for each variable
  inits <-
    lapply(seq_len(chains), function(cid) {
      with(sdat, {
        init <-
          # try keep H at least one
          list(zH = max(c(1 / y_max, stats::rbeta(1, hp_H[1], hp_H[2]))),
               zm = stats::rbeta(1, hp_m[1], hp_m[2]),
               zs = 0.15 + rhnorm(1),
               zr = as.array(stats::rbeta(1, hp_r[1], hp_r[2])),
               zd = as.array(stats::rnorm(1)),
               # low p values create undefined mu headaches
               zp = as.array(max(c(0.25, stats::rbeta(1, hp_p[1], hp_p[2])))),
               kap = stats::rexp(1, hp_kap),
               zg0 = numeric(),
               zg1 = numeric())
        if (use_zi) {
          # if using zero-inflation (and it is actually present)
          # should always have very high excess-zero intercept to avoid fighting
          # with kappa
          init$zg0 <- as.array(stats::qlogis(0.999))
          init$zg1 <- as.array(rhnorm(1))
        }
        if (from_data) {
          # use some data driven initial values
          init$zH <- max(c(1, unname(quantile(ypos,
                                              stats::runif(1, 0.3, 0.7))))) /
            sdat$y_max
          init$zm <- (
            unname(quantile(xwtd, stats::runif(1, 0.3, 0.7))) - sdat$x_pos_min
          ) / sdat$x_pos_range
          if (is.na(init$zm) || init$zm == 1) {
            init$zm <- rbeta(1, hp_m[1], hp_m[2])
          }
          # when r approaches bounds it can cause inf vals
          init$zr <- as.array(stats::plogis(skewness(xwtd) *
                                              stats::runif(1, 0.75, 1.25)))
          if (is.na(init$zr)) {
            init$zr <- as.array(stats::rbeta(1, hp_r[1], hp_r[2]))
          }
          # TODO: better s from sd?
        }
        init
      })
    })
  inits
}
#' @keywords internal
skewness <- function(x) {
  n <- length(x)
  x <- x - mean(x)
  r <- sqrt(n) * sum(x ^ 3) / (sum(x ^ 2) ^ (3 / 2)) *
    ((1 - 1 / n)) ^ (3 / 2)
  r
}
