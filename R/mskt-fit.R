#' Fit the specified modskurt model by sampling its posterior distribution
#'
#' @param spec the modskurt model specification returned by `mskt_spec`
#' @param train_prop the proportion of data to train the posterior on
#' @param train_seed a seed for the training sample of proportion `train_prop`
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
                     train_prop = 1,
                     train_seed = NULL,
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
  sdat <- c(spec, spec$data(train_prop, train_seed)$train)
  sdat$data <- NULL
  if (is.null(init)) {
    init <- mskt_inits(sdat, chains)
  }
  # browser()
  # if (FALSE) {
  #   dat <- do.call(rbind,
  #                  lapply(init, \(i) as.data.frame(
  #                    purrr::discard(i, ~ length(.x) == 0))
  #                  ))
  #   out <-
  #     lapply(seq_along(init), \(i) {
  #       dat <- init[[i]]
  #       pars <- list()
  #       with(sdat, {
  #         pars$H <- as.double(dat$zH * y_max)
  #         pars$m <- as.double(dat$zm * x_pos_range + x_min)
  #         pars$s <- as.double(exp(dat$zs * hp_s[2] + hp_s[1]) * x_pos_range)
  #         pars$r <- as.double(dat$zr)
  #         pars$d <- as.double(exp(dat$zd * hp_d[2] + hp_d[1]))
  #         pars$p <- as.double(dat$zp * 1.99 + 0.01)
  #         pars
  #         data.frame(H = pars$H,
  #                    m = pars$m,
  #                    s = pars$s,
  #                    r = pars$r,
  #                    d = pars$d,
  #                    p = pars$p,
  #                    chain = i,
  #                    x = sdat$x,
  #                    mu = do.call(mskt, c(list(x = sdat$x), pars)))
  #       })
  #     })
  #   out[[4]]
  #   ggplot(do.call(rbind, out), aes(x, mu)) +
  #     geom_line() +
  #     facet_wrap(~ chain, scales = 'free_y')
  # }
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
  attr(fit, 'train') <- list(prop = train_prop, seed = train_seed)
  fit
}

# internal only for initialising from prior
mskt_inits <- function(sdat, chains, from_data = TRUE) {
  # https://discourse.mc-stan.org/t/undesirable-behavior-in-initial-values-of-hierarchical-model/18961/8?u=hdrab127
  # "Another example is multimodal likelihood functions – priors can suppress
  # extraneous modes but they can’t remove them entirely. Initializations near
  # an irrelevant mode will tend to get the Markov chains stuck near that mode
  # regardless of how little probability that mode captures. Custom
  # initializations can avoid these extraneous modes but again it’s up to the
  # user to verify that they are in fact extraneous."
  ypos <- sdat$y[sdat$y > 0]
  xpos <- sdat$x[sdat$y > 0]
  xwtd <- unlist(lapply(seq_along(xpos), \(j) rep(xpos[j], ypos[j])))

  # sample from prior for each variable
  inits <-
    lapply(seq_len(chains), function(cid) {
      with(sdat, {
        init <-
          list(zH = rbeta(1, hp_H[1], hp_H[2]),
               zm = rbeta(1, hp_m[1], hp_m[2]),
               zs = rnorm(1, hp_s[1], hp_s[2]),
               zr = as.array(rbeta(1, hp_r[1], hp_r[2])),
               zd = as.array(rnorm(1, hp_d[1], hp_d[2])),
               zp = as.array(min(0.1, rbeta(1, hp_p[1], hp_p[2]))),
               kap = rexp(1, hp_kap),
               zg0 = numeric(),
               zg1 = numeric())
        if (use_zi) {
          init$zg0 <- as.array(rnorm(1, hp_g0[1], hp_g0[2]))
          init$zg1 <- as.array(qnorm(runif(1, 0.5, 1), 0, hp_g1))
        }
        if (from_data) {
          # use some data driven initial values
          init$zH <- unname(quantile(ypos, runif(1, 0.2, 0.8))) / sdat$y_max
          init$zm <- (
            unname(quantile(xwtd, runif(1, 0.2, 0.8))) - sdat$x_pos_min
          ) / sdat$x_pos_range
          # when r approaches bounds it can cause inf vals
          init$zr <- as.array(stats::plogis(skewness(xwtd) *
                                              runif(1, 0.75, 1.25)))
          # TODO: better s from sd?
        }
        init
      })
    })
  inits
}
skewness <- function(x) {
  n <- length(x)
  x <- x - mean(x)
  r <- sqrt(n) * sum(x ^ 3) / (sum(x ^ 2) ^ (3 / 2)) *
    ((1 - 1 / n)) ^ (3 / 2)
  r
}
