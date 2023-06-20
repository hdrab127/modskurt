#' Fit the specified ModSkurt model by sampling its posterior distribution
#'
#' @param spec the ModSkurt model specification returned by `modskurt_spec`
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
modskurt_fit <- function(spec,
                         seed = 1234,
                         init = NULL,
                         chains = 6L,
                         parallel_chains = 6L,
                         iter_warmup = NULL,
                         iter_sampling = NULL,
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
    rlang::abort('ModSkurt models not compiled, see `?compile_stanmodels()')
  }
  csm <-
    cmdstanr::cmdstan_model(stan_file = file.path(standir, 'discrete.stan'),
                            exe_file = stanexe,
                            include_paths = file.path(standir, '/blocks'))
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (is.null(init)) {
    init <- modskurt_inits(spec, chains)
  }
  fit <-
    csm$sample(data = spec,
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
  fit
}

# internal only for initialising from prior
modskurt_inits <- function(spec, chains) {
  # https://discourse.mc-stan.org/t/undesirable-behavior-in-initial-values-of-hierarchical-model/18961/8?u=hdrab127
  # "Another example is multimodal likelihood functions – priors can suppress
  # extraneous modes but they can’t remove them entirely. Initializations near
  # an irrelevant mode will tend to get the Markov chains stuck near that mode
  # regardless of how little probability that mode captures. Custom
  # initializations can avoid these extraneous modes but again it’s up to the
  # user to verify that they are in fact extraneous."
  ypos <- spec$y[spec$y > 0]
  xpos <- spec$x[spec$y > 0]
  xwtd <- unlist(lapply(seq_along(xpos), \(j) rep(xpos[j], ypos[j])))

  # sample from prior for each variable
  inits <-
    lapply(seq_len(chains), function(cid) {
      with(spec, {
        init <-
          list(zH = rbeta(1, hp_H[1], hp_H[2]),
               zm = rbeta(1, hp_m[1], hp_m[2]),
               zs = rnorm(1, hp_s[1], hp_s[2]),
               zr = as.array(rbeta(1, hp_r[1], hp_r[2])),
               zd = as.array(rnorm(1, hp_d[1], hp_d[2])),
               zp = as.array(rbeta(1, hp_p[1], hp_p[2])),
               kap = rexp(1, hp_kap),
               g0 = numeric(),
               zg1 = numeric())
        if (use_zi) {
          init$g0 <- as.array(rnorm(1, hp_g0[1], hp_g0[2]))
          init$g1 <- as.array(qnorm(runif(1, 0.5, 1), 0, hp_g1))
        }
        init
      })
    })
  inits
}
