#' Check marginal prior densities are encoded correctly
#'
#' @param spec the modskurt model specification returned by `mskt_spec`
#' @param pars parameters to plot densities for, default is all
#' @param use_subset whether to use specified subset or all data
#' @param prob_inner inner probability (high density interval) to highlight
#' @param prob_outer outer probability range to plot the density between
#' @param prob_annotate optional, add markers at certain probability interval
#' @param plot boolean, plot together, or return list of individual ggplots
#'
#' @importFrom rlang .data
#'
#' @return a patchwork of ggplots
#' @export
check_prior_dens <- function(spec,
                             pars = attr(spec, 'pars'),
                             use_subset = TRUE,
                             prob_inner = 0.5,
                             prob_outer = 0.99,
                             prob_annotate = 0.95,
                             plot = TRUE) {
  # outer prob (op), inner prob (ip)
  oplow <- (1 - prob_outer) / 2
  opupp <- 1 - oplow
  if (is.null(prob_annotate)) {
    prob_annotate <- prob_outer
  }
  aplow <- (1 - prob_annotate) / 2
  apupp <- 1 - aplow
  iplow <- (1 - prob_inner) / 2
  ipupp <- 1 - iplow
  ps <- sort(unique(c(aplow, iplow, 0.5, ipupp, apupp,
                      seq(oplow, opupp, length.out = 100))))
  ap_idx <- which(ps == aplow | ps == apupp)
  ip_idx <- which(ps >= iplow & ps < ipupp)
  med_idx <- which(ps == 0.5)

  if (use_subset) {
    sdat <- c(spec, spec$subset()$train)
  } else {
    sdat <- c(spec, spec$full()$train)
  }

  # calc full prob dens for each par
  cfgs <- mskt_pars(sdat, pars)
  plots <- lapply(cfgs, function(cfg) {
    xs <- do.call(get(paste0('q', cfg$pr)), c(list(p = ps), cfg$hp))
    ys <- do.call(get(paste0('d', cfg$pr)), c(list(x = xs), cfg$hp))
    xs <- cfg$tr(xs)
    y_nm <- NULL
    # if (cfg$nm == 'italic(r)') {
    #   y_nm <- 'Prior density'
    # }
    dens <- data.frame(x = xs, y = ys, mark = FALSE, inner = FALSE, med = FALSE)
    dens$mark[ap_idx] <- TRUE
    dens$inner[ip_idx] <- TRUE
    dens$med[med_idx] <- TRUE
    gg <-
      ggplot2::ggplot(dens, ggplot2::aes(.data[['x']], .data[['y']])) +
      ggplot2::geom_area(ggplot2::aes(fill = 'inner'),
                         data = dens[dens$inner, ],
                         # fill = '#f1f9ff',
                         outline.type = 'lower') +
      ggplot2::geom_line(ggplot2::aes(colour = 'density')) +
      ggplot2::geom_segment(ggplot2::aes(x = .data[['x']],
                                         xend = .data[['x']],
                                         y = 0,
                                         yend = .data[['y']]),
                            colour = '#6497b1',
                            data = dens[dens$med, ],
                            linewidth = 1) +
      # trick a vertical legend entry
      ggplot2::geom_vline(ggplot2::aes(xintercept = .data[['x']],
                                       linetype = 'med'),
                          na.rm = TRUE,
                          data = data.frame(x = NA_real_)) +
      ggplot2::geom_point(ggplot2::aes(shape = 'mark'),
                          data = dens[dens$mark, ],
                          colour = '#cc0000',
                          size = 2) +
      ggplot2::scale_colour_manual(NULL,
                                   labels = paste0('Density function'),
                                   values = '#03386c') +
      ggplot2::scale_fill_manual(NULL,
                                 labels = paste0(100 * prob_inner, '% HDI'),
                                 values = '#f1f9ff') +
      ggplot2::scale_linetype_manual(NULL,
                                     labels = 'Median',
                                     values = 'solid') +
      ggplot2::scale_shape_manual(NULL,
                                  labels = paste0('Edges of ',
                                                  100 * prob_annotate, '% HDI'),
                                  values = 1) +
      ggplot2::guides(colour = ggplot2::guide_legend(order = 1),
                      linetype = ggplot2::guide_legend(order = 2,
                                                       override.aes = list(
                                                         linewidth = 1,
                                                         colour = '#6497b1'
                                                       )),
                      fill = ggplot2::guide_legend(order = 3),
                      shape = ggplot2::guide_legend(order = 4)) +
      ggplot2::labs(x = parse(text = cfg$nm),
                    y = y_nm,
                    subtitle = cfg$title)

    cfg$geoms(gg)
  })
  if (plot) {
    patchwork::wrap_plots(plots, guides = 'collect')
  } else {
    plots
  }
}

#' @keywords internal
qhnorm <- function(p, sd = 1) {
  stats::qnorm((p + 1) / 2, 0, sd)
}
#' @keywords internal
dhnorm <- function(x, sd = 1, log = FALSE) {
  ds <- 2 * stats::dnorm(x, 0, sd)
  if (log) {
    ds <- log(ds)
  }
  ds
}
#' @keywords internal
rhnorm <- function(n, sd = 1) {
  us <- stats::runif(n, 0.5, 1.0)
  xs <- stats::qnorm(us, 0, sd)
  xs
}
# local({
#   sd <- 1
#   curve(dhnorm(x, sd), 0, sd * 3)
#   lines(density(rhnorm(1e5, sd)), col = 2)
# })
# qphi <- function(p, rate = 1) {
#   (rate / log(1 / (1 - p))) ^ 2
# }
# dphi <- function(x, rate = 1) {
#   (rate / (2 * x ^ (3 / 2))) * exp(-rate / sqrt(x))
# }

# local({
#   N <- 1e6
#   lambda <- 0.5
#
#   # Probability density function of Y
#   fy <- function(y) (lambda / (2 * y ^ (3 / 2))) * exp(-lambda / sqrt(y))
#
#   # Quantile function of Y
#   Qy <- function(p) (lambda / log(1 / (1 - p))) ^ 2
#
#   y <- seq(0, 2, length.out = 100)
#
#   fys <- fy(y)
#   y_approx <- rexp(N, lambda)
#   fys_approx = density(y_approx, from = 0, to = 1 / sqrt(2))
#   fys_approx$x <- fys_approx$x ^ -2
#
#   # approximate by sampling
#   plot(fys_approx, xlim = c(0, 100))
#   abline(v = median(y_approx))
#
#   # exact by derivation
#   lines(y, fys, col = 2)
#   abline(v = Qy(0.5), col = 2)
# })

#' Check prior predictions of the modskurt mean cover all possible shapes
#'
#' @param spec the modskurt model specification returned by `mskt_spec`
#' @param use_subset whether to use specified subset or all data
#' @param n_draws number of draws from the joint prior to use
#'
#' @importFrom rlang .data
#'
#' @return a ggplot object
#' @export
check_prior_mskt <- function(spec,
                             use_subset = TRUE,
                             n_draws = 50) {
  # the mean function has at most 6 parameters
  # so we could look at combinating different prob outcomes
  # from the prior
  # all median => 6 lines
  # 90% quantiles + median => 6 * 6 * 6 = 216 lines
  # vs 216 random draws which would show distributions better?

  nms <- attr(spec, 'nms')

  if (use_subset) {
    sdat <- c(spec, spec$subset()$train)
  } else {
    sdat <- c(spec, spec$full()$train)
  }

  pars <- mskt_pars(sdat, setdiff(attr(spec, 'pars'), c('kap', 'g0', 'g1')))
  draws <- lapply(setNames(names(pars), names(pars)), function(nm) {
    cfg <- pars[[nm]]
    raws <- do.call(get(paste0('r', cfg$pr)),
                    c(list(n = n_draws), cfg$hp))
    log_prob <- do.call(get(paste0('d', cfg$pr)),
                        c(list(x = raws, log = TRUE), cfg$hp))
    list(val = cfg$tr(raws), lp = log_prob)
  })
  vals <- do.call(cbind, lapply(draws, '[[', 'val'))
  prob <- exp(rowSums(do.call(cbind, lapply(draws, '[[', 'lp'))))
  mus <- do.call(rbind, lapply(seq_len(n_draws), function(n) {
    data.frame(draw = n,
               x = sdat$xrep,
               mu = do.call(mskt, c(list(x = sdat$xrep), vals[n, ])),
               jp = prob[n])
  }))
  ggplot2::ggplot(mus) +
    ggplot2::geom_line(ggplot2::aes(x = .data[['x']],
                                    y = .data[['mu']],
                                    group = .data[['draw']],
                                    alpha = .data[['jp']]),
                       colour = '#DD3300') +
    ggplot2::scale_alpha_identity(guide = guide_legend('Relative joint\nprior probability')) +
    # scale_x_continuous(expand = c(0.01, 0.01)) +
    # scale_y_continuous(sec.axis = sec_axis(~ .), limits = c(0, sdat$y_max)) +
    ggplot2::labs(y = nms$y, x = nms$x)
}

#' Check prior distribution statistics for potential misspecification clues
#'
#' @param spec the modskurt model specification returned by `mskt_spec`
#' @param use_subset whether to use specified subset or all data
#' @param n_draws number of draws from the joint prior to use
#'
#' @importFrom rlang .data
#'
#' @export
check_prior_dist <- function(spec,
                             use_subset = TRUE,
                             n_draws = 50) {
  nms <- attr(spec, 'nms')

  if (use_subset) {
    sdat <- c(spec, spec$subset()$train)
  } else {
    sdat <- c(spec, spec$full()$train)
  }

  pars <- mskt_pars(sdat, attr(spec, 'pars'))
  draws <- lapply(setNames(names(pars), names(pars)), function(nm) {
    cfg <- pars[[nm]]
    raws <- do.call(get(paste0('r', cfg$pr)), c(list(n = n_draws), cfg$hp))
    cfg$tr(raws)
  })
  vals <- do.call(cbind, draws)
  yrep <- do.call(rbind, lapply(seq_len(n_draws), function(n) {
    y <- do.call(mskt_sim_nb, c(list(x = sdat$xrep), vals[n, ]))$y
    cbind(round(vals[n, , drop = FALSE], 2),
          data.frame(name = 'rep',
                     draw = n,
                     max_y = max(y),
                     prop_0 = mean(y == 0),
                     sd_y = stats::sd(y)))
  }))
  yobs <- data.frame(name = 'obs',
                     draw = 0,
                     max_y = max(sdat$y),
                     prop_0 = mean(sdat$y == 0),
                     sd_y = stats::sd(sdat$y))
  plots <- lapply(c('max_y', 'prop_0', 'sd_y'), function(stat) {
    ggplot2::ggplot() +
      ggplot2::geom_histogram(ggplot2::aes(x = .data[[stat]],
                                           fill = .data[['name']]),
                              bins = n_draws / 2,
                              na.rm = TRUE,
                              colour = '#88CC88',
                              data = yrep) +
      ggplot2::geom_vline(ggplot2::aes(xintercept = .data[[stat]],
                                       colour = .data[['name']]),
                          linewidth = 1,
                          data = yobs) +
      ggplot2::scale_colour_manual(nms$y,
                                   values = c('obs' = '#005500'),
                                   guide = ggplot2::guide_legend(order = 1),
                                   labels = 'Observed') +
      ggplot2::scale_fill_manual(NULL,
                                 values = c('rep' = '#BBFFAA'),
                                 guide = ggplot2::guide_legend(order = 2),
                                 labels = 'Replicated') +
      ggplot2::labs(y = NULL) +
      ggplot2::theme(legend.title = ggplot2::element_text(
        margin = ggplot2::margin(b = 6, unit = 'pt')),
        legend.margin = ggplot2::margin(b = -6, unit = 'pt'))
  })
  # yrep$name <- NULL
  # print(yrep)
  (plots[[1]] + ggplot2::labs(x = 'max(y)', y = '% of reps')) +
    (plots[[2]] + ggplot2::labs(x = 'sum(y = 0) / N') + ggplot2::xlim(0, 1)) +
    (plots[[3]] + ggplot2::labs(x = 'sd(y)')) +
    patchwork::plot_layout(guides = 'collect')
}

#' Check prior predictions of excess zero probabilities cover all possible relationships
#'
#' @param spec the modskurt model specification returned by `mskt_spec`
#' @param use_subset whether to use specified subset or all data
#' @param n_draws number of draws from the joint prior to use
#'
#' Plot prob zero for median mu_n? or just for a bell curve to show how it works?
#'
#' @importFrom rlang .data
#'
#' @export
check_prior_zero <- function(spec,
                             use_subset = TRUE,
                             n_draws = 50) {
  nms <- attr(spec, 'nms')

  if (use_subset) {
    sdat <- c(spec, spec$subset()$train)
  } else {
    sdat <- c(spec, spec$full()$train)
  }
  H <- 100
  pars <- mskt_pars(sdat, intersect(attr(spec, 'pars'), c('kap', 'g0', 'g1')))
  draws <- lapply(setNames(names(pars), names(pars)), function(nm) {
    cfg <- pars[[nm]]
    raws <- do.call(get(paste0('r', cfg$pr)),
                    c(list(n = n_draws), cfg$hp))
    log_prob <- do.call(get(paste0('d', cfg$pr)),
                        c(list(x = raws, log = TRUE), cfg$hp))
    # don't need to transform any dist pars
    list(val = raws, lp = log_prob)
  })
  vals <- do.call(cbind, lapply(draws, '[[', 'val'))
  mus <- data.frame(x = sdat$xrep,
                    mu = mskt(sdat$xrep,
                              H,
                              min(sdat$xrep) + (
                                max(sdat$xrep) - min(sdat$xrep)
                              ) / 2,
                              sdat$x_pos_range / 3,
                              d = 2))
  mus$mu_std <- mus$mu / H
  probs <- list(
    nb = exp(draws$kap$lp),
    zi = exp(draws$g0$lp + draws$g1$lp),
    zinbl = exp(draws$kap$lp + draws$g0$lp + draws$g1$lp)
  )
  probs <- lapply(probs, \(ps) 0.5 * ps / max(ps))
  zis <- do.call(rbind, lapply(seq_len(n_draws), function(n) {
    data.frame(draw = n,
               x = mus$x,
               nb = do.call(prob_nb_zero, c(list(mu = mus$mu), vals[n, 1])),
               nb_jp = probs$nb[n],
               zi = do.call(zilink, c(list(mu = mus$mu, H = H), vals[n, -1])),
               zi_jp = probs$zi[n],
               zinbl = do.call(prob_zinbl_zero,
                               c(list(mu = mus$mu, H = H), vals[n, ])),
               zinbl_jp = probs$zinbl[n])
  }))
  plots <- lapply(c('nb', 'zi', 'zinbl'), function(stat) {
    ggplot2::ggplot() +
      ggplot2::geom_line(ggplot2::aes(x = .data[['x']],
                                      y = .data[['mu_std']]),
                         data = mus) +
      ggplot2::geom_line(ggplot2::aes(x = .data[['x']],
                                      y = .data[[stat]],
                                      group = .data[['draw']],
                                      alpha = .data[[paste0(stat, '_jp')]]),
                         colour = '#0033DD',
                         data = zis) +
      ggplot2::scale_alpha_identity(guide = 'none') +
      ggplot2::labs(y = NULL, x = NULL) +
      ggplot2::ylim(0, 1)
  })
  (plots[[1]] + ggplot2::labs(y = 'Prior Pr(y = 0)',
                              subtitle = 'Negative Binomial')) +
    (plots[[2]] + ggplot2::labs(x = paste0('Grid of ', nms$x),
                                subtitle = expression('Excess zero,'~pi))) +
    (plots[[3]] + ggplot2::labs(subtitle = 'ZINBL'))
}

# check_prior_zero(spec) # zi link

mskt_pars <- function(sdat, pars) {
  with(sdat, {
    list(
      'H' = list(title = substitute(
        italic(H) %~% 'Beta' * group('(', list(alpha, beta), ')') %.%
          'max' ~ italic(y),
        list(alpha = hp_H[1], beta = hp_H[2])),
        nm = 'italic(H)',
        pr = 'beta',
        hp = list(shape1 = hp_H[1], shape2 = hp_H[2]),
        tr = \(z) z * y_max,
        geoms = \(gg) {
          gg +
            ggplot2::geom_vline(xintercept = y_max,
                                alpha = 0.5,
                                lty = 2) +
            ggplot2::xlim(0, y_max)# +
          # geom_text(aes(label = label),
          #           parse = TRUE,
          #           hjust = 0,
          #           vjust = -0.25,
          #           angle = 270,
          #           family = 'serif',
          #           size = 5,
          #           data = ~ data.frame(
          #             label = '"max " ~ italic(y)',
          #             x = y_max,
          #             y = max(.x$y)))
        }),
      'm' = list(title = substitute(
        italic(m) %~% 'Beta' * group('(', list(alpha, beta), ')') %.%
          'range' ~ italic(x)[italic(y)>0] + 'min' ~ italic(x)[italic(y)>0],
        list(alpha = hp_m[1], beta = hp_m[2])),
        nm = 'italic(m)',
        pr = 'beta',
        hp = list(shape1 = hp_m[1], shape2 = hp_m[2]),
        tr = \(z) z * x_pos_range + x_pos_min,
        geoms = \(gg) {
          gg +
            ggplot2::geom_vline(xintercept = x_pos_min,
                                alpha = 0.5,
                                lty = 2) +
            ggplot2::geom_vline(xintercept = x_pos_range + x_pos_min,
                                alpha = 0.5,
                                lty = 2)# +
          # geom_text(aes(label = label),
          #           parse = TRUE,
          #           hjust = 0,
          #           vjust = -0.25,
          #           angle = 90,
          #           family = 'serif',
          #           size = 5,
          #           data = data.frame(
          #             label = '"min" ~ italic(x)[italic(y)>0]',
          #             x = x_pos_min,
          #             y = 0)) +
          # geom_text(aes(label = label),
          #           parse = TRUE,
          #           hjust = 1,
          #           vjust = -0.25,
          #           angle = 270,
          #           family = 'serif',
          #           size = 5,
          #           data = data.frame(
          #             label = '"max" ~ italic(x)[italic(y)>0]',
          #             x = x_pos_range + x_pos_min,
          #             y = 0))
        }),
      's' = list(title = substitute(
        'log' ~ italic(s) %~% 'Normal' * group('(', list(mu, sigma), ')') %.%
          'range' ~ italic(x)[italic(y)>0],
        list(mu = hp_s[1], sigma = hp_s[2])),
        nm = 'italic(s)',
        pr = 'norm',
        hp = list(mean = hp_s[1], sd = hp_s[2]),
        tr = \(z) exp(z) * x_pos_range,
        geoms = \(gg) {
          gg +
            ggplot2::geom_vline(xintercept = x_pos_range,
                                alpha = 0.5,
                                lty = 2)# +
          # geom_text(aes(label = label),
          #           parse = TRUE,
          #           hjust = -0.1,
          #           family = 'serif',
          #           size = 5,
          #           data = ~ data.frame(
          #             label = '"range" ~ italic(x)[italic(y)>0]',
          #             x = x_pos_range,
          #             y = max(.x$y)))
        }),
      'r' = list(title = substitute(
        italic(r) %~% 'Beta' * group('(', list(alpha, beta), ')'),
        list(alpha = hp_r[1], beta = hp_r[2])),
        nm = 'italic(r)',
        pr = 'beta',
        hp = list(shape1 = hp_r[1], shape2 = hp_r[2]),
        tr = identity,
        geoms = identity),
      'd' = list(title = substitute(
        'log' ~ italic(d) %~% 'N' * group('(', list(mu, sigma), ')'),
        list(mu = hp_d[1], sigma = hp_d[2])),
        nm = 'italic(d)',
        pr = 'norm',
        hp = list(mean = hp_d[1], sd = hp_d[2]),
        tr = \(z) exp(z),
        geoms = identity),
      'p' = list(title = substitute(
        italic(p) %~% 'Beta' * group('(', list(alpha, beta), ')') %.%
          1.95 + 0.05,
        list(alpha = hp_p[1], beta = hp_p[2])),
        nm = 'italic(p)',
        pr = 'beta',
        hp = list(shape1 = hp_p[1], shape2 = hp_p[2]),
        tr = \(z) z * 1.95 + 0.05,
        geoms = identity),
      'kap' = list(title = substitute(
        kappa == 1 / sqrt(phi) %~% 'Exp' * group('(', list(lambda), ')'),
        list(lambda = hp_kap)),
        nm = 'kappa',
        pr = 'exp',
        hp = list(rate = hp_kap),
        tr = identity,
        geoms = identity),
      # 'phi' = list(title = expression(
      #   phi %~% kappa ^ -2),
      #   pr = 'phi',
      #   hp = list(rate = hp_kap),
      #   tr = identity,
      #   geoms = \(gg) gg + xlim(0, 1)),
      'g0' = list(title = substitute(
        gamma[0] %~% 'Normal' * group('(', list(mu, sigma), ')'),
        list(mu = hp_g0[1], sigma = hp_g0[2])),
        nm = 'gamma[0]',
        # TODO: do this in geoms only...?
        # nm = 'logit^-1~gamma[0]',
        # tr = stats::plogis,
        pr = 'norm',
        hp = list(mean = hp_g0[1], sd = hp_g0[2]),
        tr = identity,
        geoms = identity),
      # geoms = \(gg) {
      #   gg +
      #     scale_x_continuous(expression(logit^-1~gamma[0]),
      #                        # limits = c(0.01, 0.99),
      #                        # breaks = c(0.01, seq(0.1, 0.9, 0.1), 0.99),
      #                        # labels = c('~0', seq(0.1, 0.9, 0.1), '~1'),
      #                        trans = scales::logit_trans())
      # }),
      'g1' = list(title = substitute(
        gamma[1] %~% 'Half-Normal' * group('(', list(sigma), ')'),
        list(sigma = hp_g1)),
        nm = 'gamma[1]',
        pr = 'hnorm',
        hp = list(sd = hp_g1),
        tr = identity,
        geoms = \(gg) {
          gg +
            xlim(0, NA)
        }))[pars]
  })
}

