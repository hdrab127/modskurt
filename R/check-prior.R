#' Check marginal prior densities are encoded correctly
#'
#' @param spec the modskurt model specification returned by `mskt_spec`
#' @param train_prop the proportion of data to train the posterior on
#' @param train_seed a seed for the training sample of proportion `train_prop`
#' @param prob_inner inner probability (high density interval) to highlight
#' @param prob_outer outer probability range to plot the density between
#'
#' @return a patchwork of ggplots
#' @export
#'
check_prior_dens <- function(spec,
                             train_prop = 1,
                             train_seed = NULL,
                             prob_inner = 0.5,
                             prob_outer = 0.99) {
  # outer prob (op), inner prob (ip)
  oplow <- (1 - prob_outer) / 2
  opupp <- 1 - oplow
  iplow <- (1 - prob_inner) / 2
  ipupp <- 1 - iplow
  ps <- sort(unique(c(c(iplow, 0.5, ipupp),
                      seq(oplow, opupp, length.out = 100))))
  ip_idx <- which(ps >= iplow & ps < ipupp)
  med_idx <- which(ps == 0.5)

  # TODO: this should prob use same seed as fit
  sdat <- c(spec, spec$data(train_prop, train_seed)$train)

  # calc full prob dens for each par
  plots <-
    with(sdat, {
      pars <- list(
        'italic(H)' = list(title = expression(
          italic(H) %~% 'Beta' * group('(', list(alpha, beta), ')') %.%
            'max' ~ italic(y)),
          pr = 'beta',
          hp = list(shape1 = hp_H[1], shape2 = hp_H[2]),
          tr = \(z) z * y_max,
          geoms = \(gg) {
            gg +
              geom_vline(xintercept = y_max,
                         alpha = 0.5,
                         lty = 2) +
              xlim(0, y_max)# +
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
        'italic(m)' = list(title = expression(
          italic(m) %~% 'Beta' * group('(', list(alpha, beta), ')') %.%
            'range' ~ italic(x)[italic(y)>0] + 'min' ~ italic(x)[italic(y)>0]),
          pr = 'beta',
          hp = list(shape1 = hp_m[1], shape2 = hp_m[2]),
          tr = \(z) z * x_pos_range + x_pos_min,
          geoms = \(gg) {
            gg +
              geom_vline(xintercept = x_pos_min,
                         alpha = 0.5,
                         lty = 2) +
              geom_vline(xintercept = x_pos_range + x_pos_min,
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
        'italic(s)' = list(title = expression(
          'log' ~ italic(s) %~% 'Normal' * group('(', list(mu, sigma), ')') %.%
            'range' ~ italic(x)[italic(y)>0]),
          pr = 'norm',
          hp = list(mean = hp_s[1], sd = hp_s[2]),
          tr = \(z) exp(z) * x_pos_range,
          geoms = \(gg) {
            gg +
              geom_vline(xintercept = x_pos_range,
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
        'italic(r)' = list(title = expression(
          italic(r) %~% 'Beta' * group('(', list(alpha, beta), ')')),
          pr = 'beta',
          hp = list(shape1 = hp_r[1], shape2 = hp_r[2]),
          tr = identity,
          geoms = identity),
        'italic(d)' = list(title = expression(
          'log' ~ italic(d) %~% 'Normal' * group('(', list(mu, sigma), ')')),
          pr = 'norm',
          hp = list(mean = hp_d[1], sd = hp_d[2]),
          tr = \(z) exp(z),
          geoms = identity),
        'italic(p)'   = list(title = expression(
          italic(p) %~% 'Beta' * group('(', list(alpha, beta), ')') %.%
            1.99 + 0.01),
          pr = 'beta',
          hp = list(shape1 = hp_p[1], shape2 = hp_p[2]),
          tr = \(z) z * 1.95 + 0.05,
          geoms = identity),
        'kappa' = list(title = expression(
          kappa == 1 / sqrt(phi) %~% 'Exp' * group('(', list(lambda), ')')),
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
        'gamma[0]' = list(title = expression(
          gamma[0] %~% 'Normal' * group('(', list(mu, sigma), ')')),
          pr = 'norm',
          hp = list(mean = hp_g0[1], sd = hp_g0[2]),
          tr = identity,
          geoms = identity),
        'gamma[1]' = list(title = expression(
          gamma[1] %~% 'Half-Normal' * group('(', list(sigma), ')')),
          pr = 'hnorm',
          hp = list(sd = hp_g1),
          tr = \(z) z,
          geoms = \(gg) {
            gg +
              xlim(0, NA)
          }))
      if (spec$use_zi == 0) {
        pars[['gamma[0]']] <- NULL
        pars[['gamma[1]']] <- NULL
      }
      lapply(names(pars), function(nm) {
        cfg <- pars[[nm]]
        xs <- do.call(get(paste0('q', cfg$pr)), c(list(p = ps), cfg$hp))
        ys <- do.call(get(paste0('d', cfg$pr)), c(list(x = xs), cfg$hp))
        xs <- cfg$tr(xs)
        y_nm <- NULL
        if (nm == 'italic(r)') {
          y_nm <- 'Prior density'
        }
        dens <-
          list(outer = data.frame(x = xs, y = ys),
               inner = data.frame(x = xs[ip_idx], y = ys[ip_idx]),
               med = data.frame(x = xs[med_idx], y = ys[med_idx]))
        gg <-
          ggplot(dens$outer, aes(x, y)) +
          geom_line(colour = '#03386c') +
          geom_area(data = dens$inner,
                    fill = '#d1e1ec',
                    alpha = 0.5,
                    outline.type = 'lower') +
          geom_segment(aes(x = x, xend = x, y = 0, yend = y),
                       colour = '#6497b1',
                       data = dens$med,
                       lwd = 1) +
          labs(x = parse(text = nm),
               y = y_nm,
               subtitle = cfg$title) +
          theme_classic()

        cfg$geoms(gg)
      })
    })
  patchwork::wrap_plots(plots, guides = 'collect')
}

qhnorm <- function(p, sd = 1) {
  stats::qnorm((p + 1) / 2, 0, sd)
}
dhnorm <- function(x, sd = 1) {
  2 * stats::dnorm(x, 0, sd)
}
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

# check_prior_mskt(spec) # preds (cover all possible)

# check_prior_dist(spec) # ymax (severe tests)

# check_prior_zero(spec) # zi link
