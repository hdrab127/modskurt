#' visual display of chain mixing and prior-data conflicts
#'
#' @param fit fitted model object from `mskt_fit`
#' @param pars parameters to plot densities for, default is all
#' @param prob_inner inner probability (high density interval) to highlight
#' @param prob_outer outer probability range to plot the density between
#' @param by_chain plot chains separately to identify multiple modes and mixing
#'   issues
#'
#' @importFrom rlang .data
#'
#' @export
check_post_dens <- function(fit,
                            pars = attr(spec, 'pars'),
                            prob_inner = 0.5,
                            prob_outer = 0.99,
                            by_chain = FALSE) {
  oplow <- (1 - prob_outer) / 2
  opupp <- 1 - oplow
  iplow <- (1 - prob_inner) / 2
  ipupp <- 1 - iplow
  ps <- sort(unique(c(c(iplow, 0.5, ipupp),
                      seq(oplow, opupp, length.out = 100))))
  spec <- attr(fit, 'spec')
  if (attr(fit, 'is_subset')) {
    sdat <- c(spec, spec$subset()$train)
  } else {
    sdat <- c(spec, spec$full()$train)
  }
  cfgs <- mskt_pars(sdat, pars)
  pars[pars %in% c('g0', 'g1')] <- paste0(pars[pars %in% c('g0', 'g1')], '[1]')
  names(cfgs)[names(cfgs) %in% c('g0', 'g1')] <-
    paste0(names(cfgs)[names(cfgs) %in% c('g0', 'g1')], '[1]')
  if (by_chain) {
    draws <- fit$draws(pars, format = 'array')
    draws <- do.call(rbind, lapply(seq_len(fit$num_chains()), function(i) {
      cdat <- as.data.frame(draws[, i, , drop = TRUE])
      cdat$chain <- i
      cdat
    }))
  } else {
    draws <- as.data.frame(fit$draws(pars, format = 'matrix'))
    draws$chain <- 1
  }

  plots <- lapply(names(cfgs), function(nm) {
    cfg <- cfgs[[nm]]
    # prior
    pr_xs <- do.call(get(paste0('q', cfg$pr)), c(list(p = ps), cfg$hp))
    pr_ys <- do.call(get(paste0('d', cfg$pr)), c(list(x = pr_xs), cfg$hp))
    pr_xs <- cfg$tr(pr_xs)
    pr_dens <- data.frame(x = pr_xs, y = pr_ys)
    # normalise for comparison
    pr_dens$y <- pr_dens$y / max(pr_dens$y)

    # posterior
    xs <- draws[[nm]]
    dens <- do.call(rbind, lapply(unique(draws$chain), function(i) {
      idx <- draws$chain == i
      dens <- with(stats::density(xs[idx], bw = 'SJ', adjust = 3, cut = 0),
                   data.frame(x, y))
      dens$Fy <- cumsum(dens$y) / sum(dens$y)
      dens <- dens[dens$Fy >= oplow & dens$Fy < opupp, ]
      dens$inner <- dens$Fy >= iplow & dens$Fy < ipupp
      dens$med <- FALSE
      dens$med[which.min(abs(dens$Fy - 0.5))] <- TRUE
      dens$chain <- i
      # could density return proper pdf?
      dens$y <- dens$y / max(dens$y)
      dens
    }))
    dens$chain <- factor(dens$chain)

    gg <-
      ggplot2::ggplot(dens, ggplot2::aes(.data[['x']], .data[['y']]))
    if (by_chain) {
      gg <-
        gg +
        ggplot2::geom_line(data = pr_dens,
                           linetype = 'dashed',
                           colour = '#AAAAAA') +
        ggplot2::geom_line(ggplot2::aes(group = .data[['chain']],
                                      colour = .data[['chain']])) +
        ggplot2::geom_segment(ggplot2::aes(x = .data[['x']],
                                         xend = .data[['x']],
                                         y = 0,
                                         yend = .data[['y']],
                                         group = .data[['chain']],
                                         colour = .data[['chain']]),
                             alpha = 0.5,
                             data = dens[dens$med, ]) +
        ggplot2::labs(colour = 'Chain #')
    } else {
      gg <-
        gg +
        ggplot2::geom_line(data = pr_dens, colour = 'red') +
        ggplot2::geom_area(fill = 'white',
                          alpha = 0.9,
                          outline.type = 'lower') +
        ggplot2::geom_area(data = dens[dens$inner, ],
                          fill = '#e1f1fc',
                          outline.type = 'lower') +
        ggplot2::geom_line(colour = '#03386c') +
        ggplot2::geom_segment(ggplot2::aes(x = .data[['x']],
                                          xend = .data[['x']],
                                          y = 0,
                                          yend = .data[['y']]),
                             colour = '#6497b1',
                             data = dens[dens$med, ],
                             linewidth = 1)
    }
    gg <-
      gg +
      ggplot2::labs(x = parse(text = cfg$nm),
                   y = NULL,
                   # TODO: change to terms like peak mean abundance (H) etc.
                   subtitle = cfg$title)

    cfg$geoms(gg)
  })
  patchwork::wrap_plots(plots, guides = 'collect')
}

#' Check prediction calibration using discrete PIT histograms
#'
#' @param fit fitted model object from `mskt_fit`
#' @param ndraws number of draws from the joint posterior to use
#'
#' @importFrom rlang .data
#'
#' @export
check_post_calibration <- function(fit, ndraws = 50) {
  # both use randomised
  # subset uses test data
  # full uses leave-one-out

  # for y ~ ZINBL with cdf F
  # PIT score = F(y-1) + U * [F(y) - F(y-1)]
  # this can be histogrammed
  # or its density estimated and corrected then plotted over 0,1

  # so we need mu_rep, phi, g0, g1 for the cdf at each x_rep
  # then interpolate x_rep, mu_rep to match x_test and get y_test
  # then calc PITs of y_test and compare to uniform

  spec <- attr(fit, 'spec')
  is_subset <- attr(fit, 'is_subset')
  if (is_subset) {
    d <- spec$subset()$all
  } else {
    d <- spec$full()$all
  }
  # use cpue (will do nothing if no effort specd)
  d$y <- d$y / d$eff
  mus <- mskt_predict(fit, ndraws)

  # TODO: interpolate or change subset fit to include x_test in xreps
  d$xrep <- vapply(d$x,
                   function(x) spec$xrep[which.min(abs(x - spec$xrep))],
                   double(1))
  ys <- merge(mus, d, by.x = 'x', by.y = 'xrep')
  ys$Fy <- pzinb2(ys$y, ys$mu, ys$kap, ys$zi)
  ys$Fym1 <- 0
  pos_idx <- ys$y > 0
  ys$Fym1[ys$y > 0] <- pzinb2(ys$y[pos_idx] - 1,
                              ys$mu[pos_idx],
                              ys$kap[pos_idx],
                              ys$zi[pos_idx])
  ys$PIT <- ys$Fym1 + stats::runif(nrow(ys)) * (ys$Fy - ys$Fym1)

  fPIT <- tapply(ys, ys$set, function(dset) {
    # boundary corrected density of dicrete PIT scores
    with(stats::density(dset$PIT, n = 199, from = -0.49, to = 1.49), {
      lows <- x < 0
      upps <- x > 1
      x[lows] <- abs(x[lows])
      x[upps] <- 2 - x[upps]
      x <- round(x, 2)
      dens <- data.frame(x, y)
      dens <- dens[order(dens$x), ]
      dens <- data.frame(x = unique(dens$x),
                         y = tapply(dens$y, dens$x, sum))
      dens$y[1] <- dens$y[2]
      dens$y[nrow(dens)] <- dens$y[nrow(dens) - 1]
      dens
    })
  })

  # ndraws from uniform to compare
  U <- do.call(rbind, lapply(seq_len(ndraws), function(i) {
    data.frame(draw = i,
               table(sample(0:100, 1000, replace = TRUE)))
  }))
  colnames(U) <- c('draw', 'x', 'y')
  U$x <- (as.double(U$x) - 1) / 100
  U$y <- U$y / 10

  ggplot2::ggplot(U, ggplot2::aes(.data[['x']], .data[['y']])) +
    ggplot2::geom_line(ggplot2::aes(group = .data[['draw']], colour = 'Uniform'),
                      alpha = 0.2) +
    ggplot2::geom_line(ggplot2::aes(colour = 'Training'),
                      linewidth = 1.2,
                      data = fPIT$train) +
    ggplot2::geom_line(ggplot2::aes(colour = 'Test'),
                      linewidth = 1.2,
                      data = fPIT$test) +
    ggplot2::lims(y = c(0, NA), x = c(0, 1)) +
    ggplot2::scale_colour_manual(values = c('#BBBBBB', '#88CC88', '#005500'),
                                breaks = c('Uniform', 'Training', 'Test')) +
    ggplot2::labs(x = 'Discrete PIT score',
                 y = 'Relative frequency',
                 colour = 'Source')
}

#' Check for posterior misspecification by assessing influence data points
#'
#' @param fit fitted model object from `mskt_fit`
#'
#' @importFrom rlang .data
#'
#' @export
check_post_influencers <- function(fit) {
  spec <- attr(fit, 'spec')
  nms <- attr(spec, 'nms')
  is_subset <- attr(fit, 'is_subset')
  if (is_subset) {
    d <- spec$subset()$all
    d <- d[d$set == 'train', ]
  } else {
    d <- spec$full()$all
  }
  # use cpue (will do nothing if no effort specd)
  d$y <- d$y / d$eff
  d$set <- NULL

  pks <- cbind(data.frame(fit$loo()[['pointwise']]), d)
  pks$obs <- 1:nrow(pks)
  pks$is_zero <- 'Absent'
  pks$is_zero[pks$y > 0] <- 'Present'
  names(pks)[names(pks) == 'influence_pareto_k'] <- 'pareto_khat'
  high_pk <- pks[pks$pareto_khat > 0.5, ]
  if (nrow(high_pk) > 0) {
    high_pk[nms$y] <- high_pk$y
    high_pk[nms$x] <- high_pk$x
    high_pk[nms$eff] <- high_pk$eff
    print(
      high_pk[, which(names(high_pk) %in% c(nms, 'pareto_khat'))],
      digits = 2
    )
  }

  ggplot2::ggplot(pks, ggplot2::aes(.data[['x']], .data[['pareto_khat']])) +
    ggplot2::geom_point(ggplot2::aes(colour = .data[['is_zero']],
                                   shape = .data[['is_zero']])) +
    ggplot2::geom_text(ggplot2::aes(label = .data[['obs']],
                                  colour = .data[['is_zero']]),
                      nudge_y = 0.05,
                      show.legend = FALSE,
                      data = high_pk) +
    ggplot2::scale_colour_manual(values = c(`Present` = '#6497b1',
                                           `Absent` = 'black')) +
    ggplot2::scale_shape_manual(values = c(`Present` = 3,
                                          `Absent` = 1)) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = 0.0),
                       lty = 3,
                       colour = '#d1e1ec') +
    ggplot2::geom_hline(ggplot2::aes(yintercept = 0.5),
                       lty = 4,
                       colour = '#7c0c00') +
    ggplot2::geom_hline(ggplot2::aes(yintercept = 0.7),
                       lty = 2,
                       colour = '#7c0c00') +
    ggplot2::geom_hline(ggplot2::aes(yintercept = 1.0),
                       lty = 1,
                       colour = '#7c0c00') +
    ggplot2::labs(x = nms$x,
                 y = expression('Pareto'~hat(k)),
                 colour = nms$y,
                 shape = nms$y)
}
