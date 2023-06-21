#' Distribution of abundance along the environmental gradient
#'
#' @param fit modskurt_fit object
#' @param seed optional seed for ndraws sample
#' @param ndraws optional subset of posterior draws
#' @param summaries character vector of distribution summarise, see details
#' @param include_zero_inflation whether to include zero-inflation in summaries
#' @param point_alpha opacity of data points, default 0.1 is good for ~ 1000 obs
#' @param line_alpha opacity of summary lines, default 0.1 is good for ~ 50 draws
#'
#' @return a ggplot object
#' @export
#'
abundance_dist <- function(fit,
                           seed = NULL,
                           ndraws = 50,
                           summaries = c('mean', 'median', 'q90'),
                           include_zero_inflation = TRUE,
                           point_alpha = 0.2,
                           line_alpha = 0.1) {
  # for some number of draws
  # obtain dist pars
  # calculate summaries
  # plot data
  # overlay spaghetti
  spec <- attr(fit, 'spec')
  nms <- attr(spec, 'nms')
  dist <- attr(spec, 'dist')
  d <- do.call(spec$data, attr(fit, 'train'))$all
  mus <- mskt_predict(fit, ndraws, seed, include_zero_inflation)

  # TODO: plot cpue instead?
  gg <-
    ggplot2::ggplot(d, ggplot2::aes(x, y))
  if (attr(fit, 'train')$prop == 1) {
    gg <-
      gg +
      ggplot2::geom_point(alpha = point_alpha)
  } else {
    gg <-
      gg +
      ggplot2::geom_point(ggplot2::aes(shape = set), alpha = point_alpha) +
      ggplot2::scale_shape_manual('Dataset',
                                  values = c('train' = 19, 'test' = 4),
                                  breaks = c('train', 'test'),
                                  labels = c('Train', 'Test')) +
      ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(
        alpha = 0.75
      )))
  }
  # add summaries
  # first expectation or nb mean
  if ('mean' %in% summaries) {
    gg <-
      gg +
      ggplot2::geom_line(ggplot2::aes(y = mu, group = draw, colour = name),
                         alpha = line_alpha,
                         data = mus)
  }
  # then quantiles
  qps <- summaries[grepl('^q', summaries)]
  if ('median' %in% summaries) {
    qps <- c('q50', qps)
  }
  if (length(qps) > 0) {
    qps <- setNames(as.double(gsub('^q', '', qps)) / 100,
                    qps)
    names(qps) <- paste0(gsub('q', 'Q[', names(qps)), ']')
    qcols <- grDevices::colorRampPalette(c('#88AAFF',
                                           '#EE77FF',
                                           '#DD0000'))(50)
    qcols <- c(qcols, rev(qcols[1:49]))
    qcols <- qcols[100 * qps]
    qcols <- c('Mean' = '#22CC99',
               setNames(qcols, names(qps)))
    for (qnm in names(qps)) {
      gg <-
        gg +
        ggplot2::geom_line(ggplot2::aes(y = y,
                                        group = draw,
                                        colour = name),
                           alpha = line_alpha,
                           data = data.frame(name = qnm,
                                             draw = mus$draw,
                                             x = mus$x,
                                             y = qzinb2(p = qps[[qnm]],
                                                        mu = mus$mu,
                                                        kappa = mus$kap,
                                                        zi = mus$zi)))
    }
    gg <-
      gg +
      ggplot2::scale_colour_manual(values = unname(qcols),
                                   labels = names(qcols)) +
      ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(
        alpha = 0.75,
        linewidth = 1
      )))
    # labels = parse(text = names(qcols)))
  }
  gg +
    ggplot2::labs(x = nms$x,
                  y = nms$y,
                  colour = substitute(y%~%distnm,
                                      list(distnm = toupper(dist))))
}

mskt_predict <- function(fit, ndraws, seed, include_zero_inflation) {
  spec <- attr(fit, 'spec')
  dist <- attr(spec, 'dist')
  niter <- fit$metadata()$iter_sampling
  nchains <- fit$num_chains()
  npost <- niter * nchains

  mus <- fit$draws('mu_rep', format = 'list')
  kap <- as.double(fit$draws('kap', format = 'matrix'))
  zis <- fit$draws('zi_rep', format = 'list')
  mus <- do.call(rbind, lapply(seq_along(mus), \(chain) {
    do.call(rbind, lapply(seq_along(mus[[chain]]), \(rep) {
      draws <-
        data.frame(chain = chain,
                   x = spec$xrep[[rep]],
                   mu = mus[[chain]][[rep]],
                   zi = zis[[chain]][[rep]])
      draws$draw <- chain * 1:nrow(draws)
      draws
    }))
  }))
  if (dist == 'nb' | !include_zero_inflation) {
    mus$zi <- 0
  }
  mus$kap <- kap

  # subset draws if needed
  if (!missing(seed) & !is.null(seed)) {
    set.seed(seed)
  }
  mus <- mus[mus$draw %in% sample(npost, ndraws),
             c('chain', 'draw', 'x', 'mu', 'kap', 'zi')]
  mus <- mus[order(mus$x), ]
  mus$name <- 'mu'
  mus
}

#' Title
#'
#' @param fit modskurt_fit object
#' @param seed optional seed for ndraws sample
#' @param ndraws optional subset of posterior draws
#' @param percent percent of total abundance density (based on summary) or
#'   percent tolerance from highest abundance summary (e.g. highest mean
#'   abundance)
#' @param containing "most_dens" for range of highest abundance density, or
#'   "high_zone" for range of x with abundance summary within highest possible,
#'   see details
#' @param region one-sided threshold ("left", "right") or central range
#'   ("centre")
#' @param based_on which distribution summary to use for density or high zone,
#'   one of "mean", "median", or any quantile proportion as "q" + pecentage
#'   (e.g. "q90" = 90\% quantile)
#' @param include_zero_inflation whether to use expectations and quantiles from
#'   nb or zinbl
#' @param plotted plot on the abundance distribution
#' @param point_alpha opacity of data points, default 0.1 is good for ~ 1000 obs
#' @param line_alpha opacity of summary lines, default 0.1 is good for ~ 50 draws
#'
#' @return
#' @export
#'
abundance_range <- function(fit,
                            seed = NULL,
                            ndraws = 50,
                            percent = 80,
                            containing = 'high_zone',
                            region = 'centre',
                            based_on = 'median',
                            include_zero_inflation = FALSE,
                            plotted = FALSE,
                            point_alpha = 0.2,
                            line_alpha = 0.1) {
  # calcs intervals and optionally plots over abundance_dist
  spec <- attr(fit, 'spec')
  nms <- attr(spec, 'nms')
  dist <- attr(spec, 'dist')
  pct <- percent / 100
  d <- do.call(spec$data, attr(fit, 'train'))$all
  mus <- mskt_predict(fit, ndraws, seed, include_zero_inflation)
  agg <- based_on
  if (agg == 'median') {
    agg <- 'q50'
  }
  if (grepl('^q', agg)) {
    mus$dy <- qzinb2(p = as.double(gsub('^q', '', agg)) / 100,
                     mu = mus$mu,
                     kappa = mus$kap,
                     zi = mus$zi)
  } else if (agg == 'mean') {
    mus$dy <- mus$mu
  }
  range_fn <- get(paste0('range_', containing))
  ranges <-
    do.call(rbind, tapply(mus, mus$draw, function(draw) {
      rg <-
        range_fn(x = draw$x,
                 y = draw$dy,
                 pct = pct,
                 region = region)
      rg$draw <- draw$draw[[1]]
      rg
    }))
  ranges$name[ranges$name == 'y'] <- based_on
  aggs <-
    do.call(rbind, tapply(ranges, ranges$name, function(rg) {
      as.data.frame(apply(rg[, c('left', 'centre', 'right')],
                          2,
                          function(js) {
                            c('avg' = mean(js),
                              'se' = sd(js) / sqrt(length(js)))
                          }))
    }))
  if (plotted) {
    gg <-
      ggplot2::ggplot(d, ggplot2::aes(x, y))
    if (attr(fit, 'train')$prop == 1) {
      gg <-
        gg +
        ggplot2::geom_point(alpha = point_alpha)
    } else {
      gg <-
        gg +
        ggplot2::geom_point(ggplot2::aes(shape = set), alpha = point_alpha) +
        ggplot2::scale_shape_manual('Dataset',
                                    values = c('train' = 19, 'test' = 4),
                                    breaks = c('train', 'test'),
                                    labels = c('Train', 'Test')) +
        ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(
          alpha = 0.75
        )))
    }
    mus$name <- based_on
    gg <-
      gg +
      ggplot2::geom_line(ggplot2::aes(y = dy, group = draw, colour = name),
                         alpha = line_alpha,
                         data = mus)
    xidx <- ranges$name == 'x'
    xavgidx <- row.names(aggs) == 'x.avg'
    gg <-
      gg +
      ggplot2::geom_point(ggplot2::aes(x,
                                       y,
                                       group = group,
                                       colour = colour),
                          alpha = point_alpha * 3,
                          size = 2,
                          data = data.frame(x = c(ranges$left[xidx],
                                                  # ranges$centre[xidx],
                                                  ranges$right[xidx]),
                                            y = c(ranges$left[!xidx],
                                                  # ranges$centre[!xidx],
                                                  ranges$right[!xidx]),
                                            group = rep(ranges$draw[xidx],
                                                        2),
                                            colour = 'range')) +
      ggplot2::geom_vline(ggplot2::aes(xintercept = x, colour = colour),
                          linetype = 'dashed',
                          data = data.frame(x = c(aggs$left[xavgidx],
                                                  aggs$right[xavgidx]),
                                            colour = 'range'))

    gg <-
      gg +
      ggplot2::labs(x = nms$x,
                    y = nms$y,
                    colour = substitute(y%~%distnm,
                                        list(distnm = toupper(dist))))
    print(gg)
  }
  aggs
}

range_most_dens <- function(x, y, pct, region) {
  # the smallest region of x that has cumulative y summary >= some % of total y
  # summary, i.e. highest density interval
  # calc the cumulative sum of y along x then find quantiles
  # relies on data being sorted by x
  Fy <- cumsum(y)
  Fy <- Fy / max(Fy)
  alpha <- (1 - pct) / 2
  qps <- c(alpha, 1 - alpha)
  if (region == 'left') {
    qps <- c(0, 1 - alpha * 2)
  } else if (region == 'right') {
    qps <- c(alpha * 2, 1)
  }
  # now find x closest to quantiles
  # take either side of alpha and weighted avg between them
  lidx <- c(tail(which(Fy <= qps[1]), 1), which(Fy > qps[1])[1])
  lwt <- 1 / abs(Fy[lidx] - qps[1])
  lwt[!is.finite(lwt)] <- 1
  lx <- weighted.mean(x[lidx], 1 / abs(Fy[lidx] - qps[1]), na.rm = TRUE)
  ly <- weighted.mean(y[lidx], 1 / abs(Fy[lidx] - qps[1]), na.rm = TRUE)
  midx <- c(tail(which(Fy <= 0.5), 1), which(Fy > 0.5)[1])
  mwt <- 1 / abs(Fy[midx] - 0.5)
  mwt[!is.finite(mwt)] <- 1
  mx <- weighted.mean(x[midx], 1 / abs(Fy[midx] - qps[1]), na.rm = TRUE)
  my <- weighted.mean(y[midx], 1 / abs(Fy[midx] - qps[1]), na.rm = TRUE)
  uidx <- c(tail(which(Fy <= qps[2]), 1), which(Fy > qps[2])[1])
  uwt <- 1 / abs(Fy[uidx] - qps[2])
  uwt[!is.finite(uwt)] <- 1
  ux <- weighted.mean(x[uidx], uwt, na.rm = TRUE)
  uy <- weighted.mean(x[uidx], uwt, na.rm = TRUE)
  data.frame(name = c('x', 'y'),
             left = c(lx, ly),
             centre = c(mx, my),
             right = c(ux, uy),
             stringsAsFactors = FALSE)
}

range_high_zone <- function(x, y, pct, region) {
  # the total region that has y summary >= some % of y summary
  # e.g. could be the range of x where mean(y[x]) is at least 50% of max(mean(y))
  b <- range(which(y >= max(y) * (1 - pct)))
  # weighted mean interpolates between ties
  lx <- weighted.mean(x[(b[1] - 1):b[1]], y[(b[1] - 1):b[1]], na.rm = TRUE)
  ly <- weighted.mean(y[(b[1] - 1):b[1]], y[(b[1] - 1):b[1]], na.rm = TRUE)
  mx <- weighted.mean(x[y == max(y)], y[y == max(y)], na.rm = TRUE)
  my <- max(y)
  ux <- weighted.mean(x[b[2]:(b[2] + 1)], y[b[2]:(b[2] + 1)], na.rm = TRUE)
  uy <- weighted.mean(y[b[2]:(b[2] + 1)], y[b[2]:(b[2] + 1)], na.rm = TRUE)
  if (region == 'left') {
    lx <- min(x)
    ly <- y[x == min(x)]
  } else if (region == 'right') {
    lx <- max(x)
    ly <- y[x == max(x)]
  }
  data.frame(name = c('x', 'y'),
             left = c(lx, ly),
             centre = c(mx, my),
             right = c(ux, uy),
             stringsAsFactors = FALSE)
}
