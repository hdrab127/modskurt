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
#' @importFrom stats setNames
#' @importFrom rlang .data
#'
#' @return a ggplot object
#' @export
#'
abundance_dist <- function(fit,
                           seed = NULL,
                           ndraws = 50,
                           summaries = c('mean', 'q90'),
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
  is_subset <- attr(fit, 'is_subset')
  if (is_subset) {
    d <- spec$subset()$all
  } else {
    d <- spec$full()$all
  }
  # use cpue (will do nothing if no effort specd)
  d$y <- d$y / d$eff

  mus <- mskt_predict(fit, ndraws, seed, include_zero_inflation)


  gg <-
    ggplot2::ggplot(d, ggplot2::aes(.data$x, .data$y))
  if (!is_subset) {
    gg <-
      gg +
      ggplot2::geom_point(alpha = point_alpha)
  } else {
    gg <-
      gg +
      ggplot2::geom_point(ggplot2::aes(shape = .data$set), alpha = point_alpha) +
      ggplot2::scale_shape_manual('Dataset',
                                  values = c('train' = 19, 'test' = 4),
                                  breaks = c('train', 'test'),
                                  labels = c('Train', 'Test')) +
      ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(
        alpha = 0.75,
        order = 2
      )))
  }
  # add summaries
  # first expectation or nb mean
  if ('mean' %in% summaries) {
    gg <-
      gg +
      ggplot2::geom_line(ggplot2::aes(y = .data$mu,
                                      group = .data$draw,
                                      colour = .data$name),
                         alpha = line_alpha,
                         data = mus)
  }
  if (!include_zero_inflation | dist == 'nb') {
    mu_nm <- 'NB mean'
  } else {
    mu_nm <- 'E[y]'
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
    qcols <- setNames(qcols, names(qps))
    if ('mean' %in% summaries) {
      qcols <- c(setNames('#22CC99', mu_nm), qcols)
    }
    for (qnm in names(qps)) {
      gg <-
        gg +
        ggplot2::geom_line(ggplot2::aes(y = .data$y,
                                        group = .data$draw,
                                        colour = .data$name),
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
                                   labels = names(qcols))
    # labels = parse(text = names(qcols)))
  } else if ('mean' %in% summaries) {
    gg <-
      gg +
      ggplot2::scale_colour_manual(values = '#22CC99',
                                   labels = mu_nm)
  } else {
    stop('at least one summary must be specified')
  }
  gg +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(
      alpha = 0.75,
      linewidth = 1,
      order = 1
    ))) +
    ggplot2::labs(x = nms$x,
                  y = nms$y,
                  colour = substitute(y%~%distnm,
                                      list(distnm = toupper(dist))))
}

#' @keywords internal
mskt_predict <- function(fit, ndraws = NULL, seed = NULL, include_zero_inflation = TRUE) {
  spec <- attr(fit, 'spec')
  dist <- attr(spec, 'dist')
  niter <- fit$metadata()$iter_sampling
  nchains <- fit$num_chains()
  npost <- niter * nchains

  mus <- fit$draws('mu_rep', format = 'list')
  kap <- as.double(fit$draws('kap', format = 'matrix'))
  zis <- fit$draws('zi_rep', format = 'list')
  mus <- do.call(rbind, lapply(seq_along(mus), \(chain) {
    draw_offset <- (chain - 1) * niter
    do.call(rbind, lapply(seq_along(mus[[chain]]), \(rep) {
      draws <-
        data.frame(chain = chain,
                   x = spec$xrep[[rep]],
                   mu = mus[[chain]][[rep]],
                   zi = zis[[chain]][[rep]])
      draws$draw <- draw_offset + 1:nrow(draws)
      draws
    }))
  }))
  if (dist == 'nb' | !include_zero_inflation) {
    mus$zi <- 0
  } else {
    # zinbl mean/expectation
    mus$mu <- (1 - mus$zi) * mus$mu
  }
  mus$kap <- kap
  mus <- mus[, c('chain', 'draw', 'x', 'mu', 'kap', 'zi')]

  # subset draws if needed
  if (!missing(seed) & !is.null(seed)) {
    set.seed(seed)
  }
  if (!is.null(ndraws)) {
    mus <- mus[mus$draw %in% sample(npost, ndraws), ]
  }
  mus <- mus[order(mus$x), ]
  mus$name <- 'mu'
  mus
}

#' Assess abundance range along x using distributional measures
#'
#' @param fit modskurt_fit object
#' @param seed optional seed for ndraws sample
#' @param ndraws optional subset of posterior draws
#' @param capture_pct percent of total abundance density (most_dens) or percent
#'   tolerance from highest abundance summary (high_zone), based on "based_on"
#'   summary of abundance
#' @param using_range "ADL" for highest abundance density limit, or
#'   "HAZ" for highest abundance zone, see details
#' @param based_on which distribution summary to use for density or high zone,
#'   one of "mean", "median", or any quantile proportion as "q" + pecentage
#'   (e.g. "q90" = 90\% quantile)
#' @param region one-sided threshold ("left", "right") or central range
#'   ("centre")
#' @param include_zero_inflation whether to use expectations and quantiles from
#'   nb or zinbl
#' @param plotted plot on the abundance distribution
#' @param point_alpha opacity of data points, default 0.1 is good for ~ 1000 obs
#' @param line_alpha opacity of summary lines, default 0.1 is good for ~ 50
#'   draws
#' @param range_colour colour of range geometries
#'
#' @importFrom stats sd
#' @importFrom rlang .data
#'
#' @return a matrix of ranges, or a ggplot object if plotted is TRUE
#' @export
#'
abundance_range <- function(fit,
                            seed = NULL,
                            ndraws = 50,
                            capture_pct = 80,
                            using_range = 'HAZ',
                            based_on = 'median',
                            region = 'centre',
                            include_zero_inflation = FALSE,
                            plotted = FALSE,
                            point_alpha = 0.2,
                            line_alpha = 0.02,
                            range_colour = '#22CC99') {
  # calcs intervals and optionally plots over abundance_dist
  spec <- attr(fit, 'spec')
  nms <- attr(spec, 'nms')
  dist <- attr(spec, 'dist')
  prop <- capture_pct / 100
  is_subset <- attr(fit, 'is_subset')
  if (is_subset) {
    d <- spec$subset()$all
  } else {
    d <- spec$full()$all
  }
  # use cpue (will do nothing if no effort specd)
  d$y <- d$y / d$eff
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
  range_fn <- get(paste0('range_', using_range))
  ranges <-
    do.call(rbind, tapply(mus, mus$draw, function(draw) {
      rg <-
        range_fn(x = draw$x,
                 y = draw$dy,
                 prop = prop,
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
  dy_nm <- paste0(toupper(substring(based_on, 0, 1)),
                  substring(based_on, 2))
  if (dy_nm == 'Mean' & (!include_zero_inflation | dist == 'nb')) {
    dy_nm <- 'NB mean'
  } else {
    dy_nm <- 'E[y]'
  }
  dy_nm <- gsub('^Q(\\d*)$', 'Q[\\1]', dy_nm)
    dy_console <- dy_nm
  if (grepl('^Q', dy_nm)) {
    dy_console <-
      paste0(as.double(gsub('[^0-9]*', '', dy_nm)) / 100, ' prob. quantile of')
  }
  cat('\nSpecies range (see x.avg row) calculated as the ',
      switch(using_range,
             'HAZ' = paste0('region of x where the ',
                                  dy_console, ' abundance (',
                                  dy_nm, ') is within ',
                                  capture_pct, '% of the highest value of ',
                                  dy_nm, ' along x'),
             'ADL' = paste0('smallest region of x where the density of ',
                                  dy_console, ' abundance (',
                                  dy_nm, ') is equal to ',
                                  capture_pct, '% of the total density of ',
                                  dy_nm, ' along x')),
      ' (averaged across posterior draws):\n', sep = '')

  if (plotted) {
    gg <-
      ggplot2::ggplot(d, ggplot2::aes(.data$x, .data$y))
    if (!is_subset) {
      gg <-
        gg +
        ggplot2::geom_point(alpha = point_alpha)
    } else {
      gg <-
        gg +
        ggplot2::geom_point(ggplot2::aes(shape = .data$set), alpha = point_alpha) +
        ggplot2::scale_shape_manual('Dataset',
                                    values = c('train' = 19, 'test' = 4),
                                    breaks = c('train', 'test'),
                                    labels = c('Train', 'Test')) +
        ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(
          alpha = 0.75
        )))
    }
    gg <-
      gg +
      ggplot2::geom_line(ggplot2::aes(y = .data$dy,
                                      group = .data$draw,
                                      colour = dy_nm),
                         alpha = line_alpha,
                         data = mus)

    cont_nm <- paste0(capture_pct, '% ',
                      switch(using_range,
                             'HAZ' = 'HAZ',
                             'ADL' = 'ADL'))
    xidx <- ranges$name == 'x'
    xavgidx <- row.names(aggs) == 'x.avg'
    if (region == 'centre') {
      range_vlines <-
        data.frame(x = c(aggs$left[xavgidx], aggs$right[xavgidx]),
                   leg = 'Average')
      # ranges$centre[xidx],
      range_points <-
        data.frame(x = c(ranges$left[xidx], ranges$right[xidx]),
                   y = c(ranges$left[!xidx], ranges$right[!xidx]),
                   group = rep(ranges$draw[xidx], 2),
                   leg = 'Draws')
    } else if (region == 'left') {
      range_vlines <-
        data.frame(x = c(aggs$right[xavgidx]), leg = 'Average')
      range_points <-
        data.frame(x = c(ranges$right[xidx]),
                   y = c(ranges$right[!xidx]),
                   group = ranges$draw[xidx],
                   leg = 'Draws')
    } else if (region == 'right') {
      range_vlines <-
        data.frame(x = c(aggs$left[xavgidx]), leg = 'Average')
      range_points <-
        data.frame(x = c(ranges$left[xidx]),
                   y = c(ranges$left[!xidx]),
                   group = ranges$draw[xidx],
                   leg = 'Draws')
    }
    gg <-
      gg +
      ggplot2::geom_point(ggplot2::aes(.data$x,
                                       .data$y,
                                       group = .data$group,
                                       alpha = .data$leg),
                          size = 2,
                          shape = 4,
                          colour = range_colour,
                          data = range_points) +
      ggplot2::geom_vline(ggplot2::aes(xintercept = .data$x,
                                       alpha = .data$leg),
                          linetype = 'dashed',
                          linewidth = 1,
                          colour = range_colour,
                          data = range_vlines) +
      ggplot2::scale_colour_manual(values = range_colour) +
      ggplot2::scale_alpha_manual(cont_nm,
                                  values = c(point_alpha * 2, 1),
                                  breaks = c('Draws', 'Average')) +
      ggplot2::guides(alpha = ggplot2::guide_legend(override.aes = list(
        linetype = c('blank', 'dashed'),
        shape = c(19, NA_integer_)
      ), order = 2), colour = ggplot2::guide_legend(override.aes = list(
        alpha = 0.75,
        linewidth = 1
      ), order = 1))

    gg <-
      gg +
      ggplot2::labs(x = nms$x,
                    y = nms$y,
                    colour = 'Summary')#substitute(y%~%distnm,
                                      #  list(distnm = toupper(dist))))
    print(aggs)
    gg
  } else {
    aggs
  }
}

#' @importFrom utils tail
#' @importFrom stats weighted.mean
#' @keywords internal
range_ADL <- function(x, y, prop, region) {
  # the smallest region of x that has cumulative y summary >= some % of total y
  # summary, i.e. highest density interval
  # calc the cumulative sum of y along x then find quantiles
  # relies on data being sorted by x
  Fy <- cumsum(y)
  Fy <- Fy / max(Fy)
  alpha <- (1 - prop) / 2
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
  lx <- weighted.mean(x[lidx], lwt, na.rm = TRUE)
  ly <- weighted.mean(y[lidx], lwt, na.rm = TRUE)
  midx <- c(tail(which(Fy <= 0.5), 1), which(Fy > 0.5)[1])
  mwt <- 1 / abs(Fy[midx] - 0.5)
  mwt[!is.finite(mwt)] <- 1
  mx <- weighted.mean(x[midx], mwt, na.rm = TRUE)
  my <- weighted.mean(y[midx], mwt, na.rm = TRUE)
  uidx <- c(tail(which(Fy <= qps[2]), 1), which(Fy > qps[2])[1])
  uwt <- 1 / abs(Fy[uidx] - qps[2])
  uwt[!is.finite(uwt)] <- 1
  ux <- weighted.mean(x[uidx], uwt, na.rm = TRUE)
  uy <- weighted.mean(y[uidx], uwt, na.rm = TRUE)
  data.frame(name = c('x', 'y'),
             left = c(lx, ly),
             centre = c(mx, my),
             right = c(ux, uy),
             stringsAsFactors = FALSE)
}

#' @importFrom stats weighted.mean
#' @keywords internal
range_HAZ <- function(x, y, prop, region) {
  # the total region that has y summary >= some % of y summary
  # e.g. could be the range of x where mean(y[x]) is at least 50% of max(mean(y))
  b <- range(which(y >= max(y) * (1 - prop)))
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
