#' Check posterior computation is acceptable
#'
#' @param fit fitted model object from `mskt_fit`
#' @param show_info print information about model spec and posterior sample
#' @param hide_stats which summary stats to hide
#'
#' @return a list with ...
#' @export
#'
check_computation <- function(fit,
                              show_info = TRUE,
                              hide_stats = c()) {
  spec <- attr(fit, 'spec')
  dist <- attr(spec, 'dist')
  shape <- attr(spec, 'shape')
  if (attr(fit, 'is_subset')) {
    d <- spec$subset()
  } else {
    d <- spec$full()
  }
  niter <- fit$metadata()$iter_sampling
  nchains <- fit$num_chains()
  npost <- niter * nchains
  if (show_info) {
    cat('spec: ', dist, '[Hms', shape, ']',
        ' using ', sum(d$all$set == 'train'),
        ' obs out of ', nrow(d$all), ' (', d$prop * 100, '% sample)',
        '\n',
        sep = '')
    cat('post: ',
        nchains, ' chains each with ', niter, ' draws',
        ' (', npost, ' total)',
        '\n\n',
        sep = '')
  }
  summ <- fit$summary(attr(spec, 'pars'))
  summ$num_samples <- npost
  summ$rhat <- format(summ$rhat, digits = 3)
  summ$ess_bulk <- paste0(format(round(summ$ess_bulk)), ' (',
                          format(round(summ$ess_bulk / summ$num_samples, 1),
                                 nsmall = 1), ')')
  summ$ess_tail <- paste0(format(round(summ$ess_tail)), ' (',
                          format(round(summ$ess_tail / summ$num_samples, 1),
                                 nsmall = 1), ')')

  summ$rhat <- vctrs::new_vctr(summ$rhat, class = 'rhat')
  summ$ess_bulk <- vctrs::new_vctr(summ$ess_bulk, class = 'ess')
  summ$ess_tail <- vctrs::new_vctr(summ$ess_tail, class = 'ess')
  summ$num_samples <- NULL
  if (length(hide_stats) > 0) {
    summ <- summ[, -which(names(summ) %in% hide_stats)]
  }

  # summ$mean <- NULL
  # summ$sd <- NULL
  # summ$mad <- NULL
  # TODO: if one chain fails then below errs
  diags <-
    cbind(fit$time()$chains,
          fit$diagnostic_summary())

  list(summary = summ, diagnostics = diags)
}

#' @export
#' @aliases pillar_shaft
#' @importFrom pillar pillar_shaft
#' @method pillar_shaft rhat
#' @keywords internal
pillar_shaft.rhat <- function(x, ...) {
  highlight <- crayon::bgYellow
  xd <- vctrs::vec_data(x)
  num <- as.double(xd)
  xd[num >= 1.05] <- highlight(xd[num >= 1.05])
  pillar::new_pillar_shaft_simple(xd, align = 'right')
}
#' @export
#' @aliases pillar_shaft
#' @importFrom pillar pillar_shaft
#' @method pillar_shaft ess
#' @keywords internal
pillar_shaft.ess <- function(x, ...) {
  highlight <- crayon::bgYellow
  xd <- vctrs::vec_data(x)
  pct <- as.double(gsub('^.*\\((.*)\\)$', '\\1', xd))
  xd[pct <= 0.1] <- highlight(xd[pct <= 0.1])
  pillar::new_pillar_shaft_simple(xd, align = 'right')
}
