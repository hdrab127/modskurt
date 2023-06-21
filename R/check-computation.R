#' Check posterior computation is acceptable
#'
#' @param fit modskurt_fit object
#' @param show_info print information about model spec and posterior sample
#'
#' @return a list with ...
#' @export
#'
check_computation <- function(fit, show_info = TRUE) {
  spec <- attr(fit, 'spec')
  dist <- attr(spec, 'dist')
  shape <- attr(spec, 'shape')
  idx <- attr(fit, 'train')
  d <- do.call(spec$data, idx)$all
  niter <- fit$metadata()$iter_sampling
  nchains <- fit$num_chains()
  npost <- niter * nchains
  if (show_info) {
    cat('spec: ', dist, '[Hms', shape, ']',
        ' using ', sum(d$set == 'train'),
        ' obs out of ', nrow(d), ' (', idx$prop * 100, '% sample)',
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

  # summ$mean <- NULL
  # summ$sd <- NULL
  # summ$mad <- NULL

  diags <-
    cbind(fit$time()$chains,
          fit$diagnostic_summary())

  list(summary = summ, diagnostics = diags)
}

#' @export
#' @aliases pillar_shaft
#' @importFrom pillar pillar_shaft
#' @method pillar_shaft rhat
pillar_shaft.rhat <- function(x, ...) {
  purp <- crayon::make_style('#EE00FF')
  xd <- vctrs::vec_data(x)
  num <- as.double(xd)
  xd[num >= 1.05] <- purp(xd[num >= 1.05])
  pillar::new_pillar_shaft_simple(xd, align = 'right')
}
#' @export
#' @aliases pillar_shaft
#' @importFrom pillar pillar_shaft
#' @method pillar_shaft ess
pillar_shaft.ess <- function(x, ...) {
  purp <- crayon::make_style('#EE00FF')
  xd <- vctrs::vec_data(x)
  pct <- as.double(gsub('^.*\\((.*)\\)$', '\\1', xd))
  xd[pct <= 0.1] <- purp(xd[pct <= 0.1])
  pillar::new_pillar_shaft_simple(xd, align = 'right')
}
