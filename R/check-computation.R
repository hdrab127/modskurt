#' Check posterior computation is acceptable
#'
#' @param fit modskurt_fit object
#'
#' @return a list with ...
#' @export
#'
check_computation <- function(fit) {
  spec <- attr(fit, 'spec')
  summ <- fit$summary(attr(spec, 'pars'))
  summ$num_samples <- fit$metadata()$iter_sampling * fit$num_chains()
  summ$rhat <- format(summ$rhat, digits = 3)
  summ$ess_bulk <- paste0(format(round(summ$ess_bulk)), ' (',
                          format(round(summ$ess_bulk / summ$num_samples, 1),
                                 nsmall = 1), ')')
  summ$ess_tail <- paste0(format(round(summ$ess_tail)), ' (',
                          format(round(summ$ess_tail / summ$num_samples, 1),
                                 nsmall = 1), ')')

  summ$rhat <- vctrs::new_vctr(summ$rhat,
                               class = 'rhat',
                               inherit_base_type = TRUE)
  summ$ess_bulk <- vctrs::new_vctr(summ$ess_bulk,
                                   class = 'ess',
                                   inherit_base_type = TRUE)
  summ$ess_tail <- vctrs::new_vctr(summ$ess_tail,
                                   class = 'ess',
                                   inherit_base_type = TRUE)
  summ$num_samples <- NULL

  summ$mean <- NULL
  summ$sd <- NULL
  summ$mad <- NULL

  diags <-
    cbind(fit$time()$chains,
          fit$diagnostic_summary())

  # list(summary = summ, diagnostics = diags)
  print(summ)
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
