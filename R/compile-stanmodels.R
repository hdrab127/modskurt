#' Compile modskurt Stan Models
#'
#' @param force_recompile whether to recompile even if nothing has changed
#'
#' @export
#'
compile_stanmodels <- function(force_recompile = FALSE) {
  check <- cmdstanr::check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
  pkgdir <- path.package('modskurt')
  standir <- file.path(pkgdir, 'stan')
  # TODO: remove this from testing
  is_local <- !dir.exists(standir)
  if (is_local) {
    standir <- file.path(gsub('OneDrive - Massey University',
                              'msc',
                              pkgdir,
                              fixed = TRUE), 'inst/stan')
  }
  csms <- dir(standir, '.stan$', full.names = TRUE)
  lapply(csms, function(nm) {
    short <- gsub('^.*/', '', nm)
    if (!file.exists(gsub('stan$', 'exe', nm)) | force_recompile) {
      print(paste('Compiling now:', short))
    } else {
      print(paste('Already compiled:', short))
    }
    cmdstanr::cmdstan_model(
      stan_file = nm,
      include_paths = normalizePath(paste0(standir, '/blocks/')),
      stanc_options = list('O1'),
      pedantic = is_local,
      force_recompile = force_recompile
    )
  })
  invisible(csms)
}
