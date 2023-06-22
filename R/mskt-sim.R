#' Simulate negative binomial data with a modskurt mean curve
#'
#' @param x vector of environmental gradient values
#' @param H height
#' @param m middle
#' @param s spread
#' @param r assymetry
#' @param d flatness
#' @param p persistence
#' @param kappa negative binomial dispersion
#' @param g0 constant rate of zero-inflation
#' @param g1 decay of zero-inflation as mean increases
#'
#' @return data.frame with simulated mu, zi, and y for each x
#' @export
#'
mskt_sim_nb <- function(x,
                        H,
                        m,
                        s,
                        r = 0.5,
                        d = 0,
                        p = 1,
                        kappa,
                        g0 = 0,
                        g1 = 0) {
  N <- length(x)
  xms <- (x - m) / s
  yreps <-
    data.frame(x = x,
               mu = exp(log(H) - p * log(
                 r * exp((xms / r - d) / p) +
                   (1 - r) * exp(- (xms / (1 - r) + d) / p) -
                   exp(- d / p) + 1
               )))
  yreps$zi <- stats::plogis(stats::qlogis(g0) - g1 * yreps$mu)
  yreps$y <- as.integer(stats::rnbinom(N,
                                       mu = yreps$mu,
                                       size = kappa_to_phi(kappa)))
  yreps$y[stats::runif(N) < yreps$zi] <- 0L
  yreps
}
