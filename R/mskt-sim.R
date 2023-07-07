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
  yreps <-
    data.frame(x = x, mu = mskt(x, H, m, s, r, d, p))
  yreps$zi <- stats::plogis(stats::qlogis(g0) - g1 * yreps$mu)
  yreps$y <- as.integer(stats::rnbinom(N,
                                       mu = yreps$mu,
                                       size = kappa_to_phi(kappa)))
  yreps$y[stats::runif(N) < yreps$zi] <- 0L
  yreps
}


#' Log-Modskurt (mskt) unimodal mean function
#'
#' @param x predictor variable (real number vector)
#' @param H height of curve at mode (real number > 0)
#' @param m mode (real number)
#' @param s scale (real number > 0)
#' @param r asymmetry (real number in (0, 1))
#' @param d flatness (real number >= 0)
#' @param p tail exaggeration (real number in (0, 2))
#' @return vector of log-scale mean values (real number vector)
mskt <- function(x, H, m, s, r = 0.5, d = 0, p = 1) {
  xms <- (x - m) / s
  # this underflows to -inf on log then 0 on natural
  # mu <- log(H) - p * log(
  #   r * exp((xms / r - d) / p) +
  #     (1 - r) * exp(- (xms / (1 - r) + d) / p) -
  #     exp(- d / p) + 1
  # )
  # exp(mu)
  # this version doesn't underflow and matches stan
  invp <- 1 / p
  mdop <- - d * invp
  lmu <- log(fma(H,
                 pow(fma(r,
                         exp(fma(invp, xms / r, mdop)),
                         fma(1 - r,
                             exp(fma(-invp, xms / (1 - r), mdop)),
                             exp(- mdop) + 1)),
                     -p),
                 .Machine$double.eps))
  exp(lmu)
}
fma <- function(a, b, c) a * b + c
pow <- function(a, b) a ^ b
