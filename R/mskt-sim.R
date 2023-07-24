#' Simulate negative binomial data with a modskurt mean curve
#'
#' @param x vector of environmental gradient values
#' @param n_per_x number of y replicates to simulate for each x value
#' @param H height
#' @param m middle
#' @param s spread
#' @param r assymetry
#' @param d flatness
#' @param p persistence
#' @param kap negative binomial dispersion
#' @param g0 constant rate of zero-inflation
#' @param g1 decay of zero-inflation as mean increases
#'
#' @return data.frame with simulated mu, zi, and y for each x
#' @export
#'
mskt_sim_nb <- function(x,
                        n_per_x = 1,
                        H,
                        m,
                        s,
                        r = 0.5,
                        d = 0,
                        p = 1,
                        kap,
                        g0 = -Inf,
                        g1 = 0) {
  N <- length(x) * n_per_x
  yreps <- data.frame(x = x, mu = mskt(x, H, m, s, r, d, p))
  yreps$zi <- zilink(yreps$mu, H, g0, g1)
  yreps$y <- as.integer(stats::rnbinom(N,
                                       mu = yreps$mu,
                                       size = kappa_to_phi(kap)))
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
#' @keywords internal
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
  fma(H,
      pow(fma(r, exp(fma(invp, xms / r, mdop)),
              fma(1 - r, exp(fma(-invp, xms / (1 - r), mdop)),
                  - exp(mdop) + 1)),
          -p),
      .Machine$double.eps)
}
fma <- function(a, b, c) a * b + c
pow <- function(a, b) a ^ b

#' ZINBL zero-inflation link function
#' @param mu mean negative binomial abundance (real number vector)
#' @param H height of curve at mode (real number > 0)
#' @param g0 log-odds probability of excess zero when mu = 0
#' @param g1 rate at which log-odds probability of excess zero decreases as mu increases
#' @return vector of excess zero probabilities (real number vector in [0, 1])
#' @keywords internal
zilink <- function(mu, H, g0, g1) {
  stats::plogis(fma(-g1 / H, mu, g0))
}

#' @keywords internal
prob_nb_zero <- function(mu, kap) {
  stats::dnbinom(0, mu = mu, size = kappa_to_phi(kap))
}

#' @keywords internal
prob_zinbl_zero <- function(mu, H, kap, g0 = -Inf, g1 = 0) {
  nb_0s <- prob_nb_zero(mu, kap)
  zi_0s <- stats::plogis(fma(-g1 / H, mu, g0))
  zinbl_0 <- zi_0s + (1 - zi_0s) * nb_0s
  zinbl_0
}
