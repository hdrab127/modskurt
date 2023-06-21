#' The Negative Binomial Probability Mass Function (alternative parameterisation)
#'
#' @param x vector of (non-negative integer) quantiles
#' @param mu mean
#' @param kappa dispersion (inverse square root size)
#' @param log logical; if TRUE, probabilities p are given as log(p)
#'
#' @return the (log) probability mass at x, given mu and kappa
#' @export
dnbinom2 <- function(x, mu, kappa, log = FALSE) {
  if (length(log) != 1) {
    stop("log should be length 1.")
  }
  phi <- kappa_to_phi(kappa)
  stats::dnbinom(x = x, mu = mu, size = phi, log = log)
}

#' The Negative Binomial Quantile Function (alternative parameterisation)
#'
#' @param p vector of probabilities
#' @param mu mean
#' @param kappa dispersion (inverse square root size)
#' @param lower.tail logical; if TRUE (default), probablities are \eqn{P[X <= x]}, otherwise, \eqn{[X > x]}
#' @param log.p logical; if TRUE, probabilities are given as log(p)
#'
#' @return A vector of quantiles, each of which correspond to a probability in p
#' @export
qnbinom2 <- function(p, mu, kappa, lower.tail = TRUE, log.p = FALSE) {
  if (length(lower.tail) != 1) {
    stop("lower.tail should be length 1.")
  }
  if (length(log.p) != 1) {
    stop("log.p should be length 1.")
  }
  phi <- kappa_to_phi(kappa)
  stats::qnbinom(p = p,
                 size = phi,
                 mu = mu,
                 lower.tail = lower.tail,
                 log.p = log.p)
}


#' Zero-inflated Negative Binomial Probability Mass
#'
#' @param x vector of (non-negative integer) values
#' @param mu mean of the distribution
#' @param kappa dispersion (inverse square root size)
#' @param zi vector of (real lying in \[0, 1\]) zero-inflation parameters
#' @param log logical; if TRUE, probabilities, p, are given as log(p)
#'
#' @return The (log) probability mass at x, given mu, kappa and zi
#' @export
dzinb2 <- function(x, mu, kappa, zi, log = FALSE) {
  if (length(log) != 1) {
    stop("log should be length 1.")
  }
  if(any(zi < 0 | zi > 1)) {
    stop(("argument zi contains values that lie outside the  interval [0, 1]"))
  }
  probabilityMass = zi * (x == 0) +
    (1 - zi) * dnbinom2(x = x, mu = mu, kappa = kappa)

  if (log) {
    return(log(probabilityMass))
  } else {
    return(probabilityMass)
  }
}

#' Zero-inflated Negative Binomial Quantile Function
#'
#' @param p vector of quantiles
#' @param mu mean of the distribution
#' @param kappa dispersion (inverse square root size)
#' @param zi vector of (in \eqn{[0, 1]}) zero-inflation parameters
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[X <= x]}, otherwise, \eqn{P[X > x]}
#' @param log.p logical; if TRUE, probabilities p. This doesn't affect zi
#'
#' @return A vector of quantiles, each of which coincide with the respective probability in p
#' @export
qzinb2 <- function(p, mu, kappa, zi, lower.tail = TRUE, log.p = FALSE) {
  if (length(lower.tail) != 1) {
    stop("lower.tail should be length 1.")
  }
  if (length(log.p) != 1) {
    stop("log.p should be length 1.")
  }
  if (any(zi < 0 | zi > 1)) {
    stop("argument zi contains values that lie outside the  interval [0, 1]")
  }
  p_linear <- p
  if (log.p) {
    p_linear <- exp(p)
  }
  if (any(p_linear < 0 | p_linear > 1)) {
    stop("argument p contains values that represent probabilities outside the interval [0, 1]")
  }
  p_lower <-  p_linear
  if (! lower.tail) {
    p_lower <- 1 - p_linear
  }

  output_length <- max(length(p), length(mu), length(kappa), length(zi))
  q <- numeric(output_length)
  for (i in 1:output_length) {
    p_lower_i <- p_lower[[(i - 1) %% length(p_lower) + 1]]
    pi_i <- zi[[(i - 1) %% length(zi) + 1]]
    mu_i <- mu[[(i - 1) %% length(mu) + 1]]
    kap_i <- kappa[[(i - 1) %% length(kappa) + 1]]

    if (p_lower_i <= pi_i) {
      q[[i]] <- 0
    } else {
      q[[i]] <- qnbinom2(p = (p_lower_i - pi_i) / (1.0 - pi_i),
                         mu = mu_i,
                         kappa = kap_i,
                         lower.tail = TRUE,
                         log.p = FALSE)
    }
  }

  q
}

#' TODO: document
kappa_to_phi <- function(kappa) {
  # see ...
  1 / kappa ^ 2
}
