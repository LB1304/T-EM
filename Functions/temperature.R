## Authors: Luca Brusa, Francesco Bartolucci, Fulvia Pennoni ##

temperature <- function(h, algorithm, profile_pars = list(alpha = NULL, beta = NULL, ro = NULL, T0 = NULL)) {

  T0 <- profile_pars$T0
  ro <- profile_pars$ro
  alpha <- profile_pars$a
  beta <- profile_pars$b

  if (algorithm == 1) {
    if (is.null(T0) || is.null(ro)) {
      stop("Provide all the correct tempering parameters.")
    }
    temp <- T0 * ro^(h-1)
  } else if (algorithm == 2) {
    if (is.null(alpha) || is.null(beta)) {
      stop("Provide all the correct tempering parameters.")
    }
    temp <- 1 + exp(beta - h/alpha)
  } else if (algorithm == 3) {
    if (is.null(T0) || is.null(ro) || is.null(alpha) || is.null(beta)) {
      stop("Provide all the correct tempering parameters.")
    }
    temp <- tanh(h/(2*ro)) + (T0 - beta * 2 * sqrt(2)/(3 * pi)) * alpha^(h/ro) + beta * phonTools::sinc(3 * pi/4 + h/ro, normalized = TRUE)
    if (temp < 1) temp <- 1
  } else {
    stop("Specify a correct tempering profile.")
  }

  return(temp)
}
