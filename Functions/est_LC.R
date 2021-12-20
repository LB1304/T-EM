## Authors: Luca Brusa, Francesco Bartolucci, Fulvia Pennoni ##

est_LC <- function (S, yv, k, tol = 1e-8, maxit = 1e6, sv = NULL,
                    algorithm = 0, profile_pars = list(alpha = NULL, beta = NULL, ro = NULL, T0 = NULL)) {

  # Starting Values
  if (is.null(sv)) {
    sv <- LC_sv(S = S, k = k)
  }

  # EM algorithm
  if (algorithm == 0) {
    out <- LC_em(S = S, yv = yv, k = k, piv = sv$piv, Piv = sv$Piv, Psi = sv$Psi, Phi = sv$Phi, tol = tol, maxit = maxit)
  } else if (algorithm %in% c(1, 2, 3)) {
    out <- LC_tem(S = S, yv = yv, k = k, piv = sv$piv, Piv = sv$Piv, Psi = sv$Psi, Phi = sv$Phi, tol = tol, maxit = maxit,
                  algorithm = algorithm, profile_pars = profile_pars)
  } else {
    stop("Specify an available algorithm.")
  }

  return(out)
}
