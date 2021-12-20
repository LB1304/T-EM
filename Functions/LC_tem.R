## Authors: Luca Brusa, Francesco Bartolucci, Fulvia Pennoni ##

LC_tem <- function(S, yv, k, piv, Piv, Psi, Phi, tol, maxit, algorithm, profile_pars) {

  lk <- sum(yv * log(rowSums(Psi * Piv)))
  lko <- lk - 10^10
  lkv <- vector(mode = "numeric")
  it <- 0

  C = max(S) + 1
  ns = nrow(S)
  r = ncol(S)
  n = sum(yv)

  while (((abs(lk - lko)/abs(lko) > tol) && it < maxit) || it < 2) {
    it = it + 1

    temp <- temperature(h = it, algorithm = algorithm, profile_pars = profile_pars)

    step_i <- LC_tem_step(S = S, yv = yv, C = C, ns = ns, r = r, n = n, k = k, piv = piv, Piv = Piv, Psi = Psi, Phi = Phi, temp = temp)
    piv <- step_i$piv
    Piv <- step_i$Piv
    Psi <- step_i$Psi
    Phi <- step_i$Phi
    lko <- lk
    lk <- sum(yv * log(rowSums(Psi * Piv)))
    lkv[it] <- lk
  }

  it = it + 1
  step_i <- LC_em_step(S = S, yv = yv, C = C, ns = ns, r = r, n = n, k = k, piv = piv, Piv = Piv, Psi = Psi, Phi = Phi)
  piv <- step_i$piv
  Piv <- step_i$Piv
  Psi <- step_i$Psi
  Phi <- step_i$Phi
  lko <- lk
  lk <- sum(yv * log(rowSums(Psi * Piv)))
  lkv[it] <- lk

  np = k * r * (C - 1) + k - 1
  aic = -2 * lk + 2 * np
  bic = -2 * lk + np * log(n)

  return(list(lk = lk, lkv = lkv, it = it, piv = piv, Piv = Piv, Psi = Psi, Phi = Phi, k = k, np = np, aic = aic, bic = bic,
              algorithm = algorithm, profile_pars = profile_pars))
}
