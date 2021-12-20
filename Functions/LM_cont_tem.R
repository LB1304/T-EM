## Authors: Luca Brusa, Francesco Bartolucci, Fulvia Pennoni ##

LM_cont_tem <- function (Y, k, modBasic, tol, maxit, piv, Pi, Mu, Si, algorithm, profile_pars) {

  sY = dim(Y)
  n = sY[1]
  TT = sY[2]
  if (length(sY) == 2) {
    r = 1
    if (is.matrix(Y))
      Y = array(Y, c(dim(Y), 1))
  } else {
    r = sY[3]
  }

  out = complk_cont(Y, piv, Pi, Mu, Si, k)
  lk = out$lk
  Psi = out$Psi
  L = out$L

  it = 0
  lko = lk - 10^10
  lkv = NULL

  while (abs((lk - lko)/max(lk, lko)) > tol & it < maxit) {
    it = it + 1

    temp <- temperature(h = it, algorithm = algorithm, profile_pars = profile_pars)

    step_i <- LM_cont_tem_step(Y = Y, n = n, r = r, TT = TT, k = k, L = L, Psi = Psi, piv = piv, Pi = Pi, Mu = Mu, Si = Si, modBasic = modBasic, temp = temp)

    piv <- step_i$piv
    Pi <- step_i$Pi
    Mu <- step_i$Mu
    Si <- step_i$Si
    V <- step_i$V

    lko = lk
    out = complk_cont(Y, piv, Pi, Mu, Si, k)
    lk = out$lk
    Psi = out$Psi
    L = out$L
    lkv = c(lkv, lk)
  }

  it <- it + 1
  step_i <- LM_cont_em_step(Y = Y, n = n, r = r, TT = TT, k = k, L = L, Psi = Psi, piv = piv, Pi = Pi, Mu = Mu, Si = Si, modBasic = modBasic)

  piv <- step_i$piv
  Pi <- step_i$Pi
  Mu <- step_i$Mu
  Si <- step_i$Si
  V <- step_i$V

  lko = lk
  out = complk_cont(Y, piv, Pi, Mu, Si, k)
  lk = out$lk
  Psi = out$Psi
  L = out$L
  lkv = c(lkv, lk)

  np = (k - 1) + k * r + r * (r + 1)/2
  if (modBasic == 0)
    np = np + (TT - 1) * k * (k - 1)
  if (modBasic == 1)
    np = np + k * (k - 1)
  if (modBasic > 1)
    np = np + 2 * k * (k - 1)
  aic = -2 * lk + np * 2
  bic = -2 * lk + np * log(n)

  lk = as.vector(lk)
  dimnames(Pi) = list(state = 1:k, state = 1:k, time = 1:TT)
  nameY <- dimnames(Y)[[3]]
  dimnames(Mu) <- list(nameY, state = 1:k)
  if (r == 1)
    colnames(Si) <- nameY
  else dimnames(Si) <- list(nameY, nameY)
  out = list(lk = lk, lkv = lkv, it = it, piv = piv, Pi = Pi, Mu = Mu, Si = Si, V = V, k = k, np = np, modBasic = modBasic, aic = aic, bic = bic,
             algorithm = algorithm, profile_pars = profile_pars)
  out$Y = Y

  return(out)
}
