## Authors: Luca Brusa, Francesco Bartolucci, Fulvia Pennoni ##

LM_em <- function (S, yv, k, modBasic, tol, maxit, piv, Pi, Phi) {

  n = sum(yv)
  sS = dim(S)
  ns = sS[1]
  TT = sS[2]
  if (length(sS) == 2) {
    r = 1
    if (is.matrix(S))
      S = array(S, c(dim(S), 1))
  } else {
    r = sS[3]
  }
  Sv = matrix(S, ns * TT, r)
  bv = apply(Sv, 2, max)
  b = max(bv)

  out = complk(S, yv, piv, Pi, Phi, k)
  lk = out$lk
  Psi = out$Psi
  L = out$L
  pv = out$pv

  it = 0
  lko = lk - 10^10
  lkv = NULL

  while ((lk - lko)/abs(lk) > tol & it < maxit) {
    it = it + 1
    step_i <- LM_em_step(Sv = Sv, yv = yv, n = n, ns = ns, r = r, TT = TT, k = k,
                         bv = bv, b = b, piv = piv, Pi = Pi, Phi = Phi, Psi = Psi, pv = pv, L = L, mod = modBasic)
    piv <- step_i$piv
    Pi <- step_i$Pi
    Phi <- step_i$Phi
    V <- step_i$V

    lko = lk
    out = complk(S, yv, piv, Pi, Phi, k)
    lk = out$lk
    Psi = out$Psi
    L = out$L
    pv = out$pv
    lkv = c(lkv, lk)
  }

  np = (k - 1) + k * sum(bv)
  if (modBasic == 0)
    np = np + (TT - 1) * k * (k - 1)
  if (modBasic == 1)
    np = np + k * (k - 1)
  if (modBasic > 1)
    np = np + 2 * k * (k - 1)
  aic = -2 * lk + np * 2
  bic = -2 * lk + np * log(n)

  if (any(yv != 1))
    V = V/yv
  lk = as.vector(lk)
  dimnames(Pi) = list(state = 1:k, state = 1:k, time = 1:TT)
  dimnames(Phi) = list(category = 0:b, state = 1:k, item = 1:r)
  out = list(lk = lk, lkv = lkv, it = it, piv = piv, Pi = Pi, Phi = Phi, V = V, k = k, np = np, modBasic = modBasic, aic = aic, bic = bic)
  return(out)
}
