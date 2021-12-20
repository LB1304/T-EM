## Authors: Luca Brusa, Francesco Bartolucci, Fulvia Pennoni ##

complk <- function (S, yv, piv, Pi, Phi, k) {

  sS = dim(S)
  ns = sS[1]
  TT = sS[2]
  if (length(sS) == 2)
    r = 1
  else r = sS[3]
  if (r == 1) {
    if (is.matrix(S))
      S = array(S, c(dim(S), 1))
  }

  Psi = array(1, c(ns, k, TT))
  L = array(0, c(ns, k, TT))
  for (j in 1:r)
    Psi[, , 1] = Psi[, , 1] * Phi[S[, 1, j] + 1, , j]
  L[, , 1] = Psi[, , 1] %*% diag(piv)
  for (t in 2:TT) {
    for (j in 1:r)
      Psi[, , t] = Psi[, , t] * Phi[S[, t, j] + 1, , j]
    L[, , t] = Psi[, , t] * (L[, , t - 1] %*% Pi[, , t])
  }

  if (ns == 1)
    pv = sum(L[1, , TT])
  else
    pv = rowSums(L[, , TT])

  lk = sum(yv * log(pv))

  return(list(lk = lk, Psi = Psi, L = L, pv = pv))
}
