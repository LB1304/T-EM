## Authors: Luca Brusa, Francesco Bartolucci, Fulvia Pennoni ##

LC_tem_step <- function (S, yv, C, ns, r, n, k, piv, Piv, Phi, Psi, temp) {

  Num <- (Piv * Psi)^(1/temp)
  Den <- rowSums(Num) %o% rep(1, k)
  freq <- yv %o% rep(1, k)
  V <- freq * Num/Den
  b <- colSums(V)

  piv <- b/n
  Piv = rep(1, ns) %o% piv
  a <- array(dim = c(C, r, k))
  for (u in 1:k) {
    for (c in 1:C) {
      for (j in 1:r) {
        ind <- (S[, j] == c - 1)
        a[c, j, u] <- sum(V[ind, u])
        Phi[c, j, u] <- a[c, j, u]/b[u]
      }
    }
  }
  Psi <- matrix(1, ns, k)
  for (j in 1:r) for (u in 1:k) Psi[, u] = Psi[, u] * Phi[S[, j] + 1, j, u]

  return(list(piv = piv, Piv = Piv, Psi = Psi, Phi = Phi))
}
