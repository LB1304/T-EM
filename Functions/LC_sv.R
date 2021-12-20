## Authors: Luca Brusa, Francesco Bartolucci, Fulvia Pennoni ##

LC_sv <- function(S, k) {
  
  c = max(S) + 1
  ns = nrow(S)
  r = ncol(S)
  
  piv = stats::runif(k)
  piv = piv/sum(piv)
  Piv = rep(1, ns) %o% piv
  Phi = array(stats::runif(c * r * k), c(c, r, k))
  for (u in 1:k) for (j in 1:r) Phi[, j, u] = Phi[, j, u]/sum(Phi[, j, u])
  Psi = matrix(1, ns, k)
  for (u in 1:k) {
    for (j in 1:r) {
      Psi[, u] = Psi[, u] * Phi[S[, j] + 1, j, u]
    }
  }
  
  return(list(piv = piv, Piv = Piv, Phi = Phi, Psi = Psi))
}
