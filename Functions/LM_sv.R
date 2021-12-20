## Authors: Luca Brusa, Francesco Bartolucci, Fulvia Pennoni ##

LM_sv <- function(S, k) {
  
  sS = dim(S)
  ns = sS[1]
  TT = sS[2]
  if (length(sS) == 2) {
    r = 1
    if (is.matrix(S)) 
      S = array(S, c(sS, 1))
  } else {
    r = sS[3]
  }
  Sv = matrix(S, ns * TT, r)
  bv = apply(Sv, 2, max)
  b = max(bv)
  
  piv = stats::runif(k)
  piv = piv/sum(piv)
  
  Pi = array(stats::runif(k^2 * TT), c(k, k, TT))
  for (t in 2:TT) Pi[, , t] = diag(1/rowSums(Pi[, , t])) %*% Pi[, , t]
  Pi[, , 1] = 0
  
  Phi = array(NA, c(b + 1, k, r))
  for (j in 1:r) {
    Phi[1:(bv[j] + 1), , j] = matrix(stats::runif((bv[j] + 1) * k), bv[j] + 1, k)
    for (u in 1:k) Phi[1:(bv[j] + 1), u, j] = Phi[1:(bv[j] + 1), u, j]/sum(Phi[1:(bv[j] + 1), u, j])
  }
  
  return(list(piv = piv, Pi = Pi, Phi = Phi))
}
