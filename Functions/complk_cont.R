## Authors: Luca Brusa, Francesco Bartolucci, Fulvia Pennoni ##

complk_cont <- function (Y, piv, Pi, Mu, Si, k) {
  
  sY = dim(Y)
  n = as.integer(sY[1])
  TT = as.integer(sY[2])
  if (length(sY) == 2) {
    r = 1
  } else {
    r = sY[3]
  }
  r = as.integer(r)
  if (r == 1) {
    if (is.matrix(Y)) 
      Y = array(Y, c(dim(Y), 1))
  }
  
  Psi = array(1, c(n, k, TT))
  L = array(0, c(n, k, TT))
  
  for (u in 1:k) for (t in 1:TT) Psi[, u, t] = pmax(mvtnorm:::dmvnorm(matrix(Y[, t, ], n, r), Mu[, u], Si), 0.1^300)
  
  L[, , 1] = Psi[, , 1] %*% diag(piv)
  if (n == 1) 
    Lt = sum(L[, , 1])
  else Lt = rowSums(L[, , 1])
  lk = sum(log(Lt))
  L[, , 1] = L[, , 1]/Lt
  for (t in 2:TT) {
    L[, , t] = Psi[, , t] * (L[, , t - 1] %*% Pi[, , t])
    if (n == 1) 
      Lt = sum(L[1, , t])
    else Lt = rowSums(L[, , t])
    lk = lk + sum(log(Lt))
    L[, , t] = L[, , t]/Lt
  }
  
  out = list(lk = lk, Psi = Psi, L = L)
}