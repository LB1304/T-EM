## Authors: Luca Brusa, Francesco Bartolucci, Fulvia Pennoni ##

LM_cont_sv <- function(Y, k) {

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
  Yv <- matrix(Y, n * TT, r)

  Mu = matrix(0, r, k)
  mu = colMeans(Yv, na.rm = TRUE)
  Si = cov(Yv, use = "complete.obs")
  for (u in 1:k) Mu[, u] = mvtnorm:::rmvnorm(1, mu, Si)
  Pi = array(stats::runif(k^2 * TT), c(k, k, TT))
  for (t in 2:TT) Pi[, , t] = diag(1/rowSums(Pi[, , t])) %*% Pi[, , t]
  Pi[, , 1] = 0
  piv = stats::runif(k)
  piv = piv/sum(piv)

  return(list(piv = piv, Pi = Pi, Mu = Mu, Si = Si))
}
