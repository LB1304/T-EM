## Authors: Luca Brusa, Francesco Bartolucci, Fulvia Pennoni ##

LM_cont_em_step <- function(Y, n, r, TT, k, L, Psi, piv, Pi, Mu, Si, modBasic) {
  
  V = array(0, c(n, k, TT))
  U = array(0, c(k, k, TT))
  M = matrix(1, n, k)
  if (n == 1) 
    V[, , TT] = L[, , TT]/sum(L[1, , TT])
  else V[, , TT] = L[, , TT]/rowSums(L[, , TT])
  
  for (i in 1:n) {
    Tmp = (L[i, , TT - 1] %o% Psi[i, , TT]) * Pi[, , TT]
    U[, , TT] = U[, , TT] + Tmp/sum(Tmp)
  }
  
  if (TT > 2) {
    for (t in seq(TT - 1, 2, -1)) {
      M = (Psi[, , t + 1] * M) %*% t(Pi[, , t + 1])
      M = M/rowSums(M)
      V[, , t] = L[, , t] * M
      if (n == 1) 
        V[, , t] = V[, , t]/sum(V[1, , t])
      else V[, , t] = V[, , t]/rowSums(V[, , t])
      
      for (i in 1:n) {
        Tmp = (L[i, , t - 1] %o% (Psi[i, , t] * M[i, ])) * Pi[, , t]
        U[, , t] = U[, , t] + Tmp/sum(Tmp)
      }
    }
  }
  M = (Psi[, , 2] * M) %*% t(Pi[, , 2])
  M = M/rowSums(M)
  V[, , 1] = L[, , 1] * M
  if (n == 1) 
    V[, , 1] = V[, , 1]/sum(V[1, , 1])
  else V[, , 1] = V[, , 1]/rowSums(V[, , 1])
  Vv = matrix(aperm(V, c(1, 3, 2)), n * TT, k)
  
  Y1 = array(Y, c(n, TT, r, k))
  Var = array(0, c(n, TT, r, r))
  
  Sitmp = matrix(0, r, r)
  for (u in 1:k) {
    Yv1 = matrix(Y1[, , , u], n * TT)
    Var1 = array(Var, c(n * TT, r, r))
    Mu[, u] = (t(Yv1) %*% Vv[, u])/sum(Vv[, u])
    Tmp = Yv1 - rep(1, n * TT) %*% t(Mu[, u])
    Sitmp = Sitmp + t(Tmp) %*% (Vv[, u] * Tmp) + apply(Vv[, u] * Var1, c(2, 3), sum)
  }
  Si = Sitmp/(n * TT)
  
  piv = colSums(V[, , 1])/n
  U = pmax(U, 10^-300)
  if (modBasic == 0) 
    for (t in 2:TT) Pi[, , t] = diag(1/rowSums(U[, , t])) %*% U[, , t]
  if (modBasic == 1) {
    Ut = apply(U[, , 2:TT], c(1, 2), sum)
    Pi[, , 2:TT] = array(diag(1/rowSums(Ut)) %*% Ut, c(k, k, TT - 1))
  }
  if (modBasic > 1) {
    Ut1 = U[, , 2:modBasic]
    if (length(dim(Ut1)) > 2) 
      Ut1 = apply(Ut1, c(1, 2), sum)
    Ut2 = U[, , (modBasic + 1):TT]
    if (length(dim(Ut2)) > 2) 
      Ut2 = apply(Ut2, c(1, 2), sum)
    Pi[, , 2:modBasic] = array(diag(1/rowSums(Ut1, 2)) %*% Ut1, c(k, k, modBasic - 1))
    Pi[, , (modBasic + 1):TT] = array(diag(1/rowSums(Ut2, 2)) %*% Ut2, c(k, k, TT - modBasic))
  }
  
  return(list(piv = piv, Pi = Pi, Mu = Mu, Si = Si, V = V))
}