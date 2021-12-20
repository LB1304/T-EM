## Authors: Luca Brusa, Francesco Bartolucci, Fulvia Pennoni ##

LM_tem_step <- function(Sv, yv, n, ns, r, TT, k, bv, b, piv, Pi, Phi, Psi, pv, L, mod, temp) {

  V = array(0, c(ns, k, TT))
  U = array(0, c(k, k, TT))
  Yvp = matrix(1/pv, ns, k)
  M = matrix(1, ns, k)
  Tmp = L[, , TT]^(1/temp)
  V[, , TT] = (yv/rowSums(Tmp))*Tmp
  Tmp = (1/pv)*(Psi[, , TT]*M)
  for (i in 1:ns) {
    aux <- ((L[i, , TT - 1] %o% Tmp[i,]) * Pi[, , TT])^(1/temp)
    U[, , TT] = U[, , TT] + yv[i] * aux/sum(aux)
  }
  if (TT > 2) {
    for (t in seq(TT - 1, 2, -1)) {
      M = (Psi[, , t + 1] * M) %*% t(Pi[, , t + 1])
      Tmp = (L[, , t] * M)^(1/temp)
      V[, , t] = (yv/rowSums(Tmp))*Tmp
      Tmp = (1/pv)*Psi[, , t]*M
      for (i in 1:ns) {
        aux <- ((L[i, , t - 1] %o%Tmp[i,]) * Pi[, , t])^(1/temp)
        U[, , t] = U[, , t] + yv[i] * aux/sum(aux)
      }
    }
  }
  M = (Psi[, , 2] * M) %*% t(Pi[, , 2])
  Tmp = (L[, , 1] * M)^(1/temp)
  V[, , 1] = (yv/rowSums(Tmp))*Tmp


  Y1 = array(NA, c(b + 1, k, r))
  for (j in 1:r) Y1[1:(bv[j] + 1)] = 0
  Vv = matrix(aperm(V, c(1, 3, 2)), ns * TT, k)
  for (j in 1:r) for (jb in 0:bv[j]) {
    ind = which(Sv[, j] == jb)
    if (length(ind) == 1) {
      Y1[jb + 1, , j] = Vv[ind, ]
    }
    if (length(ind) > 1) {
      Y1[jb + 1, , j] = colSums(Vv[ind, ])
    }
  }
  for (j in 1:r) for (c in 1:k) {
    tmp = Y1[1:(bv[j] + 1), c, j]
    if (any(is.na(tmp)))
      tmp[is.na(tmp)] = 0
    tmp = pmax(tmp/sum(tmp), 10^-10)
    Phi[1:(bv[j] + 1), c, j] = tmp/sum(tmp)
  }

  # Definisco piv
  piv = colSums(V[, , 1])/n

  # Definisco Pi
  U = pmax(U, 10^-300)
  if (mod == 0)
    for (t in 2:TT) Pi[, , t] = diag(1/rowSums(U[, , t])) %*% U[, , t]
  if (mod == 1) {
    Ut = apply(U[, , 2:TT], c(1, 2), sum)
    Pi[, , 2:TT] = array(diag(1/rowSums(Ut)) %*% Ut, c(k, k, TT - 1))
  }
  if (mod > 1) {
    Ut1 = U[, , 2:mod]
    if (length(dim(Ut1)) > 2)
      Ut1 = apply(Ut1, c(1, 2), sum)
    Ut2 = U[, , (mod + 1):TT]
    if (length(dim(Ut2)) > 2)
      Ut2 = apply(Ut2, c(1, 2), sum)
    Pi[, , 2:mod] = array(diag(1/rowSums(Ut1, 2)) %*% Ut1, c(k, k, mod - 1))
    Pi[, , (mod + 1):TT] = array(diag(1/rowSums(Ut2, 2)) %*% Ut2, c(k, k, TT - mod))
  }

  return(list(piv = piv, Pi = Pi, Phi = Phi, V = V))
}
