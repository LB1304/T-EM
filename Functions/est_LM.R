## Authors: Luca Brusa, Francesco Bartolucci, Fulvia Pennoni ##

est_LM <- function (data, index, k, modBasic = 0, tol = 1e-8, maxit = 1e6, sv = NULL,
                    algorithm = 0, profile_pars = list(alpha = NULL, beta = NULL, ro = NULL, T0 = NULL)) {

  id.which <- which(names(data) == index[1])
  tv.which <- which(names(data) == index[2])
  id <- data[, id.which]
  tv <- data[, tv.which]
  Y <- data[, -c(id.which, tv.which), drop = FALSE]
  Y_names <- colnames(Y)

  tmp <- long2matrices(Y = Y, id = id, time = tv)
  Y <- tmp$Y
  dimnames(Y)[[3]] <- Y_names
  freq = tmp$freq

  # Starting Values
  if (is.null(sv)) {
    sv <- LM_sv(S = Y, k = k)
  }

  # EM algorithm
  if (algorithm == 0) {
    out <- LM_em(S = Y, yv = freq, k = k, modBasic = modBasic, tol = tol, maxit = maxit, piv = sv$piv, Pi = sv$Pi, Phi = sv$Phi)
  } else if (algorithm %in% c(1, 2, 3)) {
    out <- LM_tem(S = Y, yv = freq, k = k, modBasic = modBasic, tol = tol, maxit = maxit, piv = sv$piv, Pi = sv$Pi, Phi = sv$Phi,
                  algorithm = algorithm, profile_pars = profile_pars)
  } else {
    stop("Specify a correct tempering profile.")
  }

  return(out)
}
