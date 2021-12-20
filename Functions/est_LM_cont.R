## Authors: Luca Brusa, Francesco Bartolucci, Fulvia Pennoni ##

est_LM_cont <- function(data, index, k, modBasic = 0, tol = 1e-8, maxit = 1e6, sv = NULL,
                        algorithm = 0, profile_pars = list(alpha = NULL, beta = NULL, ro = NULL, T0 = NULL)) {
  
  id.which <- which(names(data) == index[1])
  tv.which <- which(names(data) == index[2])
  id <- data[, id.which]
  tv <- data[, tv.which]
  Y <- data[, -c(id.which, tv.which), drop = FALSE]
  Y_names <- colnames(Y)
  
  Y <- long2matrices_cont(Y = Y, id = id, time = tv)
  dimnames(Y)[[3]] <- Y_names
  
  # Starting values
  if (is.null(sv)) {
    sv <- LM_cont_sv(Y = Y, k = k)
  }
  
  # EM algorithm
  if (algorithm == 0) {
    out <- LM_cont_em(Y = Y, k = k, modBasic = modBasic, tol = tol, maxit = maxit, piv = sv$piv, Pi = sv$Pi, Mu = sv$Mu, Si = sv$Si)
  } else if (algorithm %in% c(1, 2, 3)) {
    out <- LM_cont_tem(Y = Y, k = k, modBasic = modBasic, tol = tol, maxit = maxit, piv = sv$piv, Pi = sv$Pi, Mu = sv$Mu, Si = sv$Si,
                       algorithm = algorithm, profile_pars = profile_pars)
  } else {
    stop("Specify a correct tempering profile.")
  }
  
  return(out)
}