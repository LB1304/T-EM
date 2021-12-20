## Authors: Luca Brusa, Francesco Bartolucci, Fulvia Pennoni ##

long2matrices_cont <- function (Y, id, time) {
  
  Y = as.matrix(Y)
  ny = ncol(Y)
  TT = max(time)
  idu = unique(id)
  n = length(idu)
  
  YY = array(NA, c(n, TT, ny))
  for (i in 1:n) {
    ind = which(id == idu[i])
    tmp = 0
    for (t in time[ind]) {
      tmp = tmp + 1
      YY[i, t, ] = Y[ind[tmp], ]
    }
  }
  
  return(YY)
}