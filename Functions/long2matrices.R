## Authors: Luca Brusa, Francesco Bartolucci, Fulvia Pennoni ##

long2matrices <- function (Y, id, time) {

  idu = unique(id)
  n = length(idu)
  TT = max(time)

  Y = as.matrix(Y)
  ny = ncol(Y)

  temp <- data.frame(id = id, time = time, Y, check.names = FALSE)
  temp.wide <- reshape(temp, idvar = "id", timevar = "time", direction = "wide")
  temp.wide[is.na(temp.wide)] <- 999
  aggr <- MultiLCIRT::aggr_data(temp.wide[, -1], fort = TRUE)
  temp <- data.frame(1:nrow(aggr$data_dis), aggr$data_dis)
  colnames(temp) <- c("id", attributes(temp.wide)$reshapeWide$varying)
  freq <- aggr$freq
  Y <- stats::reshape(temp, direction = "long", idvar = "id", varying = 2:ncol(temp), sep = ".")
  Y[, -c(1, 2)][Y[, -c(1, 2)] == 999] <- NA
  id <- Y[, 1]
  time <- Y[, 2]
  Y = as.matrix(Y[, -c(1, 2)])

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

  return(list(Y = YY, freq = freq))
}




