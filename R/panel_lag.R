# Create lagged series for panel data frame (Internal helper function)
# @param x data matrix (dependent variables and periods)
# @param k number of lags
panel_lag <- function(x, k) {
  # first column contains time
  # other columns the data
  # k number of lags
  res <-
    rbind(matrix(NA, 1, ncol(x) - 1), as.matrix(x[1:(nrow(x) - 1), -1]), deparse.level = 1)
  colnames(res) <- colnames(x)[-1]

  if (k > 1) {
    for (l in 2:k) {
      res <-
        cbind(res, rbind(matrix(NA, l, ncol(x) - 1), as.matrix(x[1:(nrow(x) - l), -1])))
    }
  }

  #if (k > 1){
  res[rep(is.na(x[, -1]), k)] <- NA # multiple lags
  #}

  # only necessary if single variable with multiple lags
  if(ncol(x) == 2) colnames(res) <- rep(colnames(res)[1], ncol(res))

  colnames(res) <- paste0("lag", rep(1:k, each = ncol(x)-1), "_", colnames(res))
  rownames(res) <- NULL
  res
}


# Panel forward lags
# @param x data matrix (first column contains time aka periods, other columns the data)
# @param k number of lags
#
panel_forward_lag <- function(x, k) {
  res <-
    rbind(as.matrix(x[(1 + 1):nrow(x),-1]), matrix(NA, 1, ncol(x) - 1), deparse.level = 1)
  colnames(res) <- colnames(x)[-1]

  if (k > 1) {
    for (l in 2:k) {
      res <-
        cbind(res, rbind(as.matrix(x[(1 + l):nrow(x),-1]), matrix(NA, l, ncol(x) - 1)))
    }
  }

  #if (k > 1){
  res[rep(is.na(x[,-1]), k)] <- NA # multiple lags
  #}

  # only necessary if single variable with multiple lags
  if (ncol(x) == 2)
    colnames(res) <- rep(colnames(res)[1], ncol(res))

  colnames(res) <-
    paste0("forwardlag", rep(1:k, each = ncol(x) - 1), "_", colnames(res))
  rownames(res) <- NULL
  res
}

# Helmert transformation (Internal helper function)
# @param x data
helmert <- function(x) {
  T <- length(x)
  y <- NULL
  y[1] <- NA
  for (t in 1:(T - 1)) {
    #number of future NAs at time
    nas <- sum(is.na(x[-(1:t)]))

    #Set NAs for transformed observations if all future untransformed obersvations are NAs
    if (length(x[-(1:t)]) == nas) {
      y[t + 1] <- NA
    } else{
      y[t + 1] <-
        sqrt((T - t - nas) / (T - t - nas + 1)) * (x[t] - mean(x[-(1:t)], na.rm =
                                                                 TRUE))
    }
  }
  y
}

# First-difference of forward orthogonal deviations (Internal helper function)
# @param x data
# @param transformation First-difference \code{"fd"} or forward orthogonal deviations \code{"fod"}
panel_fd_fod <- function(x, transformation) {
  if (transformation == "fd") {
    res <- rbind(NA, apply(x[, -(1:2)], 2, diff))
    colnames(res) <- paste("fd", colnames(res), sep = "_")
  }
  if (transformation == "fod") {
    res <- apply(x[, -(1:2)], 2, helmert)
    colnames(res) <- paste("fod", colnames(res), sep = "_")
  }
  rownames(res) <- NULL
  res
}

