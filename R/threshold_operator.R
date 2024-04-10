#' This function computes the thresholding sparse covariance estimator for a given threshold level.
#'
#' @param z The sample covariance matrix.
#' @param operator The choice of the thresholding operator.
#' @param delta The thresholding level.
#' @param n The sample size of data matrix.
#'
#' @return The thresholding sparse covariance estimator for a given threshold level.
#' @export thresh_op
#'
#' @examples
#' ## generate data from a block diagonal covariance matrix structure
#' n <- 50
#' p <- 30
#' data.true.cov <- block.true.cov(p)
#' data <- sampleMVN(n, data.true.cov, sparse=TRUE)
#' ## compute the sample covariance
#' z <- Rfast::cova(data) *(n-1)/n
#' ## get the sparse covariance matrix estimator for a given threshold level
#' s <- thresh_op(z, operator='soft', delta=1, n=n)
#' s[1:9,1:9]

thresh_op <- function(z, operator, delta, n){
  if(operator == 'hard'){
    s_method <- s_hard
  }else if(operator == 'soft'){
    s_method <- s_soft
  }else if(operator == 'scad'){
    s_method <- s_scad
  }else if(operator == 'al'){
    s_method <- s_al
  }else{
    stop('Please specify a valid thresholding operator.')
  }
  s_method(z, delta, n)
}

# Operator 1: Hard Thresholding
s_hard <- function(z, delta, n){
  
  p<-dim(z)[2]
  lambda <- sqrt(log(p)/n)*delta
  output <- (z>lambda)*z
  diag(output) <- diag(z)
  return(output)
}

# Operator 2: Soft Thresholding
s_soft <- function(z, delta, n){
  p <- dim(z)[1]
  lambda <- sqrt(log(p)/n)*delta
  z0 <- abs(z)-lambda
  output <- sign(z)*(z0>0)*z0
  diag(output) <- diag(z)
  return(output)
}


# Operator 3: SCAD (smoothly clipped absolute deviation)
s_scad <- function(z, delta, n, a=3.7){
  p <- dim(z)[1]
  lambda <- sqrt(log(p)/n)*delta
  
  output <- matrix(NA, dim(z)[1], dim(z)[2])
  
  index1 <- which(abs(z)<=2*lambda)
  z0 <- abs(z)-lambda
  output[index1] <- sign(z[index1])*(z0[index1]>0)*z0[index1]
  
  index2 <- which(abs(z)>2*lambda & abs(z)<=a*lambda)
  output[index2] <- ((a-1)*z[index2]-sign(z[index2])*a*lambda)/(a-2)
  
  index3 <- which(abs(z)>a*lambda)
  output[index3] <- z[index3]
  
  diag(output) <- diag(z)
  return(output)
}


# Operator 4: Adaptive lasso
s_al <- function(z, delta, n){
  p <- dim(z)[1]
  lambda <- sqrt(log(p)/n)*delta
  
  eta <- 3
  z0 <- abs(z) - lambda^(eta+1)*abs(z)^(-eta)
  output <- sign(z)*(z0>0)*z0
  diag(output) <- diag(z)
  return(output)
}



# # variation of sample variance
# theta_est <- function(data){
#   n <- dim(data)[1]
#   p <- dim(data)[2]
#   
#   theta <- matrix(0, p, p)
#   for(i in 1:n){
#     v <- data[i,] - colMeans(data)
#     theta <- (Matrix::tcrossprod(v,v) - s)^2 + theta
#   }
#   theta <- theta/n
#   theta
# }

