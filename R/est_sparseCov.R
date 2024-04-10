#' This function computes the thresholding sparse covariance/correlation estimator 
#' with the optimal threshold level.
#'
#' @param data The data matrix.
#' @param method The choice of method to select the optimal threshold level.
#' @param operator The choice of the thresholding operator.
#' @param corr The indicator of computing correlation or covariance matrix.
#'
#' @return The thresholding sparse covariance/correlation estimator.
#' @export est_sparseCov
#' @importFrom Rfast cova
#' @import Matrix
#' 
#' 
#' @examples
#' ## generate data from a block diagonal covariance matrix structure
#' n <- 50
#' p <- 30
#' data.true.cov <- block.true.cov(p)
#' data <- sampleMVN(n, data.true.cov, sparse=TRUE)
#' ## compute the thresholding sparse covariance/correlation estimator
#' s <- est_sparseCov(data, method='cv', operator='scad', corr=FALSE)

est_sparseCov <- function(data, 
                         method=c('cv', 'qiu'),
                         operator=c('hard', 'soft', 'scad', 'al'),
                         corr=TRUE){
  p <- dim(data)[1]
  n <- dim(data)[2]
  
  # sample covariance
  z <- Rfast::cova(data) *(n-1)/n
  
  # select the optimal thresholding level
  delta <- est_delta(data, method=method, operator=operator)
  s <- thresh_op(z, operator=operator, delta=delta, n=n)
  
  # Modify s to make it psd
  tol <- 1e-6
  ev <- eigen(s, symmetric=TRUE, only.values = TRUE)$values
  s1 <- s + (tol-min(ev))*diag(dim(s)[1]) 
  
  if(corr){
    # make corr
    s1_corr <- Matrix::cov2cor(s1)
    output <- s1_corr
  }else{
    output <- s1
  }
  
  return(output)
}