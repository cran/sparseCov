#' This function samples MVN based on a given covariance matrix
#'
#' @param n The sample size.
#' @param Sigma The covariance matrix.
#' @param sparse The indicator of sparse sampling or not.
#' @param n_cores The number of cores used.
#' @param fastmvn The indicator of fast sampling or not.
#'
#' @return The data matrix sampled from the covariance matrix.
#' @export sampleMVN
#' @import Matrix
#' @import sparseMVN
#' @import mvnfast
#' @import methods
#' 
#' @examples
#' ## generate data from a block diagonal covariance matrix structure
#' n <- 50
#' p <- 30
#' data.true.cov <- block.true.cov(p)
#' data <- sampleMVN(n, data.true.cov, sparse=TRUE)
#' data[1:10, 1:10]


sampleMVN <- function(n,
                      Sigma, 
                      sparse=TRUE,
                      n_cores = 1, 
                      fastmvn = FALSE) {
  if(sparse){
    mvnrv <- sampleMVN.sparse(n, Sigma)
  }else{
    if(fastmvn) {
      mvnrv <- mvnfast::rmvn(n, mu = rep(0, dim(Sigma)[1]), sigma = Sigma, ncores = n_cores)
    } else {
      mvnrv <-
        rmvnorm(n, mean = rep(0, dim(Sigma)[1]), sigma = Sigma, checkSymmetry = FALSE, method="eigen")
    }
  }
  
  #mvnrvq <- apply(mvnrv, 2, stats::pnorm)
  
  return(mvnrv)
}


## Sparse version
sampleMVN.sparse <- function(n, Sigma){
  CH <- Cholesky(as(Sigma, 'dsCMatrix'))
  p <- dim(Sigma)[1]
  mvnrv <- rmvn.sparse(n=n, rep(0,p), CH, prec=FALSE) 
  return(mvnrv)
}



## fix integer overflow issue when # of gene* # of cell is too larger in rnorm(n * ncol(sigma))
## This function comes from package Mvnorm.
rmvnorm <- function(n, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)),
                    method=c("eigen", "svd", "chol"), pre0.9_9994 = FALSE, checkSymmetry = TRUE)
{
  if (checkSymmetry && !isSymmetric(sigma, tol = sqrt(.Machine$double.eps),
                                    check.attributes = FALSE)) {
    stop("sigma must be a symmetric matrix")
  }
  if (length(mean) != nrow(sigma))
    stop("mean and sigma have non-conforming size")
  
  method <- match.arg(method)
  
  R <- if(method == "eigen") {
    ev <- eigen(sigma, symmetric = TRUE)
    if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))){
      warning("sigma is numerically not positive semidefinite")
    }
    ## ev$vectors %*% diag(sqrt(ev$values), length(ev$values)) %*% t(ev$vectors)
    ## faster for large  nrow(sigma):
    t(ev$vectors %*% (t(ev$vectors) * sqrt(pmax(ev$values, 0))))
  }
  else if(method == "svd"){
    s. <- svd(sigma)
    if (!all(s.$d >= -sqrt(.Machine$double.eps) * abs(s.$d[1]))){
      warning("sigma is numerically not positive semidefinite")
    }
    t(s.$v %*% (t(s.$u) * sqrt(pmax(s.$d, 0))))
  }
  else if(method == "chol"){
    R <- chol(sigma, pivot = TRUE)
    R[, order(attr(R, "pivot"))]
  }
  
  retval <- matrix(stats::rnorm(as.double(n) * ncol(sigma)), nrow = n, byrow = !pre0.9_9994) %*%  R
  retval <- sweep(retval, 2, mean, "+")
  colnames(retval) <- names(mean)
  retval
}
