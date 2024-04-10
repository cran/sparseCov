#' This function construct a covariance matrix with a block diagonal structure.
#'
#' @param p The number of variants.
#' @param block.size The block size.
#'
#' @return A covariance matrix with a block diagonal structure.
#' @export block.true.cov
#' 
#' @import Matrix
#' @importFrom stats runif
#' @examples
#' data.true.cov <- block.true.cov(30)
#' data.true.cov[1:9,1:9]

block.true.cov <- function(p, block.size = 3){
  block.ind <- as.integer(p/block.size)
  block.list <- list()
  for(b in 1:block.ind){
    
    A <- matrix(runif(block.size^2, 0.1, 0.6), ncol=block.size) 
    diag(A) <- rep(1, block.size)
    A <- forceSymmetric(A)
    
    block.list[[b]] <- A
  }
  as.matrix(bdiag(block.list))
}

