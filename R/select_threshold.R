#' This function select the optimal thresholding level delta
#'
#' @param data The data matrix.
#' @param method The choice of method to select the optimal threshold level.
#' @param operator The choice of thresholding operator.
#'
#' @return The optimal threshold level.
#' @export est_delta
#' @importFrom Rfast cova
#' @importFrom stats pnorm optimize
#' @examples
#' ## generate data from a block diagonal covariance matrix structure
#' n <- 50
#' p <- 30
#' data.true.cov <- block.true.cov(p)
#' data <- sampleMVN(n, data.true.cov, sparse=TRUE)
#' ## select the optimal thresholding level delta
#' delta <- est_delta(data, method='cv', operator='scad')

est_delta <- function(data, 
                     method=c('cv', 'qiu'),
                     operator=c('hard', 'soft', 'scad', 'al')){
  n <- dim(data)[2]
  if((method=='qiu') ){
    s <- cova(data) *(n-1)/n
    delta <- qiu.select(data, s)
  }else if(method=='cv'){
    delta <- cv.min(data, operator)
  }else{
    stop('Please specify a valid thresholding method and an operator function.')
  }
  return(delta)
}


### Qiu function to tune delta
qiu.select = function(data, s=NULL){
  n <- dim(data)[1]
  p <- dim(data)[2]
  
  if(is.null(s)){
    s <- Rfast::cova(data) *(n-1)/n
  }
  
  # standardized covariance of Sigma
  eta <- sqrt(n/log(p))*s 
  # select lower triangular
  ltrig <- abs(eta[lower.tri(eta, diag = FALSE)])
  
  # parameter to select the optimal thr
  a <- min(sqrt(2+log(n)/log(p)), 2)
  a1 <- 2 - a
  a0 <- (log(log(p)))^(-1/2)
  M <- sum(((a1+a0)<ltrig) & (ltrig<=2))
  q <- (p-1)*p/2
  plog <- (log(p))^(1/2)
  V <- 2*q*(pnorm(2*plog)-pnorm((a1+a0)*plog))
  N2 <- max(M-V, plog)
  
  delta <- sqrt(2*(2-log(N2*1/plog)/log(p)))
  
  return(delta)
}

#### CV to tune delta
# The tuning parameter lambda was selected by minimizing the Frobenius norm of the 
# difference between s(sample_cov(train)) and sample_cov(test) by cross validation

# cv.select: For implementation, we choose the optimal delta from a delta.list 
# instead of really solving the optimization problem.
# lambda.list ranges in [0,lambda.max], with a user specified length.

# cv.min: solve the optimization problem directly (currently use this)

cv.select <- function(data, 
                      operator, 
                      fold=5, 
                      delta.length=10, 
                      delta.max=2){
  
  n <- dim(data)[1]
  p <- dim(data)[2]
  
  
  # Create delta list
  delta.list <- seq(0, delta.max, length.out=delta.length)
  
  # Create equal size folds
  folds <- cut(seq(1, nrow(data)), breaks=fold, labels=FALSE)
  
  fold_losses <- rep(0, fold) # Record fold loss 
  cv_errors <- rep(0, delta.length) # Record average estimated loss
  
  
  ## Perform k fold cross validation
  
  # iterate over the values of lambda
  for(j in 1:delta.length){
    
    # iterate over CV folds
    for(i in 1:fold){
      # Segment the data by fold 
      testIndexes <- which(folds==i, arr.ind=TRUE) # i-th fold as the testing data
      #trainIndexes <- which(folds!=i, arr.ind=TRUE) # i-th fold as the testing data
      testData <- data[testIndexes, ]
      trainData <- data[-testIndexes, ]
      # Get the covariance matrix estimator s based on the training data
      sample.cov.train <- cova(trainData, center=TRUE, large = TRUE)
      s <- thresh_op(sample.cov.train, operator=operator, delta = delta.list[j], n=n)
      # Compute Frobenius risk = norm(distance of s and covariance of testing data)
      sample.cov.test <- cova(testData, center=TRUE, large = TRUE)
      fold_losses[i] <- norm(s - sample.cov.test, type='F')
    }
    # Get the average estimated loss of CV
    cv_errors[j] <-  mean(fold_losses)
  }
  
  # Find the lambda that minimize cv_errors
  delta <- delta.list[which.min(cv_errors)]
  
  delta
}




cv.min <- function(data, 
                   operator, 
                   fold=5,
                   delta.max=2){
  n <- dim(data)[1]
  p = dim(data)[2]
  
  
  ## Perform k fold cross validation
  cv.loss <- function(delta){
    # Create equal size folds
    folds <- cut(seq(1, nrow(data)), breaks=fold, labels=FALSE)
    fold_losses <- rep(0, fold) # Record fold loss 
    
    # iterate over CV folds
    for(i in 1:fold){
      # Segment the data by fold 
      testIndexes <- which(folds==i, arr.ind=TRUE) # i-th fold as the testing data
      testData <- data[testIndexes, ]
      trainData <- data[-testIndexes, ]
      # Get the covariance matrix estimator s based on the training data
      sample.cov.train <- cova(trainData, center=TRUE, large = TRUE)
      s <- thresh_op(sample.cov.train, operator=operator, delta = delta, n=n)
      # Compute Frobenius risk = norm(distance of s and covariance of testing data)
      sample.cov.test <- cova(testData, center=TRUE, large = TRUE)
      fold_losses[i] <- norm(s - sample.cov.test, type='F')
    }
    # Get the average estimated loss of CV
    cv_errors <-  mean(fold_losses)
    cv_errors
  }
  
  # Solve the optimization problem
  res <- optimize(cv.loss, c(0, delta.max))
  delta <- res$minimum
  
  delta
}



