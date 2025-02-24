ols <- function(Y, X, se = TRUE) {
  
  X <- as.matrix(X)
  Y <- as.vector(Y)
  n <- nrow(X)
  d <- ncol(X)
  sXX <- (1/n)*crossprod(X, X)
  sXY <- (1/n)*crossprod(X, Y)
  b <- solve(sXX, sXY)
  
  if(se){
    Omega <- matrix(0, nrow = d, ncol = d)
    for(i in 1:length(Y)){
      x <- X[i, , drop=FALSE]
      u <- Y[i] - drop(x %*% b)
      Omega <- Omega + u^2 + crossprod(x)
    }
    V <- solve(sXX) %*% Omega %*% solve(sXX) / n^2
    return(list(b=b, V=V, sXX=sXX))
  }
  
  else{
    return(b)
  }
  
}

# OLS estimator with additive bias correction
ols_bca <- function(Y, Xhat, fpr, m){
  d <- ncol(Xhat)
  b_V_sXX <- ols(Y, Xhat)
  b_ <- b_V_sXX$b
  V_ <- b_V_sXX$V
  sXX <- b_V_sXX$sXX
  
  A <- matrix(0, nrow = d, ncol = d)
  A[1,1] <- 1.0
  Γ <- solve(sXX, A)
  b <- b_ + fpr * Γ %*% b_
  I_d = diag(d)
  V = (I_d + fpr * Γ) %*% V_ %*% t(I_d + fpr * Γ) + fpr * (1.0 - fpr) * Γ %*% (V_ + b * t(b)) %*% t(Γ) / m
  return(list(b=b, V=V))
}

# OLS estimator with multiplicative bias correction
ols_bcm <- function(Y, Xhat, fpr, m){
  d = ncol(Xhat)
  b_V_sXX <- ols(Y, Xhat)
  b_ = b_V_sXX$b
  V_ = b_V_sXX$V
  sXX = b_V_sXX$sXX
  
  A = matrix(0, nrow = d, ncol = d)
  A[1,1] <- 1.0
  Γ <- solve(sXX, A)
  I_d = diag(d)
  b = solve(I_d - fpr * Γ) %*% b_
  
  V = solve(I_d - fpr * Γ) %*% V_ %*% t(solve(I_d - fpr * Γ)) + fpr * (1.0 - fpr) * Γ %*% (V_ + b %*% t(b)) %*% t(Γ) / m
  return(list(b=b, V=V))
}