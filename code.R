# Standard Non-negative Matrix Factorization (NMF) with multiplicative updates

NMF <- function(A, rank, max_iter = 1000, tol = 1e-5, eps = 1e-9) {
  N <- nrow(A)
  M <- ncol(A)
  
  # Initialize W and H with non-negative random values
  W <- matrix(runif(N * rank, min=0, max=1), N, rank)
  H <- matrix(runif(rank * M, min=0, max=1), rank, M)
  
  # Compute initial reconstruction error (Residual Sum of Squares, RSS)
  prev_rss <- sum((A - W %*% H)^2) / 2
  
  for (iter in 1:max_iter) {
    # Update rule for H (multiplicative update, avoiding division by zero with eps)
    numerator_H <- t(W) %*% A
    denominator_H <- t(W) %*% W %*% H + eps
    H <- H * (numerator_H / denominator_H)
    
    # Update rule for W (multiplicative update)
    numerator_W <- A %*% t(H)
    denominator_W <- W %*% H %*% t(H) + eps
    W <- W * (numerator_W / denominator_W)
    
    # Compute new reconstruction error (RSS)
    rss <- sum((A - W %*% H)^2) / 2
    
    # Convergence criterion: relative change in RSS below tolerance
    if ((abs(rss - prev_rss) / prev_rss) < tol) {
      cat(sprintf("NMF Converged at iteration %d with RSS = %.6f\n", iter, rss))
      break
    }
    
    prev_rss <- rss
  }
  
  # Return factor matrices and final reconstruction error
  return(list(W = W, H = H, rss = rss))
}


# NMF with fixed H: only W is updated
NMF_fixed_H <- function(A, H, max_iter = 1000, tol = 1e-5, eps = 1e-9) {
  N <- nrow(A)
  rank <- nrow(H)
  
  # Initialize W with non-negative random values
  W <- matrix(runif(N * rank, min=0, max=1), N, rank)
  
  # Compute initial reconstruction error (RSS)
  prev_rss <- sum((A - W %*% H)^2) / 2
  
  for (iter in 1:max_iter) {
    # Update rule for W while keeping H fixed
    numerator_W <- A %*% t(H)
    denominator_W <- W %*% H %*% t(H) + eps
    W <- W * (numerator_W / denominator_W)
    
    # Compute new RSS
    rss <- sum((A - W %*% H)^2) / 2
    
    # Convergence criterion: relative change in RSS below tolerance
    if ((abs(rss - prev_rss) / prev_rss) < tol) {
      cat(sprintf("NMF_fixed_H Converged at iteration %d with RSS = %.6f\n", iter, rss))
      break
    }
    prev_rss <- rss
  }
  
  # Return factor matrices and final reconstruction error
  return(list(W = W, H = H, rss = rss))
}
