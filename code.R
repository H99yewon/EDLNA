NMF <- function(A, k, max_iter = 500, eps = 1e-9) {
  N <- nrow(A); M <- ncol(A)
  W <- matrix(runif(N * k), N, k)
  H <- matrix(runif(k * M), k, M)
  
  best_W <- W; best_H <- H
  best_rss <- sum((A - W %*% H)^2)
  
  for (iter in 1:max_iter) {
    # H update
    H <- H * ((t(W) %*% A) / (t(W) %*% W %*% H + eps))
    
    # W update
    W <- W * ((A %*% t(H)) / (W %*% H %*% t(H) + eps))
    
    rss <- sum((A - W %*% H)^2)
    if (rss < best_rss) {
      best_rss <- rss
      best_W <- W
      best_H <- H
    }
  }
  list(W = best_W, H = best_H, rss = best_rss)
}

NMF_fixed_H <- function(A, H, k, max_iter = 500, eps = 1e-9) {
  N <- nrow(A); k <- nrow(H)
  W <- matrix(runif(N * k), N, k)
  
  best_W <- W
  best_rss <- sum((A - W %*% H)^2)
  
  for (iter in 1:max_iter) {
    # W update (H is fixed)
    W <- W * ((A %*% t(H)) / (W %*% H %*% t(H) + eps))
    
    rss <- sum((A - W %*% H)^2)
    if (rss < best_rss) {
      best_rss <- rss
      best_W <- W
    }
  }
  list(W = best_W, H = H, rss = best_rss)
}


# exnormal, extumor: m × n matrix
# k: latent factor 

# normal tissue NMF
nmf_normal <- NMF(exnormal, k = k)
H_normal <- nmf_normal$H  


# tumor tissue NMF (fixed H matrix)
nmf_tumor <- NMF_fixed_H(extumor, H_normal)


# correlation matrix
normal_cor <- cor(t(W_normal), method = "pearson")
tumor_cor  <- cor(t(W_tumor),  method = "pearson")
cor_diff <- abs(normal_cor - tumor_cor)
