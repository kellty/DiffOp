# F_coef, U_coef: n*p matrices, (i,j)-element = ith function's jth Fourier coefficient

fit_para <- function(F_coef, U_coef, D_mat, theta_ini, lr=1e-5, tol=1e-8) {
# D_mat: function, returns p*p matrix representation of differential operator D
  # D(phi_1,...,phi_p) = (phi_1,...,phi_p)D_mat
  theta_new <- theta_ini
  for (iter in 1:1e5) {
    theta <- theta_new
    Ftheta_coef <- U_coef %*% t(D_mat(theta))
    loss <- sum((F_coef - Ftheta_coef)^2)
    grad <- sapply(1:length(theta_ini), function(j) {
      theta_new <- theta; theta_new[j] <- theta[j] + 1e-4
      Ftheta_new_coef <- U_coef %*% t(D_mat(theta_new))
      sum((F_coef - Ftheta_new_coef)^2) - loss
    }) * 1e4
    theta_change <- lr * (- grad)
    if (max(abs(theta_change))<tol) {break}
    else { theta_new <- theta + theta_change }
  }
  return(list(Ftheta_coef=Ftheta_coef, theta=theta, iter=iter))
}

fit_RKHS <- function(F_coef, U_coef, K, KL, lam, U_new_coef=0) {
# K, KL: p^2*p^2 kernel matrices
# lam: tuning parameter
  p <- ncol(F_coef)
  U.KL <- kronecker(diag(p), U_coef) %*% KL
  U.KL.2 <- t(KL) %*% kronecker(diag(p), t(U_coef)%*%U_coef) %*% KL
  inv <- solve(U.KL.2 + n*lam*K, t(U.KL))
  reg_coef <- inv %*% c(F_coef)
  Fhat_coef <- matrix(U.KL %*% reg_coef, ncol=p)
  H <- U.KL %*% inv
  trH <- sum(diag(H)) / p
  gcv <- n * sum((F_coef - Fhat_coef)^2) / (n - trH)^2
  if (sum(U_new_coef^2)>0) {
    U_new.KL <- kronecker(diag(p), U_new_coef) %*% KL
    F_new_coef <- matrix(U_new.KL %*% reg_coef, ncol=p)
    return(list(F_new_coef=F_new_coef, Fhat_coef=Fhat_coef, H=H, gcv=gcv))
  }
  return(list(Fhat_coef=Fhat_coef, H=H, gcv=gcv))
}

# H: np*np hat/smoothing matrix applied to Fourier coefficients
# res: n*p matrix of residuals of Fourier coefficients

test.stat <- function(H, res) {
  sum((H %*% c(res))^2) / nrow(res)
}
multiplier <- function(n) {
  ifelse(runif(n)<(1+sqrt(5))/(2*sqrt(5)), (1-sqrt(5))/2, (1+sqrt(5))/2)
}
test.btsp <- function(H, res, B=1e3) {
  replicate(B, test.stat(H, diag(multiplier(n)) %*% res))
}
test.pval <- function(stat, stat_btsp) {
  1 - mean(stat>=stat_btsp)
}
