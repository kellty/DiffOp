source('func.R')
MC <- 200
library(foreach)
library(doParallel)

phi <- function(x, k) { sqrt(2) * cos(k*pi*x) }
Lphi <- function(x, k) { sqrt(2) * (k*pi)^2 * cos(k*pi*x) }
p <- 25;  p2 <- p^2

# Kx_list <- Ky_list <- LKy_list <- list()
# bw_list <- c(1e-2,2e-2,5e-2)
# for (kernel_type in 0:1) for (bw_idx in 1:length(bw_list)) {
#   idx_K <- kernel_type*length(bw_list) + bw_idx
#   bw <- bw_list[bw_idx]
#   if (kernel_type==0) {
#     K1 <- function(x,y) { dnorm(x-y, sd=bw) }
#     LK1 <- function(x,y) { ((x-y)^2/bw^4 - 1/bw^2) * dnorm(x-y, sd=bw) }
#   } else {
#     K1 <- function(x,y) { (1 + (x-y)^2/(4*bw^2))^(-2) }
#     LK1 <- function(x,y) { tmp <- (x-y)^2/(4*bw^2)
#       - bw^(-2) * (1 + tmp)^(-4) * (1 - 5*tmp) }
#   }
#   Kx_list[[idx_K]] <- sapply(1:p, function(k) sapply(1:p, function(k_) 
#     integrate(function(x1vec) sapply(x1vec, function(x1)
#       integrate(function(x2) K1(x1,x2) * Lphi(x1,k_) * Lphi(x2,k), 0,1)$value), 0,1)$value))
#   Ky_list[[idx_K]] <- sapply(1:p, function(j) sapply(1:p, function(j_)
#     integrate(function(y1vec) sapply(y1vec, function(y1)
#       integrate(function(y2) K1(y1,y2) * phi(y1,j_) * phi(y2,j), 0,1)$value), 0,1)$value))
#   LKy_list[[idx_K]] <- sapply(1:p, function(j) sapply(1:p, function(j_)
#     integrate(function(y1vec) sapply(y1vec, function(y1)
#       integrate(function(y2) LK1(y1,y2) * phi(y1,j_) * phi(y2,j), 0,1)$value), 0,1)$value))
# }
# save(Kx_list, Ky_list, LKy_list, file='ker_mat.RData')
load('ker_mat.RData')

n <- 200;  SNR <- 3;  p <- 10
lglams_diff <- seq(-3,6,1);  lglams_int <- seq(-14,-5,1)

U_sd_list <- list(sapply(1:p, function(k) 1/k^3), rep(1e-2, p), 
                  sapply(1:p, function(k) ifelse(k%%2, 1/k^2, 1/k^3)))

cl <- makeCluster(100) # detectCores()
registerDoParallel(cl)
for (align in 1:length(U_sd_list)) for (idx_K in 1:length(Kx_list)) {
Kx <- Kx_list[[idx_K]]; Ky <- Ky_list[[idx_K]]; LKy <- LKy_list[[idx_K]]
idx <- cbind(rep(1:p,each=p), rep(1:p,times=p));  p2 <- p^2
K <- KL <- K_int <- diag(p2)
for (i in 1:p2) for (i_ in 1:p2) {
  j <- idx[i,1];  k <- idx[i,2];  j_ <- idx[i_,1];  k_ <- idx[i_,2]
  K_int[i,i_] <- Ky[k,k_] * Ky[j,j_]
  K[i,i_] <- (Kx[k,k_] + phi(0,k_) * phi(0,k) + phi(1,k_) * phi(1,k)) * Ky[j,j_]
  KL[i,i_] <- (Kx[k,k_] + phi(0,k_) * phi(0,k) + phi(1,k_) * phi(1,k)) * LKy[j,j_]
}
  
U_sd <- U_sd_list[[align]][1:p]
D_vec <- ((1:p)*pi)^2

  idx <- paste0('_kernel',idx_K,'_align',align)
  print(paste(Sys.time(), idx))
result <- foreach(mc=1:MC, .combine=rbind) %dopar% {
  t0 <- proc.time()[3]
  set.seed(999+(mc+55)^2)
U_coef <- sapply(1:p, function(k) runif(n,-1,1)*sqrt(3)) %*% diag(U_sd)
F_sd <- U_sd * D_vec
F0_coef <- U_coef %*% diag(D_vec)
eps_sd <- rep(sqrt(mean(F_sd^2)) / SNR, p)
eps_coef <- sapply(1:p, function(k) rnorm(n, sd=eps_sd[k]))
F_coef <- F0_coef + eps_coef

U_new_coef <- sapply(1:p, function(k) runif(n,-1,1)*sqrt(3)) %*% diag(U_sd)
F0_new_coef <- U_new_coef %*% diag(D_vec)

tmp <- c(sum(F0_coef^2), sum(eps_coef^2), sum(F0_new_coef^2))

for (lam in 10^lglams_diff) {
  fit <- fit_RKHS(F_coef, U_coef, K, KL, lam, U_new_coef)
  tmp <- c(tmp, sum((F0_new_coef - fit$F_new_coef)^2), fit$gcv)
}
for (lam in 10^lglams_int) {
  fit <- fit_RKHS(F_coef, U_coef, K_int, K_int, lam, U_new_coef)
  tmp <- c(tmp, sum((F0_new_coef - fit$F_new_coef)^2), fit$gcv)
}

tmp <- c(proc.time()[3]-t0, tmp)
names(tmp) <- c('time', 'TSS', 'ESS', 'TSS_new', 
    paste0(c('ESS_diff', 'GCV_diff'), '|lgλ=', rep(lglams_diff, each=2)), 
    paste0(c('ESS_int', 'GCV_int'), '|lgλ=', rep(lglams_int, each=2)))
return(tmp)
}
save(result, file=paste0('result',idx,'.RData'))
print(round(apply(result[,2:4], 2, function(r) c(mean(r), sqrt(var(r)/MC))),5))
tmp <- matrix(nrow=10, ncol=length(lglams_diff))
rownames(tmp) <- c('lgλ_diff', 'ESS_diff', 'ESS_diff-sd', 'GCV_diff', 'GCV_diff-sd', 
                   'lgλ_int', 'ESS_int', 'ESS_int-sd', 'GCV_int', 'GCV_int-sd')
tmp['lgλ_diff',] <- lglams_diff
tmp['lgλ_int',] <- lglams_int
for (v in c('ESS_diff', 'GCV_diff')) {
  tmp[c(v, paste0(v,'-sd')),] <- apply(result[,paste0(v,'|lgλ=',lglams_diff)], 2, 
                                       function(r) c(mean(r), sqrt(var(r)/MC)))
}
for (v in c('ESS_int', 'GCV_int')) {
  tmp[c(v, paste0(v,'-sd')),] <- apply(result[,paste0(v,'|lgλ=',lglams_int)], 2, 
                                       function(r) c(mean(r), sqrt(var(r)/MC)))
}
print(round(tmp,3))
}


n <- 200;  SNR <- 3;  p <- 10;  p2 <- p^2
lglams_diff <- seq(-3,6,1);  lglams_int <- seq(-14,-5,1)
Kx <- Kx_list[[1]]; Ky <- Ky_list[[1]]; LKy <- LKy_list[[1]]
idx <- cbind(rep(1:p,each=p), rep(1:p,times=p))
K <- KL <- K_int <- diag(p2)
for (i in 1:p2) for (i_ in 1:p2) {
  j <- idx[i,1];  k <- idx[i,2];  j_ <- idx[i_,1];  k_ <- idx[i_,2]
  K_int[i,i_] <- Ky[k,k_] * Ky[j,j_]
  K[i,i_] <- (Kx[k,k_] + phi(0,k_) * phi(0,k) + phi(1,k_) * phi(1,k)) * Ky[j,j_]
  KL[i,i_] <- (Kx[k,k_] + phi(0,k_) * phi(0,k) + phi(1,k_) * phi(1,k)) * LKy[j,j_]
}
U_sd <- U_sd_list[[1]][1:p]
D_vec <- ((1:p)*pi)^2
for (m in c(15,30,50)) {
  idx <- paste0('_disc_m',m)
  print(paste(Sys.time(), idx))
phi_val <- sapply(((1:m)-0.5)/m, function(x) sapply(1:p, function(k) phi(x,k)))
phi_reg <- solve(phi_val%*%t(phi_val), phi_val)
result <- foreach(mc=1:MC, .combine=rbind) %dopar% {
  t0 <- proc.time()[3]
  set.seed(999+(mc+55)^2)
U_coef <- sapply(1:p, function(k) runif(n,-1,1)*sqrt(3)) %*% diag(U_sd)
F_sd <- U_sd * D_vec
F0_coef <- U_coef %*% diag(D_vec)
eps_sd <- rep(sqrt(mean(F_sd^2)) / SNR, p)
eps_coef <- sapply(1:p, function(k) rnorm(n, sd=eps_sd[k]))
F_coef <- F0_coef + eps_coef

U_new_coef <- sapply(1:p, function(k) runif(n,-1,1)*sqrt(3)) %*% diag(U_sd)
F0_new_coef <- U_new_coef %*% diag(D_vec)

tmp <- c(sum(F0_coef^2), sum(eps_coef^2), sum(F0_new_coef^2))

U_val <- U_coef %*% phi_val + matrix(rnorm(n*m, sd=1e-3), nrow=n)
U_coef <- U_val %*% t(phi_reg)
F_val <- F_coef %*% phi_val + matrix(rnorm(n*m, sd=1e-3), nrow=n)
F_coef <- F_val %*% t(phi_reg)

for (lam in 10^lglams_diff) {
  fit <- fit_RKHS(F_coef, U_coef, K, KL, lam, U_new_coef)
  tmp <- c(tmp, sum((F0_new_coef - fit$F_new_coef)^2), fit$gcv)
}
for (lam in 10^lglams_int) {
  fit <- fit_RKHS(F_coef, U_coef, K_int, K_int, lam, U_new_coef)
  tmp <- c(tmp, sum((F0_new_coef - fit$F_new_coef)^2), fit$gcv)
}

tmp <- c(proc.time()[3]-t0, tmp)
names(tmp) <- c('time', 'TSS', 'ESS', 'TSS_new', 
    paste0(c('ESS_diff', 'GCV_diff'), '|lgλ=', rep(lglams_diff, each=2)), 
    paste0(c('ESS_int', 'GCV_int'), '|lgλ=', rep(lglams_int, each=2)))
return(tmp)
}
save(result, file=paste0('result',idx,'.RData'))
print(round(apply(result[,2:4], 2, function(r) c(mean(r), sqrt(var(r)/MC))),5))
tmp <- matrix(nrow=10, ncol=length(lglams_diff))
rownames(tmp) <- c('lgλ_diff', 'ESS_diff', 'ESS_diff-sd', 'GCV_diff', 'GCV_diff-sd', 
                   'lgλ_int', 'ESS_int', 'ESS_int-sd', 'GCV_int', 'GCV_int-sd')
tmp['lgλ_diff',] <- lglams_diff
tmp['lgλ_int',] <- lglams_int
for (v in c('ESS_diff', 'GCV_diff')) {
  tmp[c(v, paste0(v,'-sd')),] <- apply(result[,paste0(v,'|lgλ=',lglams_diff)], 2, 
                                       function(r) c(mean(r), sqrt(var(r)/MC)))
}
for (v in c('ESS_int', 'GCV_int')) {
  tmp[c(v, paste0(v,'-sd')),] <- apply(result[,paste0(v,'|lgλ=',lglams_int)], 2, 
                                       function(r) c(mean(r), sqrt(var(r)/MC)))
}
print(round(tmp,3))
}

stopCluster(cl)
