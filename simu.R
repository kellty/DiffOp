source('func.R')
B <- 1e3
MC <- 1e3
library(foreach)
library(doParallel)
cl <- makeCluster(100) # detectCores()
registerDoParallel(cl)

bw <- 1e-2
K1 <- function(x,y) { dnorm(x-y, sd=bw) }
LK1 <- function(x,y) { ((x-y)^2/bw^4 - 1/bw^2) * dnorm(x-y, sd=bw) }
# K2 <- function(x1,y1,x2,y2) { K1(x1,x2) * K1(y1,y2) }
# K2bd <- function(z1,y1,z2,y2) { (z1==z2) * K1(y1,y2) }
# LK2 <- function(x1,y1,x2,y2) { K1(x1,x2) * LK1(y1,y2) }
# LK2bd <- function(z1,y1,z2,y2) { (z1==z2) * LK1(y1,y2) }

p <- 10;  p2 <- p^2
phi <- function(x, k) { sqrt(2) * cos(k*pi*x) }
Lphi <- function(x, k) { sqrt(2) * (k*pi)^2 * cos(k*pi*x) }

idx <- cbind(rep(1:p,each=p), rep(1:p,times=p))
K <- KL <- diag(p2)
# for (i in 1:p2) {
#   j <- idx[i,1];  k <- idx[i,2]
#   K[i,] <- foreach(i_=1:p2, .combine=c, .packages='cubature') %dopar% {
#     j_ <- idx[i_,1];  k_ <- idx[i_,2]
#     hcubature(function(xyxy) {
#       x1 <- xyxy[1];  y1 <- xyxy[2]; x2 <- xyxy[3]; y2 <- xyxy[4]
#       K1(x1,x2) * K1(y1,y2) * Lphi(x1,k_) * phi(y1,j_) * Lphi(x2,k) * phi(y2,j)
#     }, rep(0,4), rep(1,4))$integral + 
#       (phi(0,k_) * phi(0,k) + phi(1,k_) * phi(1,k)) * hcubature(function(yy) {
#         y1 <- yy[1];  y2 <- yy[2]
#         K1(y1,y2) * phi(y1,j_) * phi(y2,j)
#       }, c(0,0), c(1,1))$integral
#   }
#   KL[i,] <- foreach(i_=1:p2, .combine=c, .packages='cubature') %dopar% {
#     j_ <- idx[i_,1];  k_ <- idx[i_,2]
#     hcubature(function(xyxy) {
#       x1 <- xyxy[1];  y1 <- xyxy[2]; x2 <- xyxy[3]; y2 <- xyxy[4]
#       K1(x1,x2) * LK1(y1,y2) * Lphi(x1,k_) * phi(y1,j_) * Lphi(x2,k) * phi(y2,j)
#     }, rep(0,4), rep(1,4))$integral + 
#       (phi(0,k_) * phi(0,k) + phi(1,k_) * phi(1,k)) * hcubature(function(yy) {
#         y1 <- yy[1];  y2 <- yy[2]
#         LK1(y1,y2) * phi(y1,j_) * phi(y2,j)
#       }, c(0,0), c(1,1))$integral
#   }
# }
Kx <- sapply(1:p, function(k) sapply(1:p, function(k_) 
  integrate(function(x1vec) sapply(x1vec, function(x1)
    integrate(function(x2) K1(x1,x2) * Lphi(x1,k_) * Lphi(x2,k), 0,1)$value), 0,1)$value))
Ky <- sapply(1:p, function(j) sapply(1:p, function(j_) 
  integrate(function(y1vec) sapply(y1vec, function(y1)
    integrate(function(y2) K1(y1,y2) * phi(y1,j_) * phi(y2,j), 0,1)$value), 0,1)$value))
LKy <- sapply(1:p, function(j) sapply(1:p, function(j_) 
  integrate(function(y1vec) sapply(y1vec, function(y1)
    integrate(function(y2) LK1(y1,y2) * phi(y1,j_) * phi(y2,j), 0,1)$value), 0,1)$value))
for (i in 1:p2) for (i_ in 1:p2) {
  j <- idx[i,1];  k <- idx[i,2];  j_ <- idx[i_,1];  k_ <- idx[i_,2]
  K[i,i_] <- (Kx[k,k_] + phi(0,k_) * phi(0,k) + phi(1,k_) * phi(1,k)) * Ky[j,j_]
  KL[i,i_] <- (Kx[k,k_] + phi(0,k_) * phi(0,k) + phi(1,k_) * phi(1,k)) * LKy[j,j_]
}

D1_eig <- function(omega) { omega^2 + ((1:p)*pi)^2 }
D_eig <- function(theta) { theta * ((1:p)*pi)^2 }
D_mat <- function(theta) { diag(D_eig(theta)) }
theta_ini = 1;  theta_len <- length(theta_ini)

for(n in c(200,400)) for(SNR in c(8,3,1)) for(case in c(0:3,4,6)) {
  idx <- paste0('_n',n,'_SNR',SNR,'_Case',case)
if (n==200 & SNR==8) { lglams <- 2.3+(-2:2)*0.2;  omega_d <- 0.25 }
if (n==200 & SNR==3) { lglams <- c(0:2,3+(-2:2)*0.2,4:6);  omega_d <- 0.42 }
if (n==200 & SNR==1) { lglams <- 3.7+(-2:2)*0.2;  omega_d <- 0.83 }
if (n==400 & SNR==8) { lglams <- 2+(-2:2)*0.2;  omega_d <- 0.21 }
if (n==400 & SNR==3) { lglams <- 2.7+(-2:2)*0.2;  omega_d <- 0.35 }
if (n==400 & SNR==1) { lglams <- 3.5+(-2:2)*0.2;  omega_d <- 0.66 }
  omega <- case * omega_d;  print(paste(Sys.time(), idx, 'ω:',omega))
result <- foreach(mc=1:MC, .combine=rbind) %dopar% {
  t0 <- proc.time()[3]
  set.seed(999+(mc+55)^2)
U_sd <- sapply(1:p, function(k) 1/k^3)
U_coef <- sapply(1:p, function(k) runif(n,-1,1)*sqrt(3)) %*% diag(U_sd)
F_sd <- U_sd * D1_eig(omega)
F0_coef <- U_coef %*% diag(D1_eig(omega))
eps_sd <- rep(sqrt(mean(F_sd^2)) / SNR, p)
eps_coef <- sapply(1:p, function(k) rnorm(n, sd=eps_sd[k]))
F_coef <- F0_coef + eps_coef

fit0 <- fit_para(F_coef, U_coef, D_mat, theta_ini)
res0 <- F_coef - fit0$Ftheta_coef
tmp <- c(sum(F0_coef^2), sum(eps_coef^2), 
         sum(res0^2), sum((F0_coef - fit0$Ftheta_coef)^2), fit0$iter, fit0$theta)

for (lam in 10^lglams) {
fit1 <- fit_RKHS(F_coef, U_coef, K, KL, lam)
res1 <- F_coef - fit1$Fhat_coef

stat <- test.stat(fit1$H, res0)
stat_btsp1 <- test.btsp(fit1$H, res1, B=B)
stat_btsp0 <- test.btsp(fit1$H, res0, B=B)
tmp <- c(tmp, sum(res1^2), sum((F0_coef - fit1$Fhat_coef)^2), fit1$gcv, 
         test.pval(stat,stat_btsp1), test.pval(stat,stat_btsp0), 
         test.pval(stat,c(stat_btsp0,rep(stat_btsp1,2))), stat, stat_btsp1, stat_btsp0)
}
tmp <- c(proc.time()[3]-t0, tmp)
names(tmp) <- c('time', 'TSS', 'ESS', 'RSS_para', 'ESS_para', 'iter', paste0('θ^',1:theta_len),
    paste0(c('RSS_nonp', 'ESS_nonp', 'GCV', 'pval_nonp', 'pval_para', 'pval_mix', 'test_stat',
             paste0(rep(c('b_nonp_','b_para_'),each=B), 1:B)), '|lgλ=', rep(lglams, each=7+2*B) ))
return(tmp)
}
save(result, file=paste0('result',idx,'.RData'))
print(round(apply(result[,c(2:5,6:(6+theta_len))], 2, function(r) c(mean(r), sqrt(var(r)/MC))),5))
tmp <- matrix(nrow=10, ncol=length(lglams))
rownames(tmp) <- c('lgλ', 'ESS_nonp', 'ESS_nonp-sd', 'GCV', 'GCV-sd', 
                   'RSS_nonp', 'RSS_nonp-sd', 'rej_nonp', 'rej_para', 'rej_mix')
tmp['lgλ',] <- lglams
for (v in c('nonp','para','mix')) {
  tmp[paste0('rej_',v),] <- 100*colMeans(result[,paste0('pval_',v,'|lgλ=',lglams)]<0.05)
}
for (v in c('ESS_nonp', 'RSS_nonp', 'GCV')) {
  tmp[c(v, paste0(v,'-sd')),] <- apply(result[,paste0(v,'|lgλ=',lglams)], 2, 
                                       function(r) c(mean(r), sqrt(var(r)/MC)))
}
print(round(tmp,3))
}

stopCluster(cl)

pdf(paste0('stat.pdf'), height=3*6, width=3*4, family = "Times New Roman")
par(mfrow=c(6,4))
for(n in c(200,400)) for(SNR in c(8,3,1)) {
  if (n==200 & SNR==8) { lglam <- 2.3;  omega_d <- 0.25 }
  if (n==200 & SNR==3) { lglam <- 3;  omega_d <- 0.42 }
  if (n==200 & SNR==1) { lglam <- 3.7;  omega_d <- 0.83 }
  if (n==400 & SNR==8) { lglam <- 2;  omega_d <- 0.21 }
  if (n==400 & SNR==3) { lglam <- 2.7;  omega_d <- 0.35 }
  if (n==400 & SNR==1) { lglam <- 3.5;  omega_d <- 0.66 }
for (case in 0:3) {
  load(paste0('result','_n',n,'_SNR',SNR,'_Case',case,'.RData'))
  plot(density(result[,paste0('test_stat|lgλ=',lglam)]), lwd=2, 
       main=bquote(n == .(n) ~','~ SNR == .(SNR) ~','~ omega == .(case*omega_d)))
  for (mc in c(233+(1:5))) {
    lines(density(result[mc,paste0('b_nonp_',1:B,'|lgλ=',lglam)]), lwd=2, lty=2, col=2)
    lines(density(result[mc,paste0('b_para_',1:B,'|lgλ=',lglam)]), lwd=2, lty=4, col=4)
  }
  legend('topright', c("test stat", "btsp (res1)", "btsp (res0)"), lwd=2, lty=c(1,2,4), col=c(1,2,4))
  }
}
dev.off()

library(plot3D)
G1 <- function(x,y,omega) {
  if (omega==0) { pmin(x,y) - x*y }
  else { sin(omega*pmin(x,y)) * sin(omega*(1-pmax(x,y))) / (omega * sin(omega)) }
}
m <- 100
x <- rep(((1:m) - 0.5) / m, times=m)
y <- rep(((1:m) - 0.5) / m, each=m)
z0 <- G1(x,y,0)
pdf("GreenFun.pdf", height=6, width=12, family = "Times New Roman")
par(mfcol=c(2,5),mar=c(1,1,2,3),oma=c(1,1,0,2))
for (omega in (0:4)*0.5) {
  z <- G1(x,y,omega)
  z_change <- z - z0
  # print(summary(z_change))
  scatter3D(x, y, z, pch='.', theta=-45, phi=0, zlab='', zlim=c(0,0.3), main='G',
            surf=list(x=matrix(x,nrow=m),
                      y=matrix(y,nrow=m),
                      z=matrix(z,nrow=m)))
  legend('bottomright', legend=bquote(omega==.(omega)), box.col = "white")
  scatter3D(x, y, z_change, pch='.', theta=-45, phi=0, zlab='',
            zlim=c(0,0.07), main=expression(Delta~'G'),
            surf=list(x=matrix(x,nrow=m),
                      y=matrix(y,nrow=m),
                      z=matrix(z_change,nrow=m)))
  legend('bottomright', legend=bquote(omega==.(omega)), box.col = "white")
}
dev.off()
