load('TX1.RData')
xmin <- 6.3;  xmax <- 6.9;  xrange <- xmax - xmin
ngrid <- 100; dx <- xrange/ngrid; xgrid <- xmin + dx * ((1:ngrid)-0.5)
T_real_val <- T_pot_val <- NULL;  n <- 0
for (i in 1:length(mydata)) {
  TX1 <- mydata[[i]]
  T_pot <- smooth.spline(log(TX1$p),TX1$TH)
  if (mean(abs(predict(T_pot, xgrid, deriv=1)$y)) < 50) {
    n <- n + 1
    T_pot_val <- rbind(T_pot_val, predict(T_pot, xgrid)$y)
    T_real_val <- rbind(T_real_val, predict(smooth.spline(log(TX1$p),TX1$T), xgrid)$y)
  }
}

library(RColorBrewer)
cols <- brewer.pal(12, "Set3")
pdf('TX1.pdf', height=4, width=4*2)
par(mfrow=c(1,2))
i0 <- 23
plot(xgrid, colMeans(T_real_val), type = "l", lty=2, lwd=2, main="real temperature",
     ylim=range(T_real_val[i0+(1:12),]), ylab=expression(T[real]), xlab=expression(log(p)))
for (i in 1:12) { lines(xgrid, T_real_val[i0+i,], col=cols[i]) }
plot(xgrid, colMeans(T_pot_val), type = "l", lty=2, lwd=2, main="potential temperature",
     ylim=range(T_pot_val[i0+(1:12),]), ylab=expression(T[potential]), xlab=expression(log(p)))
for (i in 1:12) { lines(xgrid, T_pot_val[i0+i,], col=cols[i]) }
dev.off()


source('func.R')

p <- 11;  p2 <- p^2
phi <- function(x, k) {
  if (k==1) return(cos(0*x) / sqrt(xrange))
  if (k%%2) return(sqrt(2/xrange) * cos((k-1)*pi*(x-xmin)/xrange))
  else return(sqrt(2/xrange) * sin(k*pi*(x-xmin)/xrange))
}
Lphi <- function(x, k) {
  if (k==1) return(0*x)
  if (k%%2) return(- sqrt(2/xrange) * (k-1)*pi/xrange * sin((k-1)*pi*(x-xmin)/xrange))
  else return(sqrt(2/xrange) * k*pi/xrange * cos(k*pi*(x-xmin)/xrange))
}
D0mat <- matrix(0, nrow=p, ncol=p)
for (j in 1:(p%/%2)) {
  D0mat[2*j+1,2*j] <- j*2*pi/xrange
  D0mat[2*j,2*j+1] <- - j*2*pi/xrange
}

phi_val <- sapply(xgrid, function(x) sapply(1:p, function(k) phi(x,k)))
k_pois <- 0.286
u_coef <- T_pot_val %*% t(phi_val * dx)
f_coef <- T_real_val %*% diag((1000/exp(xgrid))^k_pois) %*% t(phi_val * dx) %*% t(D0mat)
u_coef <- scale(u_coef, center=T,scale=F)
f_coef <- scale(f_coef, center=T,scale=F)

D_mat <- function(theta) { theta * D0mat }
theta_ini = 1;  theta_len <- length(theta_ini)
fit0 <- fit_para(f_coef, u_coef, D_mat, theta_ini, lr=1e-6, tol=1e-9)
res0 <- f_coef - fit0$Ftheta_coef
res00 <- f_coef - u_coef %*% t(D_mat(1))
print(c(sum(f_coef^2), sum(res00^2), sum(res0^2), fit0$iter, fit0$theta))

bw <- 1e-2
K1 <- function(x,y) { dnorm(x-y, sd=bw) }
LK1 <- function(x,y) { ((x-y)^2/bw^4 - 1/bw^2) * dnorm(x-y, sd=bw) }
# K2 <- function(x1,y1,x2,y2) { K1(x1,x2) * K1(y1,y2) }
# K2bd <- function(z1,y1,z2,y2) { (z1==z2) * K1(y1,y2) }
# LK2 <- function(x1,y1,x2,y2) { K1(x1,x2) * LK1(y1,y2) }
# LK2bd <- function(z1,y1,z2,y2) { (z1==z2) * LK1(y1,y2) }

idx <- cbind(rep(1:p,each=p), rep(1:p,times=p))
K <- KL <- diag(p2)
# for (i in 1:p2) {
#   j <- idx[i,1];  k <- idx[i,2]
#   K[i,] <- foreach(i_=1:p2, .combine=c, .packages='cubature') %dopar% {
#     j_ <- idx[i_,1];  k_ <- idx[i_,2]
#     hcubature(function(xyxy) {
#       x1 <- xyxy[1];  y1 <- xyxy[2]; x2 <- xyxy[3]; y2 <- xyxy[4]
#       K1(x1,x2) * K1(y1,y2) * Lphi(x1,k_) * phi(y1,j_) * Lphi(x2,k) * phi(y2,j)
#     }, rep(xmin,4), rep(xmax,4))$integral + 
#       (phi(xmin,k_) * phi(xmin,k)) * hcubature(function(yy) {
#         y1 <- yy[1];  y2 <- yy[2]
#         K1(y1,y2) * phi(y1,j_) * phi(y2,j)
#       }, rep(xmin,2), rep(xmax,2))$integral
#   }
#   KL[i,] <- foreach(i_=1:p2, .combine=c, .packages='cubature') %dopar% {
#     j_ <- idx[i_,1];  k_ <- idx[i_,2]
#     hcubature(function(xyxy) {
#       x1 <- xyxy[1];  y1 <- xyxy[2]; x2 <- xyxy[3]; y2 <- xyxy[4]
#       K1(x1,x2) * LK1(y1,y2) * Lphi(x1,k_) * phi(y1,j_) * Lphi(x2,k) * phi(y2,j)
#     }, rep(xmin,4), rep(xmax,4))$integral + 
#       (phi(xmin,k_) * phi(xmin,k)) * hcubature(function(yy) {
#         y1 <- yy[1];  y2 <- yy[2]
#         LK1(y1,y2) * phi(y1,j_) * phi(y2,j)
#       }, rep(xmin,2), rep(xmax,2))$integral
#   }
# }
Kx <- sapply(1:p, function(k) sapply(1:p, function(k_) 
  integrate(function(x1vec) sapply(x1vec, function(x1)
    integrate(function(x2) K1(x1,x2) * Lphi(x1,k_) * Lphi(x2,k), xmin,xmax)$value), xmin,xmax)$value))
Ky <- sapply(1:p, function(j) sapply(1:p, function(j_) 
  integrate(function(y1vec) sapply(y1vec, function(y1)
    integrate(function(y2) K1(y1,y2) * phi(y1,j_) * phi(y2,j), xmin,xmax)$value), xmin,xmax)$value))
LKy <- sapply(1:p, function(j) sapply(1:p, function(j_) 
  integrate(function(y1vec) sapply(y1vec, function(y1)
    integrate(function(y2) LK1(y1,y2) * phi(y1,j_) * phi(y2,j), xmin,xmax)$value), xmin,xmax)$value))
for (i in 1:p2) for (i_ in 1:p2) {
  j <- idx[i,1];  k <- idx[i,2];  j_ <- idx[i_,1];  k_ <- idx[i_,2]
  K[i,i_] <- (Kx[k,k_] + phi(0,k_) * phi(0,k) + phi(1,k_) * phi(1,k)) * Ky[j,j_]
  KL[i,i_] <- (Kx[k,k_] + phi(0,k_) * phi(0,k) + phi(1,k_) * phi(1,k)) * LKy[j,j_]
}

set.seed(1234)
B <- 1000
for (lam in 10^seq(0.9,1.9,0.2)) {
  fit1 <- fit_RKHS(f_coef, u_coef, K, KL, lam)
  res1 <- f_coef - fit1$Fhat_coef
  stat <- test.stat(fit1$H, res0)
  stat0 <- test.stat(fit1$H, res00)
  stat_btsp1 <- test.btsp(fit1$H, res1, B=B)
  stat_btsp0 <- test.btsp(fit1$H, res0, B=B)
  stat_btsp00 <- test.btsp(fit1$H, res00, B=B)
  print(c(log10(lam), sum(res1^2), fit1$gcv, 
          test.pval(stat,stat_btsp1), test.pval(stat,stat_btsp0), 
          test.pval(stat,c(rep(stat_btsp1,2),stat_btsp0)),
          test.pval(stat0,stat_btsp1), test.pval(stat0,stat_btsp00), 
          test.pval(stat0,c(rep(stat_btsp1,2),stat_btsp00))))
}
