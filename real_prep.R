library(readr)
TX1 <- read_csv("TX1day_traj_era5_v10_10_202011.csv")
m <- 121
n0 <- nrow(TX1) / m
logp_range <- matrix(nrow=n0, ncol=2)
for (i in 1:n0) {
  logp_range[i,] <- log(range(TX1$p[(1:m)+(i-1)*m]))
}
idx <- abs(logp_range[,1]-6.3)<0.01 & abs(logp_range[,2]-6.9)<0.01
idx <- (1:n0)[idx]
mydata <- list()
for (i in idx) {
  mydata[[i]] <- TX1[(1:m)+(i-1)*m, 4:6]
}
mydata <- mydata[idx]
save(mydata, file="TX1.RData")
