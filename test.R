rm(list=ls())





alpha <- 0.1
T <-60
Nsim <-1000
Tsim <- 500
rho_ls <- 0.4


rho.grid <- seq(0,1,0.05)
rho.sim <- matrix(0,nrow=length(rho.grid),ncol=Nsim)
rho_q <- rep(0,length(rho.grid))


it <- 1

for(rho_tmp in rho.grid){
  X <- matrix(rnorm(Tsim),nrow=Tsim,ncol=Nsim)
  e <- matrix(rnorm(Tsim),nrow=Tsim,ncol=Nsim)
  Y<- rho_tmp*X+ e
  Y <- Y[(Tsim-T+1):Tsim,]
  X <- X[(Tsim-T+1):Tsim,]
  #rho.sim[it,] <- apply(Y[(Tsim-T):(Tsim-1),]*Y[(Tsim-T+1):Tsim,],2,sum)/apply(Y[(Tsim-T):(Tsim-1),]*Y[(Tsim-T):(Tsim-1),],2,sum)
  rho.sim[it,] <- (apply(X*Y,2,mean)-apply(X,2,mean)*apply(Y,2,mean))/(apply(X^2,2,mean)-(apply(X,2,mean))^2)
  rho_q[it] <- quantile(rho.sim[it,],alpha/2)
  it <- it+1
}

index <- sum(rho_ls>rho_q)
rho.grid[index]






