rm(list=ls())





alpha <- 0.1
T <-60
Nsim <-1000
Tsim <- 1000
rho_ls <- 0.4


rho.grid <- seq(0,1,0.05)
rho.sim <- matrix(0,nrow=length(rho.grid),ncol=Nsim)
rho_q <- rep(0,length(rho.grid))


it <- 1

#for(rho_tmp in rho.grid){
rho_tmp <- 0.7
  Y <- matrix(0,nrow=Tsim,ncol=Nsim)
  e <- matrix(rnorm(Tsim),nrow=Tsim,ncol=Nsim)
  Y[1,]  <- e[1,]
  #Y[2:Tsim,] <- rho_tmp*Y[1:(Tsim-1),] + e[2:Tsim,]
  Y.data < 
  et <- Y[(Tsim-T+1):Tsim,]
  et_1 <- Y[(Tsim-T):(Tsim-1),]
  I <- diag(T)
  P <- matrix(1,nrow=T,ncol=T)
  #rho.sim[it,] <- apply(Y[(Tsim-T):(Tsim-1),]*Y[(Tsim-T+1):Tsim,],2,sum)/apply(Y[(Tsim-T):(Tsim-1),]*Y[(Tsim-T):(Tsim-1),],2,sum)
  rho.sim[it,] <- apply((I-P)%*%et_1*et,2,sum)/apply((I-P)%*%et_1*et_1,2,sum)
  rho_q[it] <- quantile(rho.sim[it,],alpha/2)
  it <- it+1

  
  

index <- sum(rho_ls>rho_q)
rho.grid[index]
rho_q





