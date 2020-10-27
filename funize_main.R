rm(list=ls())



                                        # parellel setup
library(acrt)
library(foreach)
library(doSNOW)
library(tcltk)



                                        #library(doParallel)
                                        #setup parallel backend to use many processors
                                        #cores <- detectCores(all.tests = FALSE, logical = FALSE)
                                        #cl <- makeCluster(5) #not to overload your computer
                                        #registerDoParallel(cl)


     
Par <- c(Tsim <- 300,
         T <- 50,
         I <- 1000, # replications
         Nsim <- 3000, # number of simulations for calculating critival values
         alpha <- 0.05,
         delta_low <- 0.01,
         beta0 <- 0.1,
         rho <- 0.9,
         mu <- 0,
         sig <- 1,
         K <- 24,
         p <- length(mu),
                                        #h <- 10,
                                        # Parameters for PP test
         ar.order.max <- 1,
         bandwidth <- floor(T/10),
         ker <- "Bartlett",
         R <- diag(length(beta0)),
         r <- beta0,
         Mp <- 500,
         M1 <- 10,
         M2 <- 2,
         N0 <- 200,
         N1 <- 1000,
         N2 <- 5000
         )

## 0.05 quantiles from Andrews(1993)
rho.andrews93 <- c(-0.211,-0.164,-0.116,-0.068,
              -0.019,0.03,0.08,0.13,0.181,
              0.233,0.285,0.338,0.392,0.446,
              0.502,0.559,0.618,0.678,0.719,
              0.784,0.813,0.846,0.891,0.91,0.944)



seed <- 872013
set.seed(seed)




                                        # compute quantile of AR coefficient
ar.rho_max <- function(rho_ls,delta,T,Nsim){
    rho.grid <- seq(0,1,0.05)
    rho.sim <- matrix(0,nrow=length(rho.grid),ncol=Nsim)
    rho_q <- rep(0,length(rho.grid))
    it <- 1
    for(rho_tmp in rho.grid){
        Y <- matrix(0,nrow=T,ncol=Nsim)
        e <- matrix(rnorm(T*Nsim),nrow = T,ncol=Nsim)
        Y[1,]  <- e[1,]
        Y[2:T,] <- rho_tmp*Y[1:(T-1),] + e[2:T,]
        rho_q[it] <- quantile(apply(Y[1:T-1,]*Y[2:T,],2,sum)/apply(Y[1:T-1,]*Y[1:T-1,],2,sum),delta/2)
      it <- it +1
    }
    index <- sum(rho_ls>rho_q)
    return(rho.grid[index])
}

cmax.fun <- function(h,alpha,delta,K,T,Nsim){
    t_start <- 0       #simulation start time
    t_end <- 1          #simuation end time
    dt <- 0.005         #time step
                                        #h = 0.1            relaxation time
    c <- 1                 #diffusion constant
    x0 <- 0               #initial value for stochastic variable x
                                        #mu = 0               #mean of stochatic process x
                                        # start_dist = -2.0;    %start of OU pdf
                                        # end_dist = 2.0;       %end of OU pdf
                                        #time
    Time <- seq(t_start,t_end,by=dt)
                                        #compute x and y
    i <- 1
    x <- matrix(0,nrow = length(Time),ncol = Nsim)
    x[1,] <- x0
    y <- t(rep(0,Nsim))
    if(h==0){
        for(t in seq((t_start+dt),t_end,by=dt)){
            i <- i+1
            r1 <- rnorm(Nsim)
            
            x[i,] <- x[i-1,] + sqrt(dt)*r1
        }
    } else {
        for(t in seq((t_start+dt),t_end,by=dt)){
            i <- i+1
            r1 <- rnorm(Nsim)
            
            x[i,] <- exp(-h*dt)*x[i-1,] + sqrt((c/h*0.5)*(1-exp(-2*h*dt)))*r1
        }
    }
    
    y <- dt*t(rep(1,dim(x)[1]-1))%*%x[2:dim(x)[1],]
    F <- F_infty.fun(x,y,K,dt)
    return(cv = quantile(F,1-delta))
}

F_infty.fun <- function(v,y,K,dt){
                                        # y: simulations of y vectors
                                        # v: simulations of v vectors
    T <- dim(v)[1]
    phi <- matrix(rep(0,T*K),nrow = T, ncol = K)
                                        # basis funcitons
    for(t in 1:T){
        for(j in 1:(K/2)){
            phi[t,2*j-1] <- sqrt(2)*cos(2*j*pi*t/T)
            phi[t,2*j] <- sqrt(2)*sin(2*j*pi*t/T)
        }
    }
    F <- (y^2)*K*((apply((t(phi)%*%v*dt)^2,2,sum))^-1) # simulations of F Nsim times
}


# generate basis functions
phi <- matrix(rep(0,T*K),nrow = T,ncol = K)
for(t in 1:T){  # basis funcitons
    for(j in 1:(K/2)){
        phi[t,2*j-1] <- sqrt(2)*cos(2*j*pi*t/T)
        phi[t,2*j] <- sqrt(2)*sin(2*j*pi*t/T)
    }
}


b <- t(seq(-1,1,by=0.05))
B <- length(b)
power <- matrix(0,nrow = 2, ncol=B)
#F <- t(rep(0,I))



sim.power <- function(sim.Y, sim.X){
  et <- sim.Y[2:T] -sim.X[2:T]*beta0
  et_1 <- sim.Y[1:T-1] -sim.X[1:T-1]*beta0
  idm <- diag(T-1)
  #X2 <- matrix(c(rep(1,T-1),sim.X[2:T]),nrow=T-1)
  #P <- X2%*%solve(t(X2)%*%X2)%*%t(X2)
  #r_ls <- t(et_1)%*%et/(t(et_1)%*%et_1)
  #grid <- c(seq(0,0.9,0.05),0.93,0.95,0.97,0.99,0.995,1)
  #id <- sum(c(r_ls)>rho.andrews93)
  #h <- T*(1-grid[id])  #ar.rho_max(rho_ls,alpha-delta,T,Nsim))
  h <- 30
  cv_max <- optimize(cmax.fun,interval=c(delta_low,alpha),h=h,alpha=alpha,K=K,T=T,Nsim=Nsim)[[2]]
  vhat <- sim.X*(sim.Y-sim.X*beta0)
  F <- (sum(vhat))^2*K*(sum((t(phi)%*%vhat)^2))^-1
  reject_cmax <- F>cv_max

  cv_pp <- critical.value(alpha, ar.order.max, bandwidth, ker, R, sim.X, N0, N1, N2, Mp, M1,
                                               M2, Eicker = FALSE,opt.method.1 = "Nelder-Mead",opt.method.2 = "Nelder-Mead",
                                               control.1 = list("reltol" = N1^(-.5), "maxit" = dim(sim.X)[1]*20),
                                               control.2 = list("reltol" = N2^(-.5), "maxit" = dim(sim.X)[1]*30),
                                               cores = 1, margin = rep(1, length = ar.order.max))[[8]]

  Tw <- as.numeric(F.type.test.statistic(sim.Y, R, r, sim.X, bandwidth, ker, Eicker = FALSE, cores = 1))
  reject_pp <- Tw>cv_pp
  output <- c(rc=as.numeric(reject_cmax),rp=as.numeric(reject_pp))
  #output <- list(reject_cmax,reject_pp)
  return(output)
}



cl <- makeSOCKcluster(2)
registerDoSNOW(cl)


# Compute power
for(n_b in c(1,21,41)){
    bb <- b[n_b]
    # DGP
    u <- matrix(0,nrow = Tsim,ncol = I)
    e <- matrix(rnorm(T*I),nrow = Tsim,ncol = I)
    u[1,]  <- e[1,]
    X <- matrix(rnorm(T*I)*sig+mu,nrow = T,ncol = I)
    for(t in 2:Tsim){
      u[t,] <- rho*u[t-1,] + e[t,]
    }
    u <- u[(Tsim+1-T):Tsim,]
    Y <- X*bb + u
    #vhat <- X*(Y-X*beta0)
    Data <- matrix(c(Y,X),nrow=T)
    
    pb <- txtProgressBar(min=1, max=I, style=3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress=progress)

    reject <- foreach(i=1:I, .packages = 'acrt', .combine = 'cbind', .options.snow=opts) %dopar% {
      setTxtProgressBar(pb, i)
      #i<- 3
      sim.power(matrix(Data[,i],nrow=T),matrix(Data[,I+i],nrow=T))
    }
    
    power[,n_b] <- apply(reject,1,mean)
}

stopCluster(cl)

save.image(file="Bonf_test.RData")
