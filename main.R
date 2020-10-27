rm(list=ls())



                                        # parellel setup
library(acrt)
library(foreach)
library(doSNOW)
library(tcltk)
cl <- makeSOCKcluster(24)
registerDoSNOW(cl)

                                        #library(doParallel)
                                        #setup parallel backend to use many processors
                                        #cores <- detectCores(all.tests = FALSE, logical = FALSE)
                                        #cl <- makeCluster(5) #not to overload your computer
                                        #registerDoParallel(cl)


                                        # dgp
Par <- c(Tsim <- 1000,
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
p <- length(mu)

seed <- 872013
set.seed(seed)


reject_cmax <- t(rep(0,I))
reject_pp <- t(rep(0,I))
power_cmax <- t(rep(0,41))
power_pp <- t(rep(0,41))
cv_pp <- t(rep(0,I))
cv_max <- t(rep(0,I))
Tw <- t(rep(0,I))
rho_hat <- t(rep(0,I))




                                        # compute quantile of AR coefficient
ar.quantile = function(delta,T,Nsim){
    rho_hat <- rep(0,Nsim)

    for(i in 1:Nsim){
        Y <- rep(0,T)
        e <- matrix(rnorm(T),nrow = T)
        Y[1]  <- e[1]
        for(t in 2:T){
            Y[t] <- rho*u[t-1] + e[t]
        }
        rho_hat[i] <- sum(Y[1:T-1]*Y[2:T])/sum(Y[1:T-1]*Y[1:T-1])
    }
    return(q_ar = quantile(rho_hat,1-delta))
}


cmax.fun = function(alpha,delta,K,T,Nsim){
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
    Time <- seq(t_start,t_end,dt)
                                        #compute x and y
    h <- T*(1-ar.quantile(alpha-delta,T,Nsim))
    i <- 1
    x <- matrix(0,nrow = length(Time),ncol = Nsim)
    x[1,] <- x0
    y <- t(rep(0,Nsim))
    if(h==0){
        for(t in seq((t_start+dt),t_end,dt)){
            i <- i+1
            r1 <- rnorm(Nsim)
            
            x[i,] <- x[i-1,] + sqrt(dt)*r1
        }
    } else {
        for(t in seq((t_start+dt),t_end,dt)){
            i <- i+1
            r1 <- rnorm(Nsim)
            
            x[i,] <- exp(-h*dt)*x[i-1,] + sqrt((c/h*0.5)*(1-exp(-2*h*dt)))*r1
        }
    }
    
    y <- dt*t(rep(1,dim(x)[1]-1))%*%x[2:dim(x)[1],]
    F <- F_infty.fun(x,y,K,dt)
    return(cv = quantile(F,1-delta))
}

F_infty.fun = function(v,y,K,dt){
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






phi <- matrix(rep(0,T*K),nrow = T,ncol = K)
for(t in 1:T){  # basis funcitons
    for(j in 1:(K/2)){
        phi[t,2*j-1] <- sqrt(2)*cos(2*j*pi*t/T)
        phi[t,2*j] <- sqrt(2)*sin(2*j*pi*t/T)
    }
}

                                        # cv_cmax = cmax(alpha,K,10,Nsim,5000)
b <- t(seq(-1,1,0.05))
B <- length(b)

F <- t(rep(0,I))



for(n_b in 1){
    bb <- b[n_b]
    pb <- txtProgressBar(min=1, max=I, style=3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress=progress)
    
                                        #foreach(i=1:I, .packages = 'acrt', .options.snow=opts) %dopar% {
    for(i in 1:I) {
        setTxtProgressBar(pb, i)
        set.seed(seed+i)
        u <- matrix(0,nrow = Tsim)
        e <- matrix(rnorm(T),nrow = Tsim)
        u[1]  <- e[1]
        X <- matrix(rnorm(T)*sig+mu,nrow = T)
        for(t in 2:Tsim){
            u[t] <- rho*u[t-1] + e[t]
        }
        u <- u[(Tsim+1-T):Tsim]
        Y <- X%*%bb + u
        
        ## cmax test
        #et <- (Y[2:T]-X[2:T]*beta0)-mean(Y-X*beta0)
        #et_1 <- (Y[1:T-1]-X[1:T-1]*beta0)-mean(Y-X*beta0)
                                        #rho_hat[i] <- sum(et*et_1)/sum(et_1*et_1)

        cv_max[i] <- optimize(cmax.fun,interval=c(delta_low,alpha),alpha=alpha,K=K,T=T,Nsim=Nsim)
        vhat <- X*(Y-X*beta0)
        F[i] <- (apply(vhat,2,sum))^2*K*(apply((t(phi)%*%vhat)^2,2,sum))^-1
        reject_cmax[i] = F[i]>cv_max[i]

        ## pp test
        cv_pp[i] <- as.numeric(unlist(critical.value(alpha, ar.order.max, bandwidth, ker, R, X, N0, N1, N2, Mp, M1,
                                                     M2, Eicker = FALSE,opt.method.1 = "Nelder-Mead",opt.method.2 = "Nelder-Mead",
                                                     control.1 = list("reltol" = N1^(-.5), "maxit" = dim(X)[1]*20),
                                                     control.2 = list("reltol" = N2^(-.5), "maxit" = dim(X)[1]*30),
                                                     cores = 1, margin = rep(1, length = ar.order.max))[8]))
        
        J <- t(1:(T-1))
        Mn <- floor(T/10)
        w <- 2*(rep(1,length(J)) - J/Mn)*(abs(J/Mn)<1)
                                        #Beta <- solve(t(Xtilde)%*%Xtilde)%*%(t(Xtilde)%*%Y)
                                        #Psi <- matrix(rep(0,(p+1)^2),nrow = p+1, ncol = p+1)
                                        #vhat_pp <- Xtilde*((Y-Xtilde%*%Beta)%*%t(rep(1,dim(beta)[1])))
                                        #for(j in 1:(T-1)){
                                        #    Psi <- Psi + 1/T*w[j]*(t(vhat_pp[(j+1):T,])%*%vhat_pp[1:(T-j),])
                                        #}
                                        #Psi <- t(Psi) + 1/T*(t(vhat_pp)%*%vhat_pp)
        Tw[i] <- as.numeric(unlist(F.type.test.statistic(Y, R, r, X, bandwidth, ker, Eicker = FALSE, cores = 1)))
                                        #(norm(Psi)~=0)*t(Beta[2:length(Beta)]-beta0)*solve(T*t(c(rep(1,p),0,0))*solve(t(Xtilde)%*%Xtilde)*Psi*solve(t(Xtilde)%*%Xtilde)*c(rep(1,p),0,0))*(Beta[2:length(Beta)]-beta0)
        reject_pp[i] <- Tw[i]>cv_pp[i]
    }
    power_cmax[n_b] <- mean(reject_cmax)
    power_pp[n_b] <- mean(reject_pp)
}

stopCluster(cl)

save.image(file="h_estimate.RData")
