# Under the stochastic block model
library(MASS)
source("fun_SBM.R")
set.seed(2024)

### Simulation 1: Balance Community size
# Settings
k0 <- 3
n <- 200*k0 #let the average size of each block be 200
Pi <- rep(1/k0,k0)
  
out.in.ratio <- 5
sparsity <- 0.05

P.in <- sparsity*(1+out.in.ratio)
P.bt <- sparsity
P <- (P.in-P.bt)*diag(k0) + P.bt*matrix(1,k0,k0)

# Simulation over 20 repititions
K.max <- 10
simul.max <- 20
K.hat.seq <- rep(0,simul.max)
com.dist.seq<- rep(0,simul.max)

simul <- 0
while (simul < simul.max){
  simul <- simul + 1
  ### Generating data ###
  Data.SBM <- Adj.Generating.SBM(n,Pi,P) 
  Adj <- Data.SBM$Adj
  
  ### KPPL algorithm ###
  get.est.flag <- FALSE
  temp <- list()
  while(!get.est.flag){
    temp <- try(Est.SBM(Adj,K.max,tun.c=10),silent=TRUE)
    if('try-error' %in% class(temp)) # judge weather error occurs
    {
      next
    }else{
      get.est.flag <- TRUE
    }
  }
  Para.hat <- temp # e.hat,K
  
  K.hat.seq[simul] <- Para.hat$K
  com.dist.seq[simul] <- com.dist(Data.SBM$clusters,Para.hat$e.hat)
}

### Metrics ###
# Community number estimation
K.prop <- mean(as.numeric(K.hat.seq == k0))
K.mean <- mean(K.hat.seq)
# Community detection
com.dist.true <- mean(com.dist.seq[K.hat.seq == k0])
com.dist.all <- mean(com.dist.seq)

simul1.df <- data.frame('k0'=k0,'K.prop'=K.prop,'K.mean'=K.mean,'com.dist.true'=com.dist.true,'com.dist.all'=com.dist.all)
simul1