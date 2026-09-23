# Under the degree corrected stochastic block model
library(MASS)
library(RSpectra)
#library(evd)
source("fun_DCSBM.R")
set.seed(1024)

### Simulation 4: Balance Community size
# Settings
k0 <- 3
n <- 200*k0 #let the average size of each block be 200
Pi <- rep(1/k0,k0)

out.in.ratio <- 5
sparsity <- 0.05

P.in <- sparsity*(1+out.in.ratio)
P.bt <- sparsity
P <- (P.in-P.bt)*diag(k0) + P.bt*matrix(1,k0,k0)

Theta <- sample(c(7/11,15/11,1),size = n,replace = T,prob = c(0.1,0.1,0.8))
flag.1 <- Theta==1
eta.n <- sum(flag.1)
eta <- runif(eta.n,0.6,1.4)
Theta[flag.1] <- eta

# Simulation over 20 repititions
K.max <- 10
simul.max <- 20
K.hat.seq <- rep(0,simul.max)
com.dist.seq<- rep(0,simul.max)

simul <- 0
while (simul < simul.max){
  simul <- simul + 1
  # Generating data
  Data.DCSBM <- Adj.Generating.DCSBM(n,Pi,P,Theta,sp = T,del_0d =T)
  Adj <- Data.DCSBM$Adj
  # table(Data.DCSBM$clusters)/n
  
  # proposed Algorithm
  get.est.flag <- FALSE
  temp <- list()
  bug <- 0
  while(!get.est.flag){
    temp <- try(Est.DCSBM(Adj,K.max),silent=TRUE)
    if('try-error' %in% class(temp)) # judge weather error occurs
    {
      bug <- bug+1
      if(bug==5){
        break
      }
      next
    }else{
      get.est.flag <- TRUE
    }
  }
  if(bug==5){
    next
  }
  Para.hat <- temp # e.hat,K
  
  K.hat.seq[simul] <- Para.hat$K
  com.dist.seq[simul] <- com.dist(Data.DCSBM$clusters,Para.hat$e.hat)

}

### Metrics ###
# Community number estimation
K.prop <- mean(as.numeric(K.hat.seq == k0))
K.mean <- mean(K.hat.seq)
# Community detection
com.dist.true <- mean(com.dist.seq[K.hat.seq == k0])
com.dist.all <- mean(com.dist.seq)

simul4.df <- data.frame('k0'=k0,'K.prop'=K.prop,'K.mean'=K.mean,'com.dist.true'=com.dist.true,'com.dist.all'=com.dist.all)
simul4.df