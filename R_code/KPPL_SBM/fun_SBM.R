### Generate adjacency matrix for SBM ###
Adj.Generating.SBM <- function(
    n,
    clusters.distribution,
    connection.matrix,
    del.0d = TRUE){
  
  num.clusters <- length(clusters.distribution)
  label.clusters <- 1:num.clusters
  clusters <- sample(label.clusters,size=n,replace=T,
                     prob=clusters.distribution)
  
  RM <- matrix(runif(n*n),n,n)
  RM[lower.tri(RM)] <- t(RM)[lower.tri(RM)] #Symmetrization
  detM <- connection.matrix[clusters,clusters]
  Adj <- (RM < detM) + 0
  diag(Adj) <- 0
  if(del.0d){
    Ind.0 <- apply(Adj,1,sum)!=0
    Adj <- Adj[Ind.0,Ind.0]
    clusters <- clusters[Ind.0] 
  }
  
  re <- list(clusters,Adj)
  names(re) <- c('clusters','Adj')
  return(re)
}

rep.nrows<-function(x,nrows){
  re <- matrix(rep(x,nrows),nrows,length(x),byrow = T)
  return(re)
}

### The proposed method ###
Est.SBM <- function(Adj,K,tun.c=10,del.0d = TRUE){
  if(del.0d){
    Ind.0 <- apply(Adj,1,sum)!=0
    Adj <- Adj[Ind.0,Ind.0]
  }
  ### Spectral clustering ###
  # the method of Spectral Clustering with Permutations (SCP)
  n<-nrow(Adj)
  rho<-sum(Adj)/(n*(n-1))
  pAdj<-Adj+as.numeric(0.025*rho)*matrix(1,nrow=n,ncol=n) #purmutated Adj see Amini et al. 2013
  # compute the Laplace Matrix
  D <- apply(pAdj,1,sum)^(0.5)
  D.inv <- D^(-1)
  G <- diag(D)
  G.inv <- diag(D.inv)
  L <- G.inv%*%pAdj%*%G.inv
  
  # Take eigenvalues to reduce dimension
  eig <- eigen(L,symmetric = T)
  
  # Sort eigenvalues by its absolute value
  eig.sort <- sort(abs(eig$values), decreasing = TRUE, method = "shell", index.return = TRUE)
  values <- eig$values[eig.sort$ix]
  vectors <- as.matrix(eig$vectors)[,eig.sort$ix]
  
  ### Initialization ###
  e.hat <- kmeans(vectors[,1:K],K)$cluster
  
  ### Outer iterations ###
  outer.iter.max = 25
  inner.iter.max = 50
  converge.thresh = 0.01
  
  Iter <- numeric(outer.iter.max+1)
  Iter[1] <- inner.iter.max
  break.outer.count <- 0
  
  for(i in 1:outer.iter.max){
    
    ### Inner iterations(EM) ###
    inner.count <- 1 #the counter of the iterations times of 
    break.flag <- FALSE #the sign of convenience
    
    K.I <- diag(K)
    e.OneHot <- as.matrix(K.I[e.hat, ]) #n*K
    n.vec <- as.matrix(apply(e.OneHot,2,sum)) # 1*K, size of communities
    
    Pi.hat <- n.vec/sum(n.vec)
    
    if (K==1){
      diag.minus <- n.vec
    }else{
      diag.minus <- diag(as.vector(n.vec))
    }
    n.mat <- n.vec %*% t(n.vec) # number of possible edges
    # diag(n.mat) <- diag(n.mat)/2
    m.mat <- t(e.OneHot)%*%Adj%*%e.OneHot #k*k
    # diag(m.mat) <- diag(m.mat)/2
    
    P.hat <- m.mat / n.mat
    # Here we accomplish initialization, we have e.hat,Pi.hat,P.hat
    
    Adj.e.mat <- Adj%*%e.OneHot #n*K
    
    while((!break.flag)&&(inner.count<=inner.iter.max)){
      ### E-step ###
      log.prod <- Adj.e.mat%*%t(log(P.hat+1e-30))+(rep.nrows(n.vec,n)-Adj.e.mat)%*%t(log(1-P.hat+1e-30)) # n*K
      log.tau0.mat <- log.prod+rep.nrows(log(Pi.hat),n)
      log.tau1.mat <- t(apply(log.tau0.mat,1,function(x) x-max(x)))
      tau0.mat <- exp(log.tau1.mat)
      tau <- t(apply(tau0.mat,1,function(x) x/(sum(x)+1e-30)))
      tau.n <- colSums(tau)
      
      ### M-step ###
      #  Calculate Pi.hat with selection of tuning parameter
      multi <- seq(from=0.1, to=2, by=0.1)
      gamma0 <- sqrt(log(n)/n)
      gamma.seq <- multi*gamma0
      
      Pi.update <- t(rep.nrows(1/(1-K*gamma.seq),K))*(rep.nrows(tau.n/n,length(multi))-gamma.seq) # length(multi) * K; Pi.update <- 1/(1-K*gamma)*(tau.n/n-gamma)
      Pi.update[Pi.update<0] <- 0 # nrows = length(multi)
      Pi.mat <- t(apply(Pi.update,1,function(x) x/sum(x)))
      
      bic1.mat <- exp(log.prod)
      bic2.vec <- log(bic1.mat %*% t(Pi.mat)) # n*length(multi)
      bic.term1 <- -2*colSums(bic2.vec)
      
      K.Pi <- rowSums(Pi.mat != 0)
      # tun.c <- 10 # small c large K
      bic.term2 <- tun.c*n*log(K.Pi)+K.Pi*(K.Pi+1)/2*log(n)
      
      cbic.gamma <- bic.term1+bic.term2
      
      # plot(bic.term1)
      # points(cbic.gamma,col='blue')
      gamma.id <- which.min(cbic.gamma)
      Pi.hat <- Pi.mat[gamma.id,]
      print(Pi.hat)
      
      # Calculate P.hat(which not dependent on gamma)
      P.num <- t(tau)%*%Adj.e.mat
      P.den <- tau.n%*%t(n.vec)
      P.update <- P.num/(P.den+1e-30)
      P.update[P.update>1] <- 1
      P.hat <- P.update
      
      # break if there exists 0 in Pi
      break.flag <- sum(Pi.hat==0)>0
      if(break.flag){
        break.outer.count <- 0
        # remove the parameter related to the removed community
        remove.id <- which(Pi.hat==0)
        tau <- tau[,-remove.id]
        Pi.hat <- Pi.hat[-remove.id]
        P.hat <- P.hat[-remove.id,-remove.id]
      }else{
        # otherwise compute the new value of likelihood
        Q.update <- bic.term1[gamma.id]
        if (inner.count > 1){
          #compute the new relative changes in likelihood
          delta <- (Q.update - Q.old)/(Q.old+1e-30)
          #break if converged
          break.flag <- delta < converge.thresh
        }
        Q.old <- Q.update
      }
      inner.count <- inner.count+1
    }
    
    # update number of community K
    K <- length(Pi.hat)
    
    # update column label e
    tau.Adj.mat <- t(Adj)%*%tau
    e.update <- tau.Adj.mat%*%log(P.hat+1e-30)+(rep.nrows(colSums(tau),n)-tau.Adj.mat)%*%log(1-P.hat+1e-30)
    e.hat <- apply(e.update,1,which.max)
    e.hat <- fix.e(K,e.hat)
    
    # record the number of iterations
    Iter[i+1] <- inner.count
    # break rule
    if (inner.count<=3) {
      break.outer.count <- break.outer.count+1
    }
    if (break.outer.count==3){break}
  }
  
  # Output
  re<-list(e.hat,K)
  names(re) <- c('e.hat','K')
  return(re)
}

fix.e <- function(K,e.hat){
  u <- unique(e.hat)
  if (length(u)<K){# exist empty community
    sq <- 1:K
    m <- match(sq,u)
    c.empty <- sq[is.na(m)] # find the empty community
    t <- table(e.hat)
    c.max <- as.numeric(names(which.max(t))) # the biggest community
    c.max.id <- which(e.hat==c.max)
    i=1
    for (k in c.empty){# move a node from the biggest community to the empty community
      move.id <- c.max.id[i]
      e.hat[move.id] <- k
    }
  }
  return(e.hat)
}

### Metrics ###
com.dist <- function(label.true,label.hat){
  n <- length(label.true)
  flag <- 0
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      f1 <- label.true[i]==label.true[j]
      f2 <- label.hat[i]==label.hat[j]
      f3 <- f1+f2
      f4 <- as.numeric(f3==1)
      flag <- flag + f4
    }
  }
  dist <- 2/(n*(n+1))*flag
  return(dist)
}

# The complete log likelihood 
com.lik <- function(Adj,e,P,Pi){# decouples row and column labels
  n <- nrow(Adj)
  K <- length(Pi)
  K.I <- diag(K)
  e.OneHot <- as.matrix(K.I[e, ]) #n*K
  n.vec <- as.matrix(apply(e.OneHot,2,sum)) # 1*K, size of communities
  
  if (K==1){
    diag.minus <- n.vec
  }else{
    diag.minus <- diag(as.vector(n.vec))
  }
  n.mat <- n.vec %*% t(n.vec) # - diag.minus # number of possible edges
  m.mat <- t(e.OneHot)%*%Adj%*%e.OneHot #k*k
  
  #cl0.mat <- Adj.e.mat%*%t(log(P+1e-30))+(rep.nrows(n.vec,n)-Adj.e.mat)%*%t(log(1-P+1e-30)) # n*K
  #cl1 <- sum(sum(cl0.mat * e.OneHot))
  cl0.mat <- m.mat*(log(P+1e-30))+(n.mat-m.mat)*(log(1-P+1e-30)) # K*K
  cl1 <- sum(sum(cl0.mat))
  cl2 <- sum(n.vec*log(Pi))
  com.lik <- cl1+cl2
  return(com.lik)
}

