### Generate adjacency matrix for DCSBM ###
Adj.Generating.DCSBM <- function(N,Clusters_distribution,
                                 P,Theta,
                                 sp = T,del_0d = T){
  num_clusters <- length(Clusters_distribution)
  label_clusters <- 1:num_clusters
  clusters <- sample(label_clusters,size=N,replace=T,
                     prob=Clusters_distribution)
  
  K_I <- diag(num_clusters)
  C_Mat <- K_I[clusters, ]
  
  Theta_Mat <- diag(Theta)
  
  P_A <- Theta_Mat%*%(C_Mat%*%P%*%t(C_Mat))%*%Theta_Mat
  
  RM <- matrix(runif(N*N),N,N)
  RM[lower.tri(RM)] <- t(RM)[lower.tri(RM)] #Symmetrization
  
  Adj <- (RM < P_A)+0
  diag(Adj) <- 0
  
  if(del_0d){
    Ind_0 <- apply(Adj,1,sum)!=0
    Adj <- Adj[Ind_0,Ind_0]
    clusters <- clusters[Ind_0] 
  }
  
  if(sp){
    library(Matrix)
    Adj <- Matrix(Adj,sparse=T)
  }
  
  re <- list(clusters,Adj)
  names(re) <- c('clusters','Adj')
  return(re)
}

rep.nrows<-function(x,nrows){
  re <- matrix(rep(x,nrows),nrows,length(x),byrow = T)
  return(re)
}

### Initialization: SCORE+ ###
ini.SCORE <- function(Adj, K){

  # get top K eigenvectors
  eig.out = RSpectra::eigs(Adj, k = K)
  eig.vec.w = eig.out$vectors %*% diag(eig.out$values) # reweight eigenvectors by eigen values
  
  # get ratio matrix 
  ratios = eig.vec.w[,2:K] / eig.vec.w[,1]
  
  # k-means
  labels = kmeans(ratios, K, nstart = 100)$cluster
  
  return(labels)
}


### The proposed method ###
Est.DCSBM <- function(Adj,K,tun.c=10,sp=T,del_0d=T){
  # small tun.c prefer large K
  if(sp){
    library(Matrix)
  }
  if(del_0d){
  Ind_0 <- apply(Adj,1,sum)!=0
  Adj <- Adj[Ind_0,Ind_0]
  }
  ### Spectral clustering ###
  # the method of Spectral Clustering with Permutations (SCP)
  n <- nrow(Adj)
  #rho <- sum(Adj)/(n*(n-1))
  #pAdj<-Adj+as.numeric(0.025*rho)*matrix(1,nrow=n,ncol=n) #purmutated Adj
  d <- rowSums(Adj)
  
  # compute the Laplace Matrix
  #c  <-  0.1
  #delta  <-  c * max(d) # tunning parameter for graph laplacian
  #d.inv  <-  1 / sqrt( delta + d )
  #L.delta  <-  t(d.inv * pAdj) * d.inv # graph laplacian with ridge regularization
  
  ### Initialization ###
  ########################## SCORE #########################
  e.hat <- ini.SCORE(Adj, K)
  
  theta.hat <- n*d/sum(d)
  
  ### Outer iterations ###
  outer.iter.max = 25
  inner.iter.max = 25
  converge.thresh = 0.01
  
  Iter <- numeric(outer.iter.max+1)
  Iter[1] <- inner.iter.max
  break.outer.count <- 0
  
  for(i in 1:outer.iter.max){
    
    ### Inner iterations(EM) ###
    inner.count <- 0 #the counter of the iterations times
    break.flag <- FALSE #the sign of convenience
    
    K.I <- diag(K)
    e.OneHot <- as.matrix(K.I[e.hat, ]) #n*K
    n.vec <- as.matrix(apply(e.OneHot,2,sum)) # 1*K, size of communities
    Pi.hat <- n.vec/sum(n.vec)
    
    # d.vec <- t(e.OneHot)%*%d
    # theta.hat <- d*e.OneHot%*%(n.vec/d.vec)

    theta.e <- t(theta.hat)%*%e.OneHot

    n.mat <- t(theta.e)%*%theta.e # number of possible edges
    m.mat <- t(e.OneHot)%*%Adj%*%e.OneHot #k*k
    P.hat <- m.mat / n.mat
    # Here we accomplish initialization, we have e.hat,Pi.hat,P.hat
    
    Adj.e.mat <- Adj%*%e.OneHot #n*K
    
    while((!break.flag)&&(inner.count<inner.iter.max)){
      inner.count <- inner.count+1
      ### E-step ###
      log.tau0.1 <- - theta.hat%*%theta.e%*%t(P.hat)
      log.tau0.2 <- t(rep.nrows(log(theta.hat)*d,K))
      log.tau0.3 <- Adj%*%log(theta.hat*e.OneHot%*%t(P.hat)+1e-30)
      log.prod <- log.tau0.1+log.tau0.2+log.tau0.3
      
      log.tau0.mat <- log.prod+rep.nrows(log(Pi.hat),n) # log numerator of tau
      log.tau1.mat <- t(apply(log.tau0.mat,1,function(x) x-max(x)))
      tau0.mat <- exp(log.tau1.mat)
      tau <- t(apply(tau0.mat,1,function(x) x/(sum(x)+1e-30)))
      
      tau.n <- colSums(tau)
      ### M-step ###
      #  Calculate Pi.hat with selection of tuning parameter
      multi <- seq(from=0.1, to=4, by=0.2)
      gamma0 <- sqrt(1/n)
      gamma.seq <- multi*gamma0

      Pi.update <- t(rep.nrows(1/(1-K*gamma.seq),K))*(rep.nrows(tau.n/n,length(multi))-gamma.seq) # length(multi) * K; Pi.update <- 1/(1-K*gamma)*(tau.n/n-gamma)
      Pi.update[Pi.update<0] <- 0 # nrows = length(multi)
      Pi.mat <- t(apply(Pi.update,1,function(x) x/sum(x)))
      
      bic1.mat <- exp(log.prod)
      bic2.vec <- log(bic1.mat %*% t(Pi.mat)) # n*length(multi)
      bic.term1 <- -2*colSums(bic2.vec)
      
      K.pi <- rowSums(Pi.mat != 0)
      # tun.c <- 10
      bic.term2 <- tun.c*n*log(K.pi)+K.pi*(K.pi+1)/2*log(n)
      
      cbic.gamma <- bic.term1+bic.term2
      
      #plot(bic.term1)
      #plot(cbic.gamma,col='blue')
      gamma.id <- which.min(cbic.gamma)
      Pi.hat <- Pi.mat[gamma.id,]
      print(Pi.hat)
      
      # Calculate P.hat (which not dependent on gamma)
      P.num <- t(tau)%*%Adj.e.mat
      P.den <- t(tau)%*%theta.hat%*%theta.e
      P.update <- P.num/(P.den+1e-30)
      P.update[P.update>1] <- 1
      P.hat <- P.update
      
      # Calculate theta.hat
      g <- tau %*% P.hat %*% t(e.OneHot)# n*n
      for (i in 1:n){
        theta.temp <- theta.hat
        theta.temp[i] <- 0
        hi <- g[i,] %*% theta.temp
        theta.hat[i]<- (- hi + sqrt(hi^2+8*d[i]*g[i,i]))/(4*g[i,i])
      }
      theta.hat <- n/sum(theta.hat)*theta.hat
      
              Q.update <- bic.term1[gamma.id]
        if (inner.count > 1){
          #compute the new relative changes in likelihood
          delta <- (Q.update - Q.old)/(Q.old+1e-30)
          #break if converged
          break.flag <- delta < converge.thresh
        }
        Q.old <- Q.update
        
      if(inner.count<2){
        # for a better estimation of theta, there is no rush to break the inner loop
        # compute the new value of likelihood
        Q.update <- bic.term1[gamma.id]
        Q.old <- Q.update
      }else{
        # break if there exists 0 in pi
        break.flag <- sum(Pi.hat==0)>0
        if(break.flag){
          # go to next outer loop
          break.outer.count <- -1
          # remove the parameter related to the removed community
          remove.id <- which(Pi.hat==0)
          tau <- tau[,-remove.id]
          Pi.hat <- Pi.hat[-remove.id]
          P.hat <- P.hat[-remove.id,-remove.id]
        }else{
          # if K converges
          # compute the new relative changes in likelihood
          delta <- (Q.update - Q.old)/(Q.old+1e-30)
          # break if converged
          break.flag <- delta < converge.thresh
        }
      }
    }
    
    # update number of community K
    K <- length(Pi.hat)

    # update column label e
    e.update <- theta.hat%*%t(theta.hat)%*%tau%*%P.hat + t(Adj)%*%tau%*%log(P.hat+1e-30)
    e.hat <- apply(e.update,1,which.max)
    e.hat <- fix.e(K,e.hat)
    # record the number of iterations
    Iter[i+1] <- inner.count
    
    # break rule
    if (inner.count<=3) {
      break.outer.count <- break.outer.count+1
    }
    if (break.outer.count==2){break}
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