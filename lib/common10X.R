
# common paramters
modalCNGuess <- 2

# vector functions
collapseVector <- function(v, n) { # sum every n adjacent elements of a vector 
  unname(tapply(v, (seq_along(v)-1) %/% n, sum))
}

# correlation-based distance of a set of already-centered Z scores
pearson.dist <- function (m) { # m is a matrix
  m <- m / sqrt (rowSums (m^2))
  m <-  tcrossprod (m)
  m <- as.dist(m)
  0.5 - m / 2
}

# Hidden Markov Model based on Poisson emissions
HMM_set_persistence <- function(persistence){
    persistence <<- persistence 
    tp <<- matrix(log((1-persistence)/(3-1)), 3, 3)
    for(i in 1:3) tp[i,i] <<- log(persistence)    
}
HMM_viterbi <- function(cellI){
  
  #message(cellI)

  # 0. emission probabilities
  maxCount <- round(2 * zLayers[,cellI,zLN$expG], 0)
  ep <- sapply(c(zLN$expL, zLN$exp0, zLN$expG), function(exp){
    log(dpois(pmin(maxCount, round(zLayers[,cellI,zLN$raw],0)), zLayers[,cellI,exp]))
  })

  # 1. initialization (observation t=1)
  T         <- nrow(ep) # length of the sequence of observations
  N         <- ncol(ep)
  delta     <- log(matrix(0, nrow=T, ncol=N))
  delta[1,] <- sapply(1:N, function(i) log(1/N) + ep[1,i])  
  phi       <-     matrix(NA, nrow=T, ncol=N)  
  
  # 2. recursion;
  # NB: these 'for' loops are faster than apply methods with array as implemented and given recursion restrictions
  for (t in 2:T){
    pt <- t - 1
    for (j in 1:N){   # j = this hs
      ep_j <- ep[t,j] 
      for (i in 1:N){ # i = prev hs
        delta_ <- delta[pt,i] + tp[i,j] + ep_j      
        if(delta[t,j] < delta_){
          delta[t,j] <- delta_
          phi[pt,j]  <- i
        }
      }
    }
  }   

  # 3. termination
  prob <- -Inf
  hsi  <- rep(-1, T)
  for (j in 1:N){
    if(prob < delta[T,j]){
      prob <- delta[T,j]
      hsi[T] <- j
    }
  }  
  
  # 4. reconstruction
  for (t in (T-1):1) hsi[t] <- phi[t,hsi[t+1]]  
  
  # return
  hsi - 2
  # list(prob=prob, hsi=hsi) # overall likelihood and ordered set of most likely states  
}

