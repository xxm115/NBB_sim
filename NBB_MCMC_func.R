
library(BayesLogit) #Package by Polson to generate from Polya Gamma dist.
library(MASS) #multivariate normal distribution
library(mvtnorm)

NBB_MCMC=function(X1, X2, mi, yi, samp_size=10000, alpha_prior_mean=rep(0,ncol(X1)), alpha_prior_sigma=10000*diag(ncol(X1)),
                  beta_prior_mean=rep(0,ncol(X2)+1),beta_prior_sigma=10000*diag(ncol(X2)+1),
                  alpha_init= rep(0,ncol(X1)), beta_init=rep(0,ncol(X2)), delta_init=0){
  
  ptm <- proc.time()
  
  #prior
  a <- alpha_prior_mean
  A <- alpha_prior_sigma
  b <- beta_prior_mean
  B <- beta_prior_sigma
  
  #data
  mi_rm0 <- mi[mi!=0]
  yi_rm0 <- yi[mi!=0]
  
  #Initial value
  alpha_samp <- alpha_init
  beta_samp <- c(beta_init, delta_init)
  
  #For efficient computation
  solveB <- solve(B)
  Bsolveb <- solve(B)%*%b
  solveA <- solve(A)
  Asolvea <- solveA%*%a
  
  index <- mi!=0
  mino0 <- sum(index)
  p_beta <- ncol(X2)
  X1_index <- X1[index,]
  X2_index <- X2[index,]
  
  #Recording
  record_alpha <- c()
  record_beta <- c()
  record_delta <- c()
  
  for (j in 1:samp_size){
    
    #update X2_full
    
    logit_psi_samp <- X1_index%*%alpha_samp
    X2_full <- cbind(X2_index,logit_psi_samp)
    
    #sample omega2
    
    omega2_samp <- c()
    omega2_samp <- rpg(mino0,as.numeric(mi_rm0),X2_full%*%beta_samp)
    Omega2_mat <- diag(omega2_samp)
    
    #sample beta
    X2Omg2 <- t(X2_full)%*%Omega2_mat
    
    Z2_omg <- (yi_rm0-0.5*mi_rm0)/omega2_samp
    V2_omg <- solve(X2Omg2%*%X2_full+solveB)
    m2_omg <- V2_omg%*%(X2Omg2%*%Z2_omg+Bsolveb)
    
    beta_samp <- mvrnorm(1,m2_omg,V2_omg)
    delta_samp <- beta_samp[p_beta+1]
    
    #sample omega1
    
    omega1_samp <- c()
    omega1_samp <- rpg(mino0,as.numeric(r+mi_rm0),X1_index%*%alpha_samp)
    Omega1_mat <- diag(omega1_samp+omega2_samp*delta_samp^2)
    
    #sample alpha
    X1Omg1 <- t(X1_index)%*%Omega1_mat
    
    Z1_omg <- ((mi_rm0-r)/2+(yi_rm0-mi_rm0/2)*delta_samp-omega2_samp*X2_index%*%beta_samp[1:p_beta]*delta_samp)/(omega1_samp+omega2_samp*delta_samp^2)
    V1_omg <- solve(X1Omg1%*%X1_index+solveA)
    m1_omg <- V1_omg%*%(X1Omg1%*%Z1_omg+Asolvea)
    
    alpha_samp <- mvrnorm(1,m1_omg,V1_omg)
    
    
    #record the sampling
    record_beta <- cbind(record_beta,beta_samp[1:p_beta])
    record_delta <- cbind(record_delta,delta_samp)
    record_alpha <- cbind(record_alpha,alpha_samp)
  }
  
  result_list <- list(
    "record_alpha" = as.matrix(record_alpha),
    "record_beta" = as.matrix(record_beta),
    "record_delta" = as.matrix(record_delta)
  )
  
  
  return(list("result_list" = result_list , "run_time" = proc.time() - ptm))
}
