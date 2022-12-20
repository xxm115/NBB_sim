##### This is to simulate data from negative binomial and binomial hierarchical model #####
##### and use STAN and Polson's MCMC method to sample from posterior distribution for parameters #####
##### and compare the result and efficiency of computation for the two methods  #####

##### This is the simulation shown in the paper #####

## First step: Simulate Some Data ##

## The model is m_i~NB(r,psi_i), and logit(psi_i)=x1_i'*alpha
##              y_i~Binomial(m_i,p_i),logit(p_i)=x2_i'*beta+eta_i*delta 

## We generate 20 sets of data (m_i's and y_i's) for set X1 X2 and covariates, and true parameter values (beta, alpha, delta)

rm(list = ls())

#sample size
N <- 100

#Some already known values are x1_i's and x_2i's

#For x1i, let's assume it's two dimensional, with x11 and x12. 
#For x2i, let's assume it's two dimensional, with x21 and x22.
set.seed(987987)

x11 <- rnorm(N,1,1)
x12 <- rnorm(N,1,1)

x21 <- rnorm(N,1,2)
x22 <- rnorm(N,1,2)

#build the design matrix
X1 <- matrix(c(rep(1,N),x11,x12),ncol=3)
X2 <- matrix(c(rep(1,N),x21,x22),ncol=3)

#Set the true parameter values
alpha_true <- c(1,-0.2,0.7)
beta_true <- c(0.5,-0.2,0.7)
delta_true <- 0.01
r <- 100

#Calculate p_i and psi_i from the model
logit_psi_i <- X1%*%alpha_true
psi_i <- 1/(1+exp(-logit_psi_i))

logit_pi <- X2%*%beta_true+logit_psi_i*delta_true
pi <- 1/(1+exp(-logit_pi))

#Generate nsample=20 sets of data
nsample <- 20

#store the data in datacombine list
datacombine <- list(NULL)
length(datacombine) <- nsample

#Simulate the observation m_i's and y_i's and store them
set.seed(978732987)

for (m in 1:nsample){
  mi <- rnbinom(N,r,prob=1-psi_i)
  yi <- rbinom(N,mi,pi)
  datacombine[[m]][[1]] <- mi
  datacombine[[m]][[2]] <- yi
  
}




########################################
# Use Stan to sample from posterior distribution of the parameters
library(rstan)
rstan_options(auto_write = TRUE)


##Now start with stan to estimate parameters
# make sure to set the correct working directory
setwd("G:/My Drive/NIH/NegativeBinomial/NBB_code_simulation")

Fit_NBB_result <- list(NULL)
length(Fit_NBB_result) <- nsample
result <- list(NULL)
length(result) <- nsample

tm <- c()

for (m in 1:nsample){

  mi <- datacombine[[m]][[1]]
  yi <- datacombine[[m]][[2]]
  
  
  #List of data we are going to use
  NBB_data <- list(
    N = length(mi), 
    Y = yi,
    M = mi,
    X1 = X1,
    X2 = X2,
    r=r,
    C1 = ncol(X1),
    C2 = ncol(X2),
    mu_prior_alpha = as.vector(rep(0,ncol(X1))),
    sigma_prior_alpha = 10000*diag(ncol(X1)),
    mu_prior_beta = as.vector(rep(0,ncol(X2))),
    sigma_prior_beta = 10000*diag(ncol(X2))
  )
  
  ptm <- proc.time()
  
  NBB_fit <- stan(file="NB_B_stan_sim.stan",
                  data=NBB_data, iter=10000, chains=1)
  
  run_time <- proc.time() - ptm
  
  ext <- rstan::extract(NBB_fit)
  
  result_list <- list(
    "record_alphas" = as.matrix(ext$alpha),
    "record_betas" = as.matrix(ext$beta),
    "record_delta" = ext$delta
  )
  
  result[[m]] <- result_list
  
  Fit_NBB_result[[m]] <- NBB_fit
  
  
  tm <- rbind(tm,run_time[3])
  
}

save.image(file="NBB_sim_STAN_1.RData")


#########################################
# Use the MCMC algorithm in the paper to sample from the posterior distribution of the parameters

result <- list(NULL)
length(result) <- nsample

# make sure to set the correct working directory
setwd("G:/My Drive/NIH/NegativeBinomial/NBB_code_simulation")
source("NBB_MCMC_func.R")

tm <- c()
for (m in 1:nsample){
  #data
  mi <- datacombine[[m]][[1]]
  yi <- datacombine[[m]][[2]]
  mi_rm0 <- mi[mi!=0]
  yi_rm0 <- yi[mi!=0]
  
  #set sampling size
  samp_size <- 10000
  
  MCMC_result <- NBB_MCMC(X1,X2,mi,yi,samp_size=samp_size)

  result[[m]] <- MCMC_result$result_list
  tm <- rbind(tm,MCMC_result$run_time[3])
}


save.image(file="NBB_sim_MCMC_1.RData")

#############################################



####################################################
# Figures #

load("NBB_sim_Stan_1.RData")

####################################################
###### Make the plot for Stan Sampling ###########
####################################################
par(mfrow = c(1, ncol(X1))) #for alpha
#############
for (j in 1:ncol(X1)){
  
  plot(density(result[[1]]$record_alphas[,j]),main=paste0("alpha_", (j-1)),ylab="Density",col = rgb(0, 0, 255, 50, maxColorValue=255))
  
  for (i in 1:nsample){
    lines(density(result[[i]]$record_alphas[,j]),col = rgb(0, 0, 255, 50, maxColorValue=255))
  }
  abline(v=alpha_true[j],col="red")
  
}

##########################
par(mfrow = c(1, ncol(X2)+1)) #for beta and delta
#############

for (j in 1:ncol(X2)){
  
  plot(density(result[[1]]$record_betas[,j]),main=paste0("beta_", (j-1)),ylab="Density",col = rgb(0, 0, 255, 50, maxColorValue=255))
  
  for (i in 1:nsample){
    lines(density(result[[i]]$record_betas[,j]),main="beta 1",ylab="Density",col = rgb(0, 0, 255, 50, maxColorValue=255))
  }
  abline(v=beta_true[j],col="red")
  
}

plot(density(result[[1]]$record_delta),main="delta",ylab="Density",col = rgb(0, 0, 255, 50, maxColorValue=255))

for (i in 1:nsample){
  lines(density(result[[i]]$record_delta),col = rgb(0, 0, 255, 50, maxColorValue=255))
}
abline(v=delta_true,col="red")
####################################################


load("NBB_sim_MCMC_1.RData")
####################################################
###### Make the plot for MCMC Sampling ###########
####################################################
par(mfrow = c(1, ncol(X1))) #for alpha
#############

for (j in 1:ncol(X1)){
  
  plot(density(result[[1]]$record_alpha[j,(samp_size/2):samp_size]),main=paste0("alpha_",j-1),ylab="Density",col = rgb(0, 0, 255, 50, maxColorValue=255))
  
  for (i in 1:nsample){
    lines(density(result[[i]]$record_alpha[j,(samp_size/2):samp_size]),col = rgb(0, 0, 255, 50, maxColorValue=255))
  }
  abline(v=alpha_true[j],col="red")
}

##########################
par(mfrow = c(1, ncol(X2)+1)) #for beta and delta
#############

for (j in 1:ncol(X2)){
  plot(density(result[[1]]$record_beta[j,(samp_size/2):samp_size]),main=paste0("beta_",j-1),ylab="Density",col = rgb(0, 0, 255, 50, maxColorValue=255))
  
  for (i in 1:nsample){
    lines(density(result[[i]][[2]][j,(samp_size/2):samp_size]),col = rgb(0, 0, 255, 50, maxColorValue=255))
  }
  abline(v=beta_true[j],col="red")
}

plot(density(result[[1]]$record_delta[(samp_size/2):samp_size]),main="delta",ylab="Density",col = rgb(0, 0, 255, 50, maxColorValue=255))

for (i in 1:nsample){
  lines(density(result[[i]]$record_delta[(samp_size/2):samp_size]),main="delta",ylab="Density",col = rgb(0, 0, 255, 50, maxColorValue=255))
}
abline(v=delta_true,col="red")

#######################################
#######################################

# traceplot for MCMC sampling #
par(mfrow = c(ncol(X1), 1)) #for alpha
#############

for (j in 1:ncol(X1)){
  
  plot(result[[1]]$record_alpha[j,(samp_size/2):samp_size],main=paste0("alpha_",j-1),ylab="Density",col = rgb(0, 0, 255, 50, maxColorValue=255),type='l')
  
  for (i in 1:nsample){
    
    lines(result[[i]]$record_alpha[j,(samp_size/2):samp_size],col = rgb(0, 0, 255, 50, maxColorValue=255),type='l')
    abline(h=alpha_true[j],col="red")
    
  }
  
}

#############
par(mfrow = c(ncol(X2)+1, 1)) #for beta and delta

#############

for (j in 1:ncol(X2)){
  
  plot(result[[1]]$record_beta[j,(samp_size/2):samp_size],main=paste0("beta_",j-1),ylab="Density",col = rgb(0, 0, 255, 50, maxColorValue=255),type='l')
  
  for (i in 1:nsample){
    
    ####First use non-informative prior for x_a S_a and x_z, S_z), without model discrepancy
    lines(result[[i]]$record_beta[j,(samp_size/2):samp_size],col = rgb(0, 0, 255, 50, maxColorValue=255),type='l')
    abline(h=beta_true[j],col="red")
    
  }
}

plot(result[[1]]$record_delta[(samp_size/2):samp_size],main="delta",ylab="Density",col = rgb(0, 0, 255, 50, maxColorValue=255),type='l')

for (i in 1:nsample){
  
  ####First use non-informative prior for x_a S_a and x_z, S_z), without model discrepancy
  lines(result[[i]]$record_delta[(samp_size/2):samp_size],col = rgb(0, 0, 255, 50, maxColorValue=255),type='l')
  abline(h=delta_true,col="red")
  
} 

#################################################################################################
#################################################################################################


#################################ESS vs ESR######################################################
rm(list = ls())

setwd("G:/My Drive/NIH/NegativeBinomial/NBB_code_simulation")

library(coda)  ## A convenient R package to analyse MCMCs

#stan
load("NBB_sim_Stan_1.RData")


library(rstan)

ESS_Stan <- c()

for (i in 1:nsample){
  ext <- rstan::extract(Fit_NBB_result[[i]],permuted=FALSE,pars=c("alpha","beta","delta"))
  dimens <- dim(ext)
  
  Use_MCMC <- list()
  
  n_chain <- dimens[2]
  for (j in 1:n_chain){
    Use_MCMC[[j]] <- cbind(ext[,j,1],ext[,j,2],ext[,j,3],ext[,j,4],ext[,j,5],ext[,j,6],ext[,j,7])
    colnames(Use_MCMC[[j]]) <- c("alpha0","alpha1","alpha2","beta0","beta1","beta2","delta")
  }
  
  
  Use_MCMC_List <- lapply(Use_MCMC,mcmc)
  
  ESS_Stan <- rbind(ESS_Stan,effectiveSize(Use_MCMC_List))
}
ESR_Stan <- ESS_Stan/c(tm)
tm


load("NBB_sim_MCMC_1.RData")

Burnin <- 5000
ESS_Polson <- c()

for (i in 1:nsample){
  
  Use_MCMC <- list()
  alpha_list <- list(t(result[[i]][[1]]))
  beta_list <- list(t(result[[i]][[2]]))
  delta_list <- list(t(result[[i]][[3]]))
  
  Use_MCMC[[1]] <- cbind(alpha_list[[1]][-(1:Burnin),],beta_list[[1]][-(1:Burnin),],delta_list[[1]][-(1:Burnin)])
  #since we only have one chain here
  
  colnames(Use_MCMC[[1]]) <- c("alpha0","alpha1","alpha2","beta0","beta1","beta2","delta")
  
  
  Use_MCMC_List <- lapply(Use_MCMC,mcmc)
  #traceplot(Use_MCMC_List)
  
  ## effective sample size
  ESS_Polson <- rbind(ESS_Polson,effectiveSize(Use_MCMC_List))
  
  #summary(mcmc.list(Use_MCMC_List))
  
  
  
}
ESR_Polson <- ESS_Polson/c(tm)
tm

library(ggplot2)

par <- c()

for (i in 1:7){
  par <- c( par,ESS_Stan[,i])
  par <- c( par,ESS_Polson[,i])
}

ESS_rebuild <- data.frame(par_val=par, method=rep(rep(c("Stan","Polson"),each=nsample),(ncol(X1)+ncol(X2)+1)),
                          par_name=rep(c("alpha0","alpha1","alpha2","beta0","beta1","beta2","delta"),each=2*nsample))


ggplot(ESS_rebuild,aes(x=par_name,y=par_val,fill=method))+
  geom_boxplot(position=position_dodge(.6),width=.5)+
  xlab("")+ylab("ESS")+ggtitle("Boxplot of ESS for parameters of Stan and Polson's method (NB-B Model)")

######################################################


#############################################################################################################
##############################################################################################################
library(ggplot2)

par <- c()

for (i in 1:7){
  par <- c( par,ESR_Stan[,i])
  par <- c( par,ESR_Polson[,i])
}

ESR_rebuild <- data.frame(par_val=par, method=rep(rep(c("Stan","Polson"),each=nsample),(ncol(X1)+ncol(X2)+1)),
                          par_name=rep(c("alpha0","alpha1","alpha2","beta0","beta1","beta2","delta"),each=2*nsample))


ggplot(ESR_rebuild,aes(x=par_name,y=par_val,fill=method))+
  geom_boxplot(position=position_dodge(.6),width=.5)+
  xlab("")+ylab("ESR")+ggtitle("Boxplot of ESR for parameters of Stan and Polson's method (NB-B Model)")

######################################################








