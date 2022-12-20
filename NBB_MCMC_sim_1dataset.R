##### This code is to show how to use the function to sample from the posterior of the parameters in the hierachical model #####
##### The simulation in the paper generates 20 sets of data #####
##### Here we just generate one set of data #####
##### This is intend to use as a tutorial of how to use the function we write for the sampler in the paper #####

rm(list = ls())

#sample size
N <- 100

#Some already known values are x1_i's and x_2i's

#For x1i, let's assume it's two dimensional, with x11 and x12. 
#For x2i, let's assume it's two dimensional, with x21 and x22.
set.seed(987987)

x11 <- rnorm(N,.5,.5)
x12 <- rnorm(N,.5,.5)

x21 <- rnorm(N,.5,.5)
x22 <- rnorm(N,.5,.5)

#build the design matrix
X1 <- matrix(c(rep(1,N),x11,x12),ncol=3)
X2 <- matrix(c(rep(1,N),x21,x22),ncol=3)

#Set the true parameter values
alpha_true <- c(1,-0.1,0.1)
beta_true <- c(0.5,-0.1,0.1)
delta_true <- 0.01
r <- 100

#Calculate p_i and psi_i from the model
logit_psi_i <- X1%*%alpha_true
psi_i <- 1/(1+exp(-logit_psi_i))

logit_pi <- X2%*%beta_true+logit_psi_i*delta_true
pi <- 1/(1+exp(-logit_pi))


#Simulate the observation m_i's and y_i's and store them
set.seed(978732987)

mi <- rnbinom(N,r,prob=1-psi_i)
yi <- rbinom(N,mi,pi)

#######################################################################################
# we have the data observations mi and yi
# and we have the covariates X1 and X2

# make sure to set the correct working directory
setwd("G:/My Drive/NIH/NegativeBinomial/NBB_code_simulation")
source("NBB_MCMC_func.R")

result <- NBB_MCMC(X1,X2,mi,yi)







  
