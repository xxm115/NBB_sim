data {
  int<lower=0> N; //number of samples
  int<lower=0> C1; //number of columns for X1
  int<lower=0> C2; //number of columns for X2
  int<lower=0> Y[N]; //number of mutants
  int<lower=0> M[N]; //number of mtDNA
  real r;//r
  matrix[N,C1] X1;//"design" matrix X1 (covariate for log(lambda))
  matrix[N,C2] X2;//Design matrix (covariate for logit(pi))
  vector[C1] mu_prior_alpha;
  matrix[C1,C1] sigma_prior_alpha;
  vector[C1] mu_prior_beta;
  matrix[C1,C1] sigma_prior_beta;
}
parameters {
  //regression coefficients:
  vector[C1] alpha;
  vector[C2] beta;
  real delta;
}
model {
    // Priors on regression coeff:
  alpha ~ multi_normal(mu_prior_alpha,sigma_prior_alpha);
  beta ~ multi_normal(mu_prior_beta,sigma_prior_beta);
  delta ~ normal(0.0,1.0E3);
  
  
  // Likelihood: Binomial with a logit link
  M ~ neg_binomial_2(r*exp(X1 * alpha),r);
  Y ~ binomial_logit(M, X2 * beta+ delta * X1 * alpha);
}
